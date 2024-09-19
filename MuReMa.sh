#!/bin/bash

# Function to display script usage
function display_usage() {
    echo "Usage: $0 -s FILE -d FILE [-r INTEGER] [-t INTEGER] [-h]"
    echo ""
    echo "Options:"
    echo "  -s FILE     Samples file OR path to the samples file"
    echo "  -d FILE     Multifasta file OR path to the multifasta file"
    echo "  -r INTEGER  Read length (default: 150)"
    echo "  -t INTEGER  Average Depth threshold (default: 1000)"
    echo "  -h          Display this help message"
    exit 0
}

# Initialize variables with default values
samples_file=""
multifasta_file=""
read_length=150
avg_depth_threshold=1000

# Parse options using getopts
while getopts "s:d:r:t:h" option; do
    case "$option" in
        s) samples_file=$PWD/$OPTARG;;
        d) multifasta_file=$PWD/$OPTARG;;
        r) read_length=$OPTARG;;
        t) avg_depth_threshold=$OPTARG;;
        h) display_usage;;
        :) printf "missing argument for -%s\n" "$OPTARG" >&2; display_usage >&2; exit 1;;
        \?) printf "illegal option: -%s\n" "$OPTARG" >&2; display_usage >&2; exit 1;;
    esac
done

# Check if required arguments are provided
if [ -z "$samples_file" ] || [ -z "$multifasta_file" ]; then
    echo "Error: At least arguments -s and -d must be provided."
    display_usage
fi

# Log file
log_file="${samples_file}.log"

# Log the initial prompt to the log file
{
    echo "Script invoked with the following arguments:"
    echo "Samples file: $samples_file"
    echo "Multifasta DB: $multifasta_file"
    echo "Read length: $read_length"
    echo "Average Depth threshold: $avg_depth_threshold"
} > "$log_file"

# Inputs validation
if [ ! -e "$samples_file" ]; then
    echo "Error: $samples_file is not a valid file." | tee -a "$log_file"
    exit 1
fi

if [ ! -e "$multifasta_file" ]; then
    echo "Error: $multifasta_file is not a valid file." | tee -a "$log_file"
    exit 1
fi

if ! [[ "$read_length" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $read_length is not a valid integer greater than 0." | tee -a "$log_file"
    exit 1
fi

if ! [[ "$avg_depth_threshold" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $avg_depth_threshold is not a valid integer greater than 0." | tee -a "$log_file"
    exit 1
fi

# Ensure required tools are installed
{
    command -v bowtie2 >/dev/null 2>&1 || { echo >&2 "bowtie2 is required but it's not installed. Aborting."; exit 1; }
    command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools is required but it's not installed. Aborting."; exit 1; }
    command -v python3 >/dev/null 2>&1 || { echo >&2 "python3 is required but it's not installed. Aborting."; exit 1; }
    command -v ivar >/dev/null 2>&1 || { echo >&2 "ivar is required but it's not installed. Aborting."; exit 1; }
} | tee -a "$log_file"

# Make the Python scripts executable (if not already)
chmod +x formater.py
chmod +x grapher.py

# Create DB directory and index the multifasta
mkdir -p DB_dir
cd DB_dir
cp "$multifasta_file" murema_DB.fasta
if [ ! -f "murema_DB_index.1.bt2" ]; then
    bowtie2-build -f murema_DB.fasta murema_DB_index
fi
cd ../

# Process each sample
while IFS= read -r sample_name; do
    mkdir -p "$sample_name"
    cd "$sample_name"
    {    
        echo ""
        echo "Processing Sample: $sample_name"
    } | tee -a "$log_file"

    # Align reads and sort bam
    bowtie2 --end-to-end --very-sensitive -x "../DB_dir/murema_DB_index" -1 "../${sample_name}_1.fastq.gz" -2 "../${sample_name}_2.fastq.gz" | samtools sort | samtools view -@ 8 -b -F 4 -q 1 -o "${sample_name}.sorted.bam"
    if [ $? -ne 0 ]; then
        {
            echo ""
            echo "Error: bowtie2 or samtools sorting failed for sample $sample_name"
        } | tee -a "$log_file"
        cd ../
        continue
    fi

    # Index sorted BAM file and generate a summary TSV file for the sample
    samtools index "${sample_name}.sorted.bam"
    samtools idxstats "${sample_name}.sorted.bam" > "${sample_name}.tsv"

    # Run formater.py script to filter and format the TSV files
    python3 ../formater.py "${sample_name}.tsv" "$sample_name" "$read_length" "$avg_depth_threshold"
    if [ $? -ne 0 ]; then
        {
            echo ""
            echo "Error: formater.py failed for sample $sample_name"
        } | tee -a "$log_file"
        cd ../
        continue
    fi
    sed -i 's/\r$//' "${sample_name}.filtered.tsv"
    sed -i 's/\r$//' "${sample_name}.refs.tsv"
    {
        echo ""
        echo "Finished creating filtered and refs reports for $sample_name"
    } | tee -a "$log_file"
    cd ../
done < "$samples_file"

# Prepare refs to create the consensus and graphs
cd DB_dir/
cat ../*/*.refs.tsv > all_refs.tsv
uniq all_refs.tsv > uniq_refs.tsv
sed -i 's/\r$//' "all_refs.tsv"
sed -i 's/\r$//' "uniq_refs.tsv"

# Read the references file and extract the sequences
while IFS= read -r ref_name; do
    awk -v RS='>' -v ref="$ref_name" '
    $1 == ref {
        print ">"$0;
        exit;
    }' murema_DB.fasta > "${ref_name}.fasta"
done < uniq_refs.tsv

# Index founded references
for f in *.fasta; do
    if [ ! -f "${f%.fasta}_index.1.bt2" ]; then
        bowtie2-build -f "$f" "${f%.fasta}_index"
    fi
done
cd ../

# Consensus & Grapher
while IFS= read -r sample_name; do
    cd "$sample_name"
    while IFS= read -r ref_name; do
        bowtie2 --end-to-end --very-sensitive -x "../DB_dir/${ref_name}_index" -1 "../${sample_name}_1.fastq.gz" -2 "../${sample_name}_2.fastq.gz" | samtools sort | samtools view -@ 8 -b -F 4 -q 1 -o "${sample_name}.${ref_name}.sorted.bam"
        samtools mpileup -A -d 6000000 -B -Q 0 -q 20 --reference ../DB_dir/${ref_name}.fasta $sample_name.$ref_name.sorted.bam | ivar consensus -p $sample_name.$ref_name.consensus -t 0.75
        {
            echo ""
            echo "Finished creating consensus sequence for $sample_name using $ref_name as reference"
        } | tee -a "$log_file"
        # Run grapher.py script to generate graphs per refs
        python3 ../grapher.py "${sample_name}.${ref_name}.sorted.bam" "$sample_name" "$ref_name" "$avg_depth_threshold"
        if [ $? -ne 0 ]; then
            {
                echo ""
                echo "Error: grapher.py failed for sample $sample_name"
            } | tee -a "$log_file"
            cd ../
            continue
        fi
    done < "${sample_name}.refs.tsv"
    {
        echo ""
        echo "Finished creating graphs for $sample_name"
        echo ""
        echo "Finished Processing Sample: $sample_name"
    } | tee -a "$log_file"
    cd ../
done < "$samples_file"

{
    echo ""
    echo "Script completed successfully." | tee -a "$log_file"
} | tee -a "$log_file"
