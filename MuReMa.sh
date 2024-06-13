#!/bin/bash

# Function to display script usage
function display_usage() {
    echo "Usage: $0 -s FILE -d FILE [-r INTEGER] [-m INTEGER] [-h]"
    echo ""
    echo "Options:"
    echo "  -s FILE     Path to the samples file"
    echo "  -d FILE     Path to the multifasta file"
    echo "  -r INTEGER  Read length (default: 150)"
    echo "  -m INTEGER  Mean Depth threshold (default: 1000)"
    echo "  -h          Display this help message"
    exit 0
}

# Initialize variables with default values
samples_file=""
multifasta_file=""
read_length=150
mean_depth_threshold=1000

# Parse options using getopts
while getopts "s:d:r:m:h" option; do
    case "$option" in
        s) samples_file=$OPTARG;;
        d) multifasta_file=$OPTARG;;
        r) read_length=$OPTARG;;
        m) mean_depth_threshold=$OPTARG;;
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
    echo "Mean Depth threshold: $mean_depth_threshold"
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

if ! [[ "$mean_depth_threshold" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $mean_depth_threshold is not a valid integer greater than 0." | tee -a "$log_file"
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
    {    
        echo ""
        echo "Processing Sample: $sample_name"
    } | tee -a "$log_file"
    mkdir -p "$sample_name"
    cd "$sample_name"

    # Align reads with bowtie2, sort, and filter with samtools
    bowtie2 --end-to-end --very-sensitive -x "../DB_dir/murema_DB_index" -1 "../${sample_name}_1.fastq.gz" -2 "../${sample_name}_2.fastq.gz" | samtools sort | samtools view -@ 8 -b -F 4 -q 1 -o "${sample_name}.sorted.bam"
    if [ $? -ne 0 ]; then
        {
            echo ""
            echo "Error: Bowtie2 or samtools sorting failed for sample $sample_name"
        } | tee -a "$log_file"
        cd ../
        continue
    fi

    # Index sorted BAM file and generate a summary TSV file for the sample
    samtools index "${sample_name}.sorted.bam"
    samtools idxstats "${sample_name}.sorted.bam" > "${sample_name}.tsv"

    ### Formater
    # Run formater.py script to filter and format the TSV files
    python3 ../formater.py "${sample_name}.tsv" "$sample_name" "$read_length" "$mean_depth_threshold"
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

    # Filter the BAM file by the different references and output to a separate BAM file for each one
    while IFS= read -r ref_name; do
        samtools view -b "${sample_name}.sorted.bam" "$ref_name" > "${sample_name}.${ref_name}.sorted.bam"
    done < "${sample_name}.refs.tsv"

    ### Grapher
    # Generate a dispersion graph for each reference
    while IFS= read -r ref_name; do
        python3 ../grapher.py "${sample_name}.${ref_name}.sorted.bam" "$sample_name" "$ref_name" "$mean_depth_threshold"
    done < "${sample_name}.refs.tsv"
    cd ../
    {
        echo ""
        echo "Finished creating filtered.tsv report and graphs for $sample_name"
    } | tee -a "$log_file"
done < "$samples_file"

### Consensus
cd DB_dir/
cat ../*/*.refs.tsv > all_refs.tsv
uniq all_refs.tsv > uniq_refs.tsv # unique references with an overall vertical depth >= $mean_depth_threshold 

# Read the references file and extract the sequences
while IFS= read -r ref_name; do
    # Extract the sequence for the reference name and save it to a separate fasta file
    awk -v RS='>' -v ref="$ref_name" '$1 == ref {print ">"$0}' murema_DB.fasta > "${ref_name}.fasta"
done < uniq_refs.tsv

# Index founded references
for f in *.fasta; do
    if [ ! -f "${f%.fasta}_index.1.bt2" ]; then
        bowtie2-build -f "$f" "${f%.fasta}_index"
    fi
done
cd ../

while IFS= read -r sample_name; do
    {
        echo ""
        echo "Preparing consensus sequence for $sample_name"
    } | tee -a "$log_file"
    cd "$sample_name"
    while IFS= read -r ref_name; do
        samtools mpileup -A -d 6000000 -B -Q 0 --reference ../DB_dir/${ref_name}.fasta $sample_name.$ref_name.sorted.bam | ivar consensus -p $sample_name.$ref_name.consensus -t 0.75
    done < "${sample_name}.refs.tsv"
    cd ../
    {
        echo ""
        echo "Finished Processing Sample: $sample_name"
    } | tee -a "$log_file"
done < "$samples_file"

{
    echo ""
    echo "Script completed successfully." | tee -a "$log_file"
} | tee -a "$log_file"
