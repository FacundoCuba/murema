#!/bin/bash

# Function to display script usage
function display_usage() {
    echo "Usage: $0 -1 FASTQ -2 FASTQ -d FILE [-b FILE] [-r INTEGER] [-t INTEGER] [-T|-U] [-m] [-v] [-h]"
    echo ""
    echo "Options:"
    echo "  -1 FASTQ   Forward read OR path to forward read"
    echo "  -2 FASTQ   Reverse read OR path to reverse read"
    echo "  -d FILE    Multifasta file OR path to the multifasta file"
    echo "  -b FILE    Bed file OR path to bed file (optional)"
    echo "  -r INTEGER Read length (default: 150)"
    echo "  -t INTEGER Average Depth threshold (default: 1000)"
    echo "  -T         Specify that the provided reads are already trimmed"
    echo "  -U         Specify that the provided reads are untrimmed (default)"
    echo "  -m         Specify that the consensus sequence will be called using the specifically aligned reads"
    echo "  -v         Display the version of the script"
    echo "  -h         Display this help message"
    exit 1
}

# Display version
function display_version() {
    echo "Script Version: 4.0"
    exit 0
}

# Initialize variables with default values
multifasta_file=""
bed_file=""
read_length=150
avg_depth_threshold=1000
sample_name=""
log_file=""
DB_log_file=""
version="4.0"
trim_reads=true
meta=false

# Parse options using getopts
while getopts "1:2:d:b:r:t:TUmvh" option; do
    case "$option" in
        1) r1="$PWD/$OPTARG";;
        2) r2="$PWD/$OPTARG";;
        d) multifasta_file="$PWD/$OPTARG";;
        b) bed_file="$PWD/$OPTARG";;
        r) read_length="$OPTARG";;
        t) avg_depth_threshold="$OPTARG";;
        T) trim_reads=false;;
        U) trim_reads=true;;
        m) meta=true;;
        v) display_version;;
        h) display_usage;;
        :) printf "Missing argument for -%s\n" "$OPTARG" >&2; display_usage >&2;;
        \?) printf "Illegal option: -%s\n" "$OPTARG" >&2; display_usage >&2;;
    esac
done

# Check if required arguments are provided
if [[ -z "$r1" || -z "$r2" || -z "$multifasta_file" ]]; then
    echo "Error: Arguments -1, -2, and -d must be provided." >&2 | tee -a "$log_file"
    display_usage
fi

# Set sample name
if [[ $trim_reads == true ]]; then
    sample_path="${r1%_S*_R1_001.fastq.gz}"
    sample_name=$(basename "$sample_path")
else
    sample_path="${r1%_*1.fastq.gz}"
    sample_name=$(basename "$sample_path")
fi

# Create sample directory and log file
mkdir -p "$sample_name"
cd "$sample_name"
log_file="${sample_name}.log"

# Log initial arguments
{
    echo "Script invoked with the following arguments:"
    echo "Sample name: $sample_name"
    echo "Forward read: $r1"
    echo "Reverse read: $r2"
    echo "Multifasta DB: $multifasta_file"
    echo "Bed file: $bed_file"
    echo "Read length: $read_length"
    echo "Average Depth threshold: $avg_depth_threshold"
    echo "Trimming reads: $trim_reads"
    echo "Metagenomic called consensus: $meta"
} > "$log_file"
cd ../

# Inputs validation
for file in "$r1" "$r2" "$multifasta_file"; do
    if [[ ! -e "$file" ]]; then
        echo "Error: $file is not a valid file." >&2 | tee -a "$log_file"
        exit 1
    fi
done

if ! [[ "$read_length" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $read_length is not a valid integer greater than 0." >&2 | tee -a "$log_file"
    exit 1
fi

if ! [[ "$avg_depth_threshold" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: $avg_depth_threshold is not a valid integer greater than 0." >&2 | tee -a "$log_file"
    exit 1
fi

# Check for required tools
for tool in trim_galore bowtie2 samtools python3 ivar; do
    command -v $tool >/dev/null 2>&1 || { echo "$tool is required but it's not installed. Aborting." >&2 | tee -a "$log_file"; exit 1; }
done

# Create DB directory, log file and index the multifasta
mkdir -p DB_dir
cp "$multifasta_file" DB_dir/murema_DB.fasta
cd DB_dir
DB_log_file="DB_dir.log"
if [ ! -f "murema_DB_index.1.bt2" ]; then
    bowtie2-build -f murema_DB.fasta murema_DB_index || { echo "Error: bowtie2-build failed." | tee -a "$DB_log_file"; exit 1; }
fi
cd ../

# Check if reads need trimming
cd "$sample_name"
if [[ $trim_reads == true ]]; then
    # Perform trimming
    echo "Trimming reads for sample $sample_name" | tee -a "$log_file"
    trim_galore -j 8 -q 30 --paired -o ./ --no_report_file "$r1" "$r2" --basename "$sample_name"
    mv "${sample_name}_val_1.fq.gz" "${sample_name}_1.fastq.gz"
    mv "${sample_name}_val_2.fq.gz" "${sample_name}_2.fastq.gz"
    r1_trimmed="${sample_name}_1.fastq.gz"
    r2_trimmed="${sample_name}_2.fastq.gz"
else
    # Skip trimming
    echo "Skipping trimming step as reads are already trimmed." | tee -a "$log_file"
    r1_trimmed="$r1"
    r2_trimmed="$r2"
fi

# Align reads and sort BAM
echo "Aligning reads and generating sorted BAM for $sample_name" | tee -a "$log_file"
bowtie2 --end-to-end --very-sensitive -p 8 -x "../DB_dir/murema_DB_index" -1 "$r1_trimmed" -2 "$r2_trimmed" | \
samtools sort | samtools view -@ 8 -b -F 4 -q 1 -o "${sample_name}.sorted.bam" || { echo "Error: bowtie2 or samtools sorting failed" | tee -a "$log_file"; exit 1; }
samtools index "${sample_name}.sorted.bam"
samtools idxstats "${sample_name}.sorted.bam" > "${sample_name}.tsv"

# Check if optional bed_file exists
if [ -n "$bed_file" ] && [ -f "$bed_file" ]; then
    echo "Trimming primers using bed file: $bed_file" | tee -a "$log_file"
    ivar trim -e -i "${sample_name}.sorted.bam" -b "${bed_file}" -p "${sample_name}.primertrim"
    samtools sort "${sample_name}.primertrim.bam" -o "${sample_name}.primertrim.sorted.bam"
else
    echo "No bed file provided or file not found; skipping primer trimming." | tee -a "$log_file"
fi

# Run formater.py
formater.py "${sample_name}.tsv" "$sample_name" "$read_length" "$avg_depth_threshold" || { echo "formater.py failed" | tee -a "$log_file"; exit 1; }

# Convert Windows line endings
sed -i 's/\r$//' "${sample_name}.filtered.tsv" "${sample_name}.refs.tsv"
cd ../

# Extract sequences from DB
cd DB_dir
while IFS= read -r ref_name; do
    awk -v RS='>' -v ref="$ref_name" '$1 == ref { print ">"$0; exit }' murema_DB.fasta > "${ref_name}.fasta"
done < "../${sample_name}/${sample_name}.refs.tsv"

# Index references
for f in *.fasta; do
    [ ! -f "${f%.fasta}_index.1.bt2" ] && bowtie2-build -f "$f" "${f%.fasta}_index" | tee -a "$DB_log_file"
done
cd ../

# Consensus & Graphing
cd "$sample_name"
while IFS= read -r ref_name; do
    bowtie2 --end-to-end --very-sensitive -p 8 -x "../DB_dir/${ref_name}_index" -1 "$r1_trimmed" -2 "$r2_trimmed" | \
    samtools sort | samtools view -@ 8 -b -F 4 -q 1 -o "${sample_name}.${ref_name}.sorted.bam"
    if [ -n "$bed_file" ] && [ -f "$bed_file" ]; then
        echo "Creating consensus sequence with primer trimming using $ref_name" | tee -a "$log_file"
        ivar trim -e -i "${sample_name}.${ref_name}.sorted.bam" -b "${bed_file}" -p "${sample_name}.${ref_name}.primertrim"
        samtools sort "${sample_name}.${ref_name}.primertrim.bam" -o "${sample_name}.${ref_name}.primertrim.sorted.bam"
        if [[ $meta == true ]]; then
            samtools view  -@ 8 -b -F 4 -q 1 "${sample_name}.primertrim.sorted.bam" "$ref_name" -o "${sample_name}.${ref_name}_filtered.primertrim.sorted.bam"
            samtools index "${sample_name}.${ref_name}_filtered.primertrim.sorted.bam"
            samtools mpileup -A -d 6000000 -B -Q 0 -q 20 --reference "../DB_dir/${ref_name}.fasta" "${sample_name}.${ref_name}_filtered.primertrim.sorted.bam" | \
                ivar consensus -p "${sample_name}.${ref_name}_filtered.primertrim.consensus" -t 0.75
            grapher.py "${sample_name}.${ref_name}_filtered.primertrim.sorted.bam" "$sample_name" "$ref_name" "$avg_depth_threshold" || { echo "grapher.py failed" | tee -a "$log_file"; exit 1; }
        else
            samtools index "${sample_name}.${ref_name}.primertrim.sorted.bam"
            samtools mpileup -A -d 6000000 -B -Q 0 -q 20 --reference "../DB_dir/${ref_name}.fasta" "${sample_name}.${ref_name}.primertrim.sorted.bam" | \
                ivar consensus -p "${sample_name}.${ref_name}.primertrim.consensus" -t 0.75
            grapher.py "${sample_name}.${ref_name}.primertrim.sorted.bam" "$sample_name" "$ref_name" "$avg_depth_threshold" || { echo "grapher.py failed" | tee -a "$log_file"; exit 1; }
        fi
    else
        echo "Creating consensus sequence without primer trimming using $ref_name" | tee -a "$log_file"
        if [[ $meta == true ]]; then
            samtools view  -@ 8 -b -F 4 -q 1 "${sample_name}.sorted.bam" "$ref_name" -o "${sample_name}.${ref_name}_filtered.sorted.bam"
            samtools index "${sample_name}.${ref_name}_filtered.sorted.bam"
            samtools mpileup -A -d 6000000 -B -Q 0 -q 20 --reference "../DB_dir/${ref_name}.fasta" "${sample_name}.${ref_name}_filtered.sorted.bam" | \
                ivar consensus -p "${sample_name}.${ref_name}_filtered.consensus" -t 0.75
            grapher.py "${sample_name}.${ref_name}_filtered.sorted.bam" "$sample_name" "$ref_name" "$avg_depth_threshold" || { echo "grapher.py failed" | tee -a "$log_file"; exit 1; }
        else
            samtools index "${sample_name}.${ref_name}.sorted.bam"
            samtools mpileup -A -d 6000000 -B -Q 0 -q 20 --reference "../DB_dir/${ref_name}.fasta" "${sample_name}.${ref_name}.sorted.bam" | \
                ivar consensus -p "${sample_name}.${ref_name}.consensus" -t 0.75
            grapher.py "${sample_name}.${ref_name}.sorted.bam" "$sample_name" "$ref_name" "$avg_depth_threshold" || { echo "grapher.py failed" | tee -a "$log_file"; exit 1; }
        fi
    fi
done < "${sample_name}.refs.tsv"
echo "Script completed successfully." | tee -a "$log_file"
cd ../
