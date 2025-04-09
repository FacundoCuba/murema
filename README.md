# MuReMa: Multi-reference Read Mapper
MuReMa is a tool designed to map paired-end reads to multiple reference genomes simultaneously, ideal for comparative genomics and metagenomics studies.

## The Spirit of MuReMa
MuReMa was originally created as an easy way to work with segmented viruses, but during development, it demonstrated a broader range of uses. I have used it to identify genotypic variants of a virus along with the percentage of those variants in the sample pool, isolate specific genes from a complex sample, and obtain consensus sequences of specific organisms in metagenomic samples, among other applications. Please test it for yourself and let me know how MuReMa can be improved.

## Table of Contents
- [Installation](#installation)
- [Setting up the Inputs](#setting-up-the-inputs)
- [How to use MuReMa](#how-to-use-murema)
- [The Outputs](#the-outputs)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

## Installation
Just download the scripts (`MuReMa.sh`, `formater.py`, and `grapher.py`).

### Software Tools Requiered
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [ivar](https://github.com/gkarthik/ivar)
- python3

Be kind and please acknowledge these great authors too!

## Setting up the Inputs
### Compatible reads
MuReMa works with Illumina paired-end reads.

### Multifasta file format
The multifasta_file is the backbone of MuReMa. The selection of references will impact the output of this tool. The header of each reference should be as short as possible and edited to remove invalid characters for a bash environment. The multifasta_file must be a concatenation of fastas as follows:
```
>reference-1
sequence of refence 1
>reference-2
sequence of refence 2
>reference-3
sequence of refence 3
...
>reference-n
sequence of refence n
```

## How to use MuReMa
### Must provide options
- `-1`, Forward read OR path to forward read.
- `-2`, Reverse read OR path to reverse read.
- `-d`, The multifasta file that contains the reference sequences to map to.
```bash
# Basic command
./MuReMa.sh -1 R1.fastq.gz -2 R2.fastq.gz -d multifasta_file.fasta
```

### Optional options
- `-b`, Bed file OR path to bed file.
- `-r`, The read length. An integer > 0. Default = 150.
- `-t`, The average depth threshold. An integer > 0. Default = 1000. The average depth threshold is a key value. If the reads aligned to a reference equals or surpasses it, MuReMa will select that reference to obtain a consensus using it as template and generate a graph for it. The threshold could be lower if you work with bacterial reads (for example, `-t 100`) or higher if you are working with amplicons (for example, `-t 3000`).
- `-T`, Specify that the provided reads are already trimmed.
- `-U`, Specify that the provided reads are untrimmed (default).
- `-v`, Display the version of the script.
- `-h`, Display help message.
```bash
# Example of a loop for multiple reads
for f in *_1.fastq.gz; do ./MuReMa.sh -1 $f -2 ${f%R1.fastq.gz}R2.fastq.gz -d multifasta_file.fasta; done
```

## The Outputs
- A directory per sample containing:
  - Filtered and Refs TSVs: `sample-1.filtered.tsv` includes metrics like mapping reference, reference length, number of mapped reads, and percentage of mapped reads. `sample-1.refs.tsv` lists references surpassing the average depth threshold.
  - BAMs: `sample-1.sorted.bam` for the original alignment and `sample-1.reference-1.sorted.bam` for each individual reference surpassing the threshold.
  - Log file: to track script behavior and identify errors `sample-1.log`.
  - Consensus: Consensus sequences in FA format for each individual reference `sample-1.reference-1.consensus.fa`.
  - Depth dispersion and Mapped Proportion graphs:
    - X-axis: Number of positions, equal to the reference length.
    - Left Y-axis: Depth dispersion represented by short vertical lines (|), one per mapped position. Color-coded as green if the depth for that position is equal to or higher than the threshold, yellow if below the threshold but higher than 0, and red if the depth is 0.
    - Right Y-axis: Coverage, calculated as the cumulative sum of positions with depth equal to or higher than the threshold, divided by the reference sequence length. Represented by a blue line showing values between 0 and 1.
- DB_dir directory:
  - A copy of the original multifasta file used as the reference sequences database and its index files.
  - The extracted individual references and their index files.

## Contact
Email me to: facundogcuba@gmail.com or raise an issue

## Acknowledgments
To my co-workers who served as alpha testers, provided ideas, and suggested features for MuReMa.
