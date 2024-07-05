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
Just download the scripts (`MuReMa.sh`, `formater.py`, and `grapher.py`) and use them in the same directory where the reads are located.

### Software Tools Requiered
- [bowtie2](https://github.com/BenLangmead/bowtie2)
- [samtools](https://github.com/samtools/samtools)
- [ivar](https://github.com/gkarthik/ivar)
- python3

Be kind and please acknowledge these great authors too!

## Setting up the Inputs
### Compatible reads and sufix
MuReMa works with trimmed paired-end reads. The reads must have the following suffixes: `_1.fastq.gz` for forward reads and `_2.fastq.gz` for reverse reads.

### Samples file format
The samples_file must be a plain text file with one sample name per line as follows:
```
sample-1
sample-2
sample-3
...
sample-n
```
TIP: To create the samples_file in the terminal, use the following commands:
```
ls *_1.fastq.gz > samples_file.txt
sed -i 's/_1.fastq.gz//g' samples_file.txt
```
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
- `-s`, The samples file that contains the names of the samples to be analyzed.
- `-d`, The multifasta file that contains the reference sequences to map to.
```bash
# Basic command
bash MuReMa.sh -s samples_file -d multifasta_file
```
### Optional options
- `-r`, The read length. An integer > 0. Default = 150.
- `-t`, The average depth threshold. An integer > 0. Default = 1000. The average depth threshold is a key value. If the reads aligned to a reference equals or surpasses it, MuReMa will select that reference to obtain a consensus using it as template and generate a graph for it. The threshold could be lower if you work with bacterial reads (for example, `-t 100`) or higher if you are working with amplicons (for example, `-t 3000`).
```bash
# Example for bacterial reads and longer reads
bash MuReMa.sh -s samples_file -d multifasta_file -r 250 -t 100
```
## The Outputs
- A directory per sample containing:
  - Filtered and Refs TSVs: `sample-1.filtered.tsv` includes metrics like mapping reference, reference length, number of mapped reads, and percentage of mapped reads. `sample-1.refs.tsv` lists references surpassing the average depth threshold.
  - BAMs: `sample-1.sorted.bam` for the original alignment and `sample-1.reference-1.sorted.bam` for each individual reference surpassing the threshold.
  - Depth dispersion and Mapped Proportion graphs:
    - X-axis: Number of positions, equal to the reference length.
    - Left Y-axis: Depth dispersion represented by short vertical lines (|), one per mapped position. Color-coded as green if the depth for that position is equal to or higher than the threshold, yellow if below the threshold but higher than 0, and red if the depth is 0.
    - Right Y-axis: Coverage, calculated as the cumulative sum of positions with depth equal to or higher than the threshold, divided by the reference sequence length. Represented by a blue line showing values between 0 and 1.
  - Consensus: Consensus sequences in FA format for each individual reference `sample-1.reference-1.consensus.fa`.
- DB_dir directory:
  - A copy of the original multifasta file used as the database reference sequences and its index files.
  - A TSV file listing references surpassing the average depth threshold.
  - The extracted individual references and their index files.
- samples_file.log: Simple logging system to track script behavior and identify errors.

## Contact
Email me to: facundogcuba@gmail.com or raise an issue

## Acknowledgments
To my co-workers who served as alpha testers, provided ideas, and suggested features for MuReMa.
