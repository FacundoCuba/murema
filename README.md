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
Just download the scripts (`MuReMa.sh`, `formater.py`, and `grapher.py`) and save them in the same directory where the reads are located.

### Software Tools
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
### Multifasta file
The multifasta_file is the backbone of MuReMa. The selection of references will impact the output of this tool. The header of each reference should be as short as possible and edited to remove invalid characters for a bash environment. The multifasta_file must be concatenation of fastas as follows:
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
- `-m`, The mean depth threshold. An integer > 0. Default = 1000.
## The Outputs
### Directory structure
### The files
#### Filtered and Refs TSVs
#### BAMs
#### Depth dispersion and Coverage Proportion graphs
#### Consensus made

## Contact
facundogcuba@gmail.com
## Acknowledgments
