# MuReMa: Multi-reference Read Mapper

MuReMa is a tool designed to map paired-end reads to multiple reference genomes simultaneously, ideal for comparative genomics and metagenomics studies.

## The Spirit of MuReMa
MuReMa was originally created as an easy way to work with segmented viruses, but during development, it demonstrated a broader range of uses. I have used it to identify genotypic variants of a virus along with the percentage of those variants in the sample pool, isolate specific genes from a complex sample, and obtain consensus sequences of specific organisms in metagenomic samples, among other applications. Please test it for yourself and let me know how MuReMa can be improved.

## Table of Contents
- [Installation](#installation)
- [Setting up the Inputs](#setting-up-the-inputs)
- [Usage](#usage)
- [The Outputs](#the-outputs)
- [Features](#features)
- [Configuration](#configuration)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
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
### Samples file format
### Multifasta file

## Usage
### Must provide options
- `-s`, 
- `-d`,
```bash
# Basic command
bash MuReMa.sh -s samples_names_file -d multifasta_file
```
### Optional options
- `-r`,
- `-m`,
## The Outputs
### Directory structure
### The files
#### Filtered and Refs TSVs
#### BAMs
#### Depth dispersion and Coverage Proportion graphs
#### Consensus made
