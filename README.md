# Spike-in_analysis
This script is designed to facilitate processing of 16S sequencing data with Spikes, allowing for the calculation of "absolute abundance".

# WARNING - THIS PAGE IS UNDER CONSTRUCTION AND THE SCRIPT IS CURRENTLY IN BETA

## Installation 

To use this script on your own data you must have anaconda/miniconda installed with python 3.7, biopython (tested with 1.7.3), pandas (tested with 0.24.2) and numpy (tested with 1.16.3). You can install miniconda [here](https://docs.conda.io/en/latest/miniconda.html#macosx-installers). 

Once you have conda, clone this git repository, and run the following to install all required packages: 

```
conda env create -f environment.yml
```

Additionally, BLAST is required to determine which OTU sequences are spikes. Instructions for installing BLAST can be found [here](https://www.ncbi.nlm.nih.gov/books/NBK52637/). Make sure to check the BLAST cli application is in your path before starting by opening a terminal and running:

```
blastn -h
``` 

If you see the help message then you are good to go! 

Once everything is installed copy your OTU table, mapping file, fasta file with OTU sequences and fast file containing the spike sequences to the same directory as the process_abs_data.py. 

## Basic Usage 

Expected input:

Sequences in fasta format and OTU and mapping files as tab delimited text. The -ng option is the amount of spike-in added to samples before DNA extraction in ng, while -w represents the column in the mapping file containing sample weights. 

```
python process_abs_data.py -f OTUs-Seqs.fasta -otu OTUs-Table.tab, -m mapping_file.tab -s spikes.fasta -ng 1 -w "sample_weights"
``` 

