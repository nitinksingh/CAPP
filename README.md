# CAPP
Coupling Analysis of Protein residues with Phi-mixing Coefficient 

This project aims to predict co-evolutionary coupling among amino acid residues.

# Installation
Download or clone this repository. This python program depends on biopython module that can be installed by
```bash
pip install biopython
```

The program needs an input FASTA file containing the multi-sequence aligned protein sequences. The other parameters i.e. number of bootstrap runs and threshold parameters can be passed as command line arguments to the script. See below.

```bash
./capp.py -h
usage: capp.py [-h] [-o OUTFILE] [-b BS] [-t THRESHOLD] [-v] fasta_file

Coupling analysis of co-evolving protein residues with phi-mixing coefficient.
This method predicts the co-evolutionary interaction between pairwise amino
acid residues of a protein given multi-alignment sequences in FASTA format.

positional arguments:
  fasta_file    Input fasta file

optional arguments:
  -h, --help    show this help message and exit
  -o OUTFILE    Output file name. default: fasta_file.csv
  -b BS         Number of bootstraps to run. default: 100
  -t THRESHOLD  Threshold to drop low confidence interactions. default: 0
  -v            show program's version number and exit
