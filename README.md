# kaks - Ka/Ks analysis tool

Performs a Ka/Ks analysis - i.e. determine the ratio of the number of non-synonymous (Ka) to the number of synonymous (Ks) mutations of a given set of sequences and a set of corresponding reference sequences. In order for this to work, the sequences and the reference need to have a common ancestor.

`kaks` requires a FASTA file containing the query (nucleotide) sequences and a set of one or more FASTA files containing the reference sequences. It can be used as a command line tool (see usage below) or as a python module (please refer to the code for the API).

`kaks` needs NCBI `blastp` to be installed and callable via the command line.

## Mode of operation

`kaks` first translates the nucleotide sequences and performs a reciprocal best blastp hit step to determine sequence correspondence between the query and reference (based on the lowest evalue on the forward and reverse blast run). It then uses the blastp alignment on the reverse step to determine the number of non-synonymous changes (based on the amino-acid sequence). It then transfers the alignment over to the nucleotide sequences to determine the number of synonymous changes. The result is collated into a table with the Ka/Ks value listed per query sequence and per reference sequence FASTA file given.

## Installation

`kaks` requires NCBI `blastp` to be installed. By far the easiest way is via conda by running `conda install -c bioconda blast`. After this you can run
```
pip install kaks
```
to install the package.

## Usage

Simply run
```
kaks query.fasta reference1.fasta reference2.fasta > output.csv
```
to perform a Ka/Ks analysis on the query sequences in `query.fasta` with the reference sequences in `reference1.fasta` and `reference2.fasta`. You can specify any number of reference sequence files and the tool will run a separate analysis for each. The results are given as the Ka/Ks ratio for each sequence in `query.fasta` as the rows and each reference sequence file as the columns. Missing results (if a given sequence has no reciprocal best blast hit in the reference) are indicated by `NaN` entries.

```
usage: kaks [-h] query reference [reference ...]

Perform a Ka/Ks analysis for the given nucleotide query and reference FASTA files. The result is given in CSV format to stdout.

positional arguments:
  query       The filename of the query FASTA.
  reference   The filename(s) of the reference FASTA.

options:
  -h, --help  show this help message and exit
  ```