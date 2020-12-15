# FindFur

Simple Python-based tool for extracting and predicting furin cleavage sites in enveloped virus proteins using HMMER3

## Project details
This project was created on macOS Mojave Version 10.14.6 and IDLE 3.7.2. 

## Requirements
To run this project, the following programs are required:
* Python 3+
* Hmmer3 installed locally (for hmmsearch)
* Pandas
* Biopython

2POSITIVE.HMM needs to be in the same directory as FindFur

## Purposes of script
FindFur_Extract.py
* Purpose: Extract all potential furin cleavage sites in a FASTA file containing the protein sequences a user would like to search with the profile HMM 
* Input: Protein sequences of interest in FASTA format and path to HMMER's hmmsearch 
* Output: Putative furin cleavage site motifs of length 20-residues, HMMER results in domain table format

FindFur_HMMERParse.py
* Purpose: Reports whether motifs of interest are furin cleavage site hits for viral substrates 
* Input: Previous HMMER results and motif FASTA file
* Optional output: Text file of HMMERParse results

## To run:
1. Run protein sequences of sequences through FindFur_Extract
2. Run outputs of FindFur_Extract into FindFur_HMMERParse
