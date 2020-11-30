# Objective: Extract potential furin cleavage site from protein sequences of interest

import os, sys, argparse, subprocess
import pandas as pd
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqIO.FastaIO import FastaIterator


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


def findIndexes(fasta):
    """Pull indexes to mark the P1, P14 and P6' location"""
    import re

    # For each record, find positions of all Rs ("P1s") in its sequence and store it
    # as a list within a dictionary
    all = {}
    for record in fasta:
        positions = [m.start() for m in re.finditer('R', str(record.seq))]
        all[record.id] = positions

    # Convert into dataframe, add index to define sequence groups
    df = pd.DataFrame(list(all.items()), columns = ['ID', 'P1'])
    num = [*range(len(df))]
    df.insert(0, "Index", num)
    df = df.explode('P1')
    df['ID'] = df.ID.astype(str) + "_" + df.groupby(level=0).cumcount().add(1).astype(str)

    # Calculate P14 and P6' based on P1 position
    df['P14'] = df['P1'].apply(lambda x: x-13 if x > 13 else 0)
    df['PP6'] = df['P1'].apply(lambda x: x+7)
    df['total'] = df['PP6'] - df['P14']
    # df_filtered = df[df['total'] == 20]

    df = df.reset_index(drop = True).dropna()
    
    return(df)

    

def generateMotifs(df, inputfile):
    """Generate motifs based on the indexes previously found."""

    output_handle = "{0}_motif.fasta".format(inputfile.split('.')[0])

    # Retake in fasta file to parse by name
    record = SeqIO.index(inputfile, "fasta")

    col = []
    with open(output_handle, 'w') as writer:
        for i in range(len(df)):
            if df['total'].values[i] == 20:
                rec_name = df['ID'].values[i]
                rec_name1 = rec_name.split('_')[0]
                rec = record[rec_name1]
                n_term = int(df['P14'].values[i])
                c_term = int(df['PP6'].values[i])
                subseq = str(rec.seq)[n_term:c_term]
                writer.write(">{0}\n{1}\n".format(rec_name, subseq))
                col.append(subseq)
            else:
                continue

    #df['Subseq'] = col
    #print(df)
    
    return(df)


def readFasta(inputfile):
    """"Read fasta file"""
    
    with open(inputfile) as handle:
        for seq_record in FastaIterator(handle):
            yield seq_record



def parseArgs():
    parser = argparse.ArgumentParser("Extract putative FCS and run program through hmmsearch")
    
    parser.add_argument('-f','--fasta', help = 'input fasta file')
    parser.add_argument('-p', '--hs_path', help = 'hmmsearch path')
    #parser.add_argument('-t','--txt', action = 'store_true', help = 'output txt file', required = False)
    args = parser.parse_args()

    return(args)


def main():
    args = parseArgs()
    wd = os.getcwd()

    inputfile = args.fasta
    hs_path = args.hs_path
    #txtfile = args.txt
    # hs_path = '/Users/##/opt/miniconda3/bin/hmmsearch'

    record = readFasta(inputfile)
    record_df = findIndexes(record)
    generateMotifs(record_df, inputfile)
    print(record_df)

    hmm_file = "{0}/2POSTIVE.HMM".format(wd)
    file = "{0}_motif.fasta".format(inputfile.split('.')[0])
    print(file)
    output = '{0}-domtblout.txt'.format(file.split('.')[0])
    print(output)

    subprocess.Popen([hs_path, '--domtblout', output, '--F1', '0.08', '--F2', '0.02', '--F3', '0.00008',  hmm_file, file])



if __name__ == "__main__":
    main()
