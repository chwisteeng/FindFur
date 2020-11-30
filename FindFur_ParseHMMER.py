# Objective: Take and parse HMMER3 domblout result
# Takes in: (1) HMMER3 domblout txt file
#           (2) Determine e-value
#           (3) Respective fasta file

from Bio import SearchIO
from Bio import SeqIO
import pandas as pd
import sys, argparse

# subprocess.call(['hmmscan', '--domtblout', 'output.txt', '/location/to/profile.hmm', 'input.proteins.fasta']


def parseData(inputfile, fastafile):
    """Parse through HMMER3 generated domblout file"""

    count = 0
    
    with open(inputfile, 'r') as input:
        for qresult in SearchIO.parse(inputfile, 'hmmsearch3-domtab'):
            query_id = qresult.id      # fasta id
            hits = qresult.hits        # hits 
            acc = qresult.hsps         # domains
            num_hits = len(hits)

            empty = []
            if num_hits > 0:
                for i in range(0,num_hits):
                    org = []
                    
                    empty.append(org)
                    
                    hit_name = hits[i].id         # hit id
                    hit_description = hits[i].description
                    hit_evalue = hits[i].evalue   # evalue
                    hit_bitscore = hits[i].bitscore  #hit bit score 
                    hmm_name = hits[i].accession  
                    acc_start = acc[i].hit_start  # hit start
                    acc_end = acc[i].hit_end      # hit

                    hooray = 'HIT!' if hit_evalue < 0.0023 else 'prob not'
                    if hooray == 'HIT!':
                        count += 1


                    #if ethreshold == 0:
                    #    hooray = "no e-value set"
                    # hooray = 'HIT!' if hit_evalue < 0.024 else 'Maybe not'
                    
                    org.extend((hit_name, hit_description, hit_evalue, hit_bitscore, acc_start, acc_end, hooray))

    df = pd.DataFrame(data=empty, columns=['Hit', 'Description', 'Evalue', 'Bitscore', 'Start', 'End', 'FCS?'])
    
    record = (SeqIO.index(fastafile, "fasta")) # load fasta file

    # collect motif sites from sequences
    motifs = []
    for i in range(len(df)):
        rec_name = df['Hit'].values[i]
        rec = record[rec_name]                # use hit id to find query
        qstart = df['Start'].values[i]        
        qend = df['End'].values[i]
        rec_seq = str(rec.seq)[qstart:qend]
        motifs.append(str(rec_seq))

    df.insert(5, "Align Region", motifs)

    return(df, record, count)

def parseArgs():
    parser = argparse.ArgumentParser("Parse HMMER domblout result and report putative FCS")
    parser.add_argument('-d', '--domblout', help = 'input HMMER3 domblout file')
    # parser.add_argument('-e','--evalue', type=float, help = 'determine e-value threshold')
    parser.add_argument('-f','--fasta', help = 'input respective fasta file')
    parser.add_argument('-t','--txt', action = 'store_true', help = 'output txt file', required = False)
    args = parser.parse_args()

    return(args)

def main():
    args = parseArgs()

    inputfile = args.domblout
    # ethreshold = args.evalue
    fastafile = args.fasta
    txtfile = args.txt

    df, record, count = parseData(inputfile, fastafile)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    pd.set_option('expand_frame_repr', False)

    print(df)
    print("\n")

    print("Number of queries: ", len(record))
    print("Number of hits: ", len(df))
    print("Number of potential FCS: ", count)

    if txtfile:
        input_handle = inputfile.split(".")[0]
        output_handle = '%s_parser' %input_handle
        df.to_csv('%s.txt' %output_handle, sep='\t')
        print('\n%s.txt is done' % output_handle)

        
if __name__ == "__main__":
    main()
