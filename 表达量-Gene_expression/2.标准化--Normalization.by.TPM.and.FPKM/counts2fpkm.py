import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="id_pick.")
parser.add_argument("-c", "--inputfile", help=" counts file .",
                    default="./counts.xls")
parser.add_argument("-l", "--lens", help="gene(all-exon) length file",
                    default="./lens.txt")
parser.add_argument("-r", "--resultfile", help="fpkm file name.",
                    default="./fpkm.xls")
args = parser.parse_args()

lens = os.path.abspath(args.lens)
counts = os.path.abspath(args.inputfile)
output = os.path.abspath(args.resultfile)

def counts2fpkm(counts, lens, output):
    counts = pd.read_csv(counts,sep='\t',index_col=0)
    anno = pd.read_csv(lens,sep='\t',index_col=0)
    anno.columns = ['len']
    d = anno.to_dict()['len']
    sum_reads = counts.sum()
    lens = pd.DataFrame(d,index=counts.columns).T
    f1 = counts/sum_reads
    fpkm = f1*1000000000/lens
    fpkm.to_csv(output,sep='\t')
    
counts2fpkm(counts, lens, output)
