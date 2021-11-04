import os
import argparse

parser = argparse.ArgumentParser(description="format_transform.")
parser.add_argument("-fa", "--fasta", help="Input .fa file [default: ./fasta]",
                    default="./fasta")
parser.add_argument("-r", "--resultfile", help="Output file with selected gene [default:./Seq_selected.txt].",
                    default="./Seq_selected.txt")
args = parser.parse_args()

fasta = os.path.abspath(args.fasta)
outputDirectory = os.path.abspath(args.resultfile)
dict={}
fa = open(fasta, 'r')
fr = open(outputDirectory, 'a')

import re
dnaSeq = ''
with open(fasta) as f:
    for line in f:
        line = line.rstrip()
        dnaSeq += line.upper()
rnaSeq1 = re.sub('T', 'U', dnaSeq)
rnaSeq2 = dnaSeq.replace('T', 'U')
print(rnaSeq1)
print(rnaSeq2)
