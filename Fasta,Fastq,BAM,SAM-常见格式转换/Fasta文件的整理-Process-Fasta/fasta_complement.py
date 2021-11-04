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
n = 0
complement = {'A':'T','G':'C','C':'G','T':'A'}

with open(outputDirectory, 'a') as f1:
    with open(fasta, 'w') as f:
       for dna_seq in f:
           n +=1
           if n < 100:
               rev_seq = ''
               if dna_seq.startswith(">"):
                   f1.write(dna_seq)
               else:
                   dna_seq = list(dna_seq.strip())
                   for i in dna_seq:
                       rev_seq += complement[i]
                       rev_seq = rev_seq[::-1]
                   f1.write(rev_seq)
                   f1.write('\n')
                   