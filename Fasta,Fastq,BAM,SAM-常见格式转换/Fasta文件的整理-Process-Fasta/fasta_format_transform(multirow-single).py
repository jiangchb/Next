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
for line in fa:
    if line.startswith('>'):
        line = line.strip('\n')
        line = line.split("  ", 1)
        #line = line.replace(' ', '')
        #line = line.replace("locus*", '')
        name = str(line[0])
        locus = str(line[1])
        dict[name] = ''
    else:
        dict[name] += line.replace('\n','')
for ID in dict.keys():
        fr.write(ID + "\n" + dict[ID] + "\n")
fr.close()
fa.close()