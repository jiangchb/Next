

import os
import argparse

parser = argparse.ArgumentParser(description="transform_fastefile_id2genename.")
parser.add_argument("-fb", "--genename2id", help="Input relation of id and genename [default: ./relation.txt]",
                    default="./relation.txt")
parser.add_argument("-fa", "--fasta", help="Input .fa file [default: ./gene.fa]",
                    default="./gene.fa")
                    
parser.add_argument("-r", "--resultfile", help="Output fasta file with the name you want [default:./gene_id.fa].",
                    default="./gene_id.fa")
args = parser.parse_args()
fasta = os.path.abspath(args.fasta)
outputDirectory = os.path.abspath(args.resultfile)
gene_list=os.path.abspath(args.genename2id)
dict2={}
fa = open(fasta, 'r')#fasta
fu = open(gene_list, 'r')#genelist
fr = open(outputDirectory, 'a')#selected_fast
for line in fu:
    line = line.strip('\n')
    line = line.split()
    #line = line.replace(' ', '')
    #line = line.replace("locus*", '')
    name2 = str(line[0])
    content = str(line[1])
    dict2[name2] = content
dict = {}
for line in fa:
    line = line.strip('\n')
    line = line.split()
    #line = line.replace(' ', '')
    #line = line.replace("locus*", '')
    name1 = str(line[1])
    id = str(line[0])
    dict[name1] = id
for name in dict2.keys():
    if dict.has_key(name):
        fr.write(dict[name] +"\t"+ dict2[name]+'\n')
    else:
        fr.write(name + "\t" + dict2[name]+"\n")
        
fr.close()
fu.close()
fa.close()