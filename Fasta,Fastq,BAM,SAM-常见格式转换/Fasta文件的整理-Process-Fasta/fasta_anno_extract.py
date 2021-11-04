#####
#根据列表gene信息，从gene.fa中，提取相应的，做成fasta文件
#####
import os
import argparse

parser = argparse.ArgumentParser(description="id_anno_pick.")
#parser.add_argument("-g", "--inputfile", help=" gene list file.[default : ./Gene_list] .",
#                    default="./Gene_list")
parser.add_argument("-fa", "--fasta", help="Input .fa file [default: ./fasta]",
                    default="./fasta")
parser.add_argument("-r", "--resultfile", help="Output file with selected gene [default:./Seq_selected.txt].",
                    default="./anno.txt")
args = parser.parse_args()

#inputDirectory = os.path.abspath(args.inputfile)
fasta = os.path.abspath(args.fasta)
outputDirectory = os.path.abspath(args.resultfile)
dict={}
dict2 = {} 
fa = open(fasta, 'r')#fasta
#fu = open(inputDirectory, 'r')#genelist
fr = open(outputDirectory, 'w')#annotation
for line in fa:
    if line.startswith('>'):
        line = line.strip('\n')
        line = line.split(" ", 1)
        #line = line.replace(' ', '')
        #line = line.replace("locus*", '')
        name = str(line[0])
        locus = str(line[1])
        fr.write(name + "\t" + locus + "\n")
        dict2[name] = locus 
        dict[name] = ''
    else:
        dict[name] += line.replace('\n','')
#print dict2
#print dict
"""
fa.close()
for line in fu:
    line = line.strip('\n')
    line = line.replace(' ', '')
    Cho = ">" + line
    STR_Cho = str(Cho)
    #print STR_Cho
    for id in dict2.keys():
        if id == STR_Cho:
            fr.write(id + "  " + dict2[id] + "\n")
            for ID in dict.keys():
                if ID == STR_Cho:
                    fr.write(dict[ID] + "\n")
fu.close()
"""

fr.close()