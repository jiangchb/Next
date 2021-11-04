import argparse
import glob
import subprocess
import shutil
import os
import configparser
import sys
import numpy as np
import pandas as pd
import pybedtools
import tempfile
from Bio import SeqIO
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

# Read the command line arguments.
parser = argparse.ArgumentParser(description="gene coverage,genome coverage")
parser.add_argument("-b1", "--inputfile1", help="Input bam file")
parser.add_argument("-b2", "--inputfile2", help="Input genome bed file.default: ./genome/exon.bed")
parser.add_argument("-b3", "--inputfile3", help="Input gene gtf file.default:./genome/transcroptome_exon.gff")
parser.add_argument("-b4", "--inputfile4", help="Input genome fasta.default:./genome/genome.fa")
parser.add_argument("-p", "--output_prefix", help="output file's prefix,eg: Sample_A")
parser.add_argument("-o", "--outputDirectory", help="Output directory ")
args = parser.parse_args()

# Process the command line arguments.
bamfile = os.path.abspath(args.inputfile1)
bedfile = os.path.abspath(args.inputfile2)
gfffile = os.path.abspath(args.inputfile3)
genomefile = os.path.abspath(args.inputfile4)
prefix = args.output_prefix
outputDirectory = os.path.abspath(args.outputDirectory)

if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

os.chdir(outputDirectory)


# define function for genome reads coverage


def transcript_coverage(bed, gff, bam, output, prefix):
    if not os.path.isfile(output + "/" + prefix + ".coverage.chart.pdf"):
        ##############1.get gene length#################
        odata = pd.read_csv(bed, header=None, sep='\t',
                              names=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"], low_memory=False)
        if not os.path.isfile(output + "/" + "transcript.length.txt"):
            outputfile = open(output + "/" + "transcript.length.txt", 'a')
            for i in range(0, len(odata)):
                exon_len = odata.iloc[i, 10].split(",")
                del exon_len[len(exon_len) - 1]
                exon_len = [int(x) for x in exon_len]
                trancript_len = sum(exon_len)
                name = odata.iloc[
                    i, 3]
                outputfile.write(str(name) + "\t" + str(trancript_len) + "\n")
            outputfile.close()

        ##############2.get genomereads depth 
        command1 = "bedtools genomecov -ibam " + bam + " -bg >" + output + "/" + prefix + ".depth.txt"
        subprocess.run(command1, shell=True, check=True)

        command2 = "bedtools intersect -wb -wa -b " + output + "/" + prefix + ".depth.txt -a " + gff + " > " + output + "/" + prefix + ".intersect.txt"
        subprocess.run(command2, shell=True, check=True)

        command3 = "awk -F'\\t'  '$3==\"exon\"{print substr($9,8),$4,$5,$11,$12}' " + output + "/" + prefix + ".intersect.txt " + "|awk '{if($3<=$5) a=$3; else a=$5;if($2<=$4) b=$4; else b=$2; print $0,a,b}'|awk ' {a[NR]=$4; b[NR]=$5}(NR>1){n=$6-$7; if(a[NR]!=b[NR-1]) n=n+1; print $0,n}'|awk '{a[$1]+=$8} END{for(i in a) print i,a[i]}'" + " > " + output + "/" + prefix + ".transcript.coverage.txt"
        command4 = "Rscript " + sys.path[
            0] + "/coverage.piechart.r  " + output + "/" + "transcript.length.txt " + output + "/" + prefix + ".transcript.coverage.txt " + prefix
        subprocess.run(command3, shell=True, check=True)
        subprocess.run(command4, shell=True, check=True)
    if os.path.isfile(output + "/" + prefix + ".coverage.chart.pdf"):
        os.remove(output + "/" + prefix + ".depth.txt")
        os.remove(output + "/" + prefix + ".intersect.txt")
        os.remove(output + "/" + prefix + ".transcript.coverage.txt")
        os.remove(output + "/" + "transcript.length.txt")


def genome_coverage(genome, bed, bam, output, prefix):
    if not os.path.isfile(output + "/" + prefix + ".genome_distribution_plot.pdf"):
        if not os.path.isfile(output + "/" + prefix + "/align.bed"):
            command5 = "bedtools bamtobed -i " + bam + " >" + output + "/" + prefix + ".align.bed"
            subprocess.run(command5, shell=True,check=True)
            command6 = "perl " + sys.path[0] + "/cal_Reads_genes.pl " + genome + " " + output + "/" + prefix + ".align.bed " + bed + " 100000 " + output + "/" + prefix
            subprocess.run(command6, shell=True,check=True)
           
    if os.path.isfile(output + "/" + prefix + ".genome_distribution_plot.pdf"):
        os.remove(output + "/" + prefix + ".align.bed")
        os.remove(output + "/" + prefix + ".reads.in.window.out")
        os.remove(output + "/" + prefix + ".gene.in.window.out")
        os.remove(output + "/" + "genome.window.100000.bed")
        os.remove(output + "/" + "genome.chrsize.txt")


# exe functions
if __name__ == "__main__":
    transcript_coverage(bedfile, gfffile, bamfile, outputDirectory, prefix)
    genome_coverage(genomefile, bedfile, bamfile, outputDirectory, prefix)
