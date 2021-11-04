# Version 1.1
# Author: xiufeng yang
# Date: 10/11/2017
# Email: yxfhenu@163.com

import argparse
import os
import shlex
import subprocess
import shutil
import sys
import glob

# Define software path 
perl="/home/fanyucai/software/perl/perl-v5.24.1/bin/perl"
bedtools = "/home/fanyucai/software/bedtools/bedtools2/bin/bedtools"
gffread="/home/fanyucai/software/Cufflinks/cufflinks-2.2.1.Linux_x86_64/gffread"
R_path="/home/fanyucai/software/R/R-v3.4.0/bin/R"
# Read the command line arguments.
parser = argparse.ArgumentParser(description=".")
parser.add_argument("-i", "--inputDirectory", help="Input directory with reconstruction gtf. (default: ./reconstruction)", default="./reconstruction")
parser.add_argument("-g", "--refgenome", help="Input genome.fa (force).")
parser.add_argument("-f", "--refgff", help="Input transcriptome_exon.gff (force).")
parser.add_argument("-b", "--refbed", help="Input exon.bed (force).")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Advanced analysis results.[default:./Advanced_analysis]", default="./Advanced_analysis")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
gff = os.path.abspath(args.refgff)
genome = os.path.abspath(args.refgenome)
bed = os.path.abspath(args.refbed)
outputDirectory = os.path.abspath(args.outputDirectory)

# Creat outputDirectory if it is not exist
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)    
os.chdir(outputDirectory)	

#define function
        
def get_candidate(input,genome,gff,bed):
    if not os.path.exists("./tmp"):
          os.makedirs("./tmp")
    os.chdir("./tmp")
    cmd1='awk \'BEGIN{FS=OFS=\"\\t\"} $3==\"mRNA\"{split($9,a,\"[;=]\"); $9=a[2]; print $0}\' %s > gene_names.gff'%(gff)
    cmd2='cut -f1-6 %s > gene.bed'%(bed)
    cmd3='%s -w cuffmerge.fa -g %s %s && awk \'{print $1}\' cuffmerge.fa >tmp1_cuffmerge.fa && mv tmp1_cuffmerge.fa cuffmerge.fa '%(gffread,genome,input+"/"+"stringtie_merged.gtf")
    cmd4='awk -F \"\\t\" \'$3==\"c\"\' %s > c.refmap'%(input+"/"+"cuffcmp.stringtie_merged.gtf.refmap")
    cmd5='awk -F \"\\t\" \'$3==\"u\"\' %s > u.tmap'%(input+"/"+"cuffcmp.stringtie_merged.gtf.tmap")
    cmd6='awk -F \"\\t\" \'$3==\"exon\"\' %s >merged.gtf ' %(input+"/"+"stringtie_merged.gtf")
    cmd7='sh %s/novel_extend_candidate.sh %s '%(sys.path[0],genome)
    print(cmd1)
    subprocess.run(cmd1,shell=True,check=True)
    print(cmd2) 
    subprocess.run(cmd2,shell=True,check=True)
    print(cmd3)
    subprocess.run(cmd3,shell=True,check=True)
    print(cmd4)
    subprocess.run(cmd4,shell=True,check=True)
    print(cmd5)
    subprocess.run(cmd5,shell=True,check=True)
    print(cmd6)
    subprocess.run(cmd6,shell=True,check=True)  
    print(cmd7)
    subprocess.run(cmd7,shell=True,check=True)  
    os.chdir(outputDirectory)
    shutil.rmtree(outputDirectory+"/tmp")

    
    
if __name__ == "__main__":

      get_candidate(inputDirectory,genome,gff,bed)
