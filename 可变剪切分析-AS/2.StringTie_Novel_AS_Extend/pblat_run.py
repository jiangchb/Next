# Version 1.1
# Author: xiufeng yang
# Date: 10/11/2017
# Email: yxfhenu@163.com

import argparse
import glob
import os
import subprocess
import shutil
import sys
import configparser

# Define software path 
pblat="/home/fanyucai/software//pblat/pblat-2.0.0/pblat"
R_path="/home/fanyucai/software/R/R-v3.4.0/bin/R"

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates pblat script.")
parser.add_argument("-i", "--inputDirectory", help="Input directory with Extend_gene candidate (default: ./Advanced_analysis/Extend_gene)", default="./Advanced_analysis/Extend_gene")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Advanced analysis results.[default:./Advanced_analysis/Extend_gene]", default="./Advanced_analysis/Extend_gene")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)

# Creat outputDirectory if it is not exist
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)    
os.chdir(outputDirectory)	

def  pblat_run(input):
    os.chdir(outputDirectory)
    cmd1='%s threads=10 %s %s output.psl'%(pblat,input+"/Extendgene.fa",input+"/Extend.ref.gene.fa") 
    cmd2='cut -f 10,11,12,13,14,15,16,17 output.psl > blat_results_tmp.txt'
    cmd3='awk -F\'\\t\' -v OFS=\'\\t\' \'NR==FNR{a[$1$3]=$1;b[$1$3]=$3;next}{$1=substr($1,1,index($1,\":\")-1);$5=substr($5,1,index($5,\":\")-1);if($1==a[$1$5]) print $0 }\'  Extendgene_candidate.txt blat_results_tmp.txt |sed \'1,5d\' >blat_results.txt' 
    cmd4='%s CMD BATCH %s'%(R_path,sys.path[0]+"/getExtendGene.r")
    print(cmd1)
    subprocess.run(cmd1,shell=True,check=True)
    print(cmd2) 
    subprocess.run(cmd2,shell=True,check=True)
    print(cmd3)
    subprocess.run(cmd3,shell=True,check=True)
    print(cmd4)
    subprocess.run(cmd4,shell=True,check=True)
if __name__ == "__main__":
    pblat_run(inputDirectory)
    
   
