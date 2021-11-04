# Version 1.1
# Author: xiufeng yang
# Date: 10/05/2017
# Email: yxfhenu@163.com

import argparse
import os
import shlex
import subprocess
import shutil
import sys

# Define software path 
bwa="/home/fanyucai/software/bwa/bwa-0.7.12/bwa"
samtools="/home/fanyucai/software/samtools/samtools-v1.4/bin/samtools"
hisat2_path="/home/fanyucai/software/HISAT2/hisat2-2.1.0"

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Building hisat index for genome.")
parser.add_argument("-g", "--refgenome", help="Input genome fasta (force).")  
parser.add_argument("-n", "--threads", help="number of threads for hisat2 [default:8]", default="8")
args = parser.parse_args()

# Process the command line arguments.
genome = os.path.abspath(args.refgenome)
thread = args.threads

#define function for hisat2

def hisat2_build(fasta_file_name,threads):
    #check if index is already build, if not build index of bwa
    if not os.path.exists('%s.6.ht2'%(fasta_file_name.replace(".fa",""))) :
        cmd = '%s/hisat2-build -p %s  %s  %s'%(hisat2_path,threads,fasta_file_name,fasta_file_name.replace(".fa",""))
        print(cmd)
        subprocess.run(cmd,shell=True,check=True)
    else: 
        print("hisat2 index already exists")

if __name__ == "__main__":

      hisat2_build(genome,thread)
