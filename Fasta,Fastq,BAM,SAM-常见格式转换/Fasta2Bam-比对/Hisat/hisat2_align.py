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
parser = argparse.ArgumentParser(description="Using hisat2 to align RNAseq fastq to reference genome,and convert sam to bam .")
parser.add_argument("-pe1","--paired_reads_5",help="5' reads fastq (force).")
parser.add_argument("-pe2","--paired_reads_3",help="3' reads fastq (force).")
parser.add_argument("-g", "--refgenome", help="Input genome fasta (force).") 
parser.add_argument("-l", "--librarytype", help="Specify strand-specific information(force).1)RF: fr-firststrand,dUTP; 2)FR: fr-secondstrand;  3)unstranded:not strand-sepcific. ('F' means a read corresponds to a transcript. 'R' means a read corresponds to the reverse complemented counterpart of a transcript)",choices=["RF", "FR", "unstranded"])
parser.add_argument("-o", "--outputDirectory", help="Output directory with bam results.")
parser.add_argument("-p", "--prefix", help="the prefix of output sam (force).") 
parser.add_argument("-n", "--threads", help="number of threads for hisat2 [default:8]", default="8")
args = parser.parse_args()

# Process the command line arguments.
pe1 = os.path.abspath(args.paired_reads_5)
pe2 = os.path.abspath(args.paired_reads_3)
genome = os.path.abspath(args.refgenome)
outputDirectory = os.path.abspath(args.outputDirectory)
library = args.librarytype
outname = args.prefix
thread = args.threads

# Creat outputDirectory if it is not exist
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)    
os.chdir(outputDirectory)	

#define function for hisat2
        
def hisat2_align(read1,read2,threads,lib,ref_name,prefix):
    #check if sam is already exsit, if not bwa with fastq
    if os.path.exists('%s.sam'%(prefix)) == False:
        if (lib=="RF") or  (lib=="FR"): 
            cmd = ' %s/hisat2 -x %s  -1 %s -2 %s -p %s --rna-strandness %s --fr -S  %s.sam  2>%s.hisat2.log  '%(hisat2_path,ref_name.replace(".fa",""),read1,read2,threads,lib,prefix,prefix)
        if lib=="unstranded":
            cmd = ' %s/hisat2 -x %s  -1 %s -2 %s -p %s --fr -S  %s.sam  2>%s.hisat2.log  '%(hisat2_path,ref_name.replace(".fa",""),read1,read2,threads,prefix,prefix)
        print(cmd)
        subprocess.run(cmd,shell=True,check=True)

    else: 
        print("output sam for "+prefix+" already exists")

def sam_to_bam(prefix):
    if os.path.exists('%s.bam'%(prefix)) == False:
        cmd1 = '%s view -uS  %s.sam | %s sort -@ 5 -o  %s.bam '%(samtools,prefix,samtools,prefix)
        print(cmd1)
        subprocess.run(cmd1,shell=True,check=True)
        cmd2 = '%s index %s.bam '%(samtools,prefix)
        print(cmd2)
        subprocess.run(cmd2,shell=True,check=True)
        if os.path.exists('%s.bam'%(prefix)) == True:
           os.remove(prefix+".sam")
    else: 
        print(" output hisat2  bam file for "+prefix+" already exists")

if __name__ == "__main__":

      hisat2_align(pe1,pe2,thread,library,genome,outname)
      sam_to_bam(outname)


