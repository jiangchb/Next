# Version 1.1
# Author: Congjia Chen
# Date: 19/08/2021
import argparse
import glob
import os
import subprocess
import shutil
import sys
import configparser
import pathlib

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates htseq-counts script.")
parser.add_argument("-b", "--input_BAM_filepath", help="Input directory with genome alignment files.",
                    default="./1.hisat_result")
parser.add_argument("-f", "--gene_gtf", help="Input gene annotation gtf file.[default:./genome/gene.gtf]",
                    default="./genome/gene.gtf")
#parser.add_argument("-g", "--ref_genome", help="Input genome fasta file.[default:./genome/genome.fa]",
#                    default="./genome/genome.fa")
#parser.add_argument("-p", "--prefix", help="Output file prefix name ")
parser.add_argument("-o", "--outputDirectory", help="Output directory with fpkm_counts files.[default:./Quantification/fpkm_counts]",
                    default="./Quantification/fpkm_counts")
args = parser.parse_args()

# Process the command line arguments.

bam = os.path.abspath(args.input_BAM_filepath)
gtf = os.path.abspath(args.gene_gtf)
#genome = os.path.abspath(args.ref_genome)
#prefix = args.prefix
outputDirectory = os.path.abspath(args.outputDirectory)

# Define software path 
# huawei cloud
#htseq = "/home/fanyucai/software/samtools/samtools-v1.4/bin/htseq-count"
#cufflinks = "/data/software/Cufflinks/2.2.1/cufflinks"
#bamtools = "/data/software/IND_wRNA/miniconda3/envs/bamtools/bin/bamtools"

htseq = "/home/fanyucai/software/python/Python-v2.7.9/bin/htseq-count"
cufflinks = "/home/fanyucai/software/Cufflinks/cufflinks-2.2.1.Linux_x86_64/cufflinks"
bamtools = "/home/fanyucai/software/bamtools/bin/bamtools"

# creat outputDirectory if it is not exist
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)

# define function
def htseq_counts(bam_file, gtf, output, prefix):
    # check if counts files is already exsit, if not htseq-count with bam
    files = pathlib.Path('%s.counts.txt' % (prefix))
    if files.exists() and os.path.getsize(files) > 0:
        print(" counts file for " + prefix + " is already existed,please check!Skipping!")
    # run htseq-count
    else:
        #if not os.path.exists('%s.bam' % (prefix)):
            #cmd1 = '%s sort -in %s -out %s.bam ' % (bamtools, bam, prefix)
            #cmd1 = 'ln -s %s %s.bam ' % (bam, prefix)
            #print(cmd1)
            #subprocess.run(cmd1, shell=True, check=True)
        cmd2 = '%s -i gene_id  -f bam -s no -r name  %s %s > %s/%s.counts.txt' % (htseq, bam_file, gtf, output, prefix)
        print(cmd2)
        subprocess.run(cmd2, shell=True, check=True)
#    if os.path.exists('%s.counts.txt' % (prefix)):
#        os.remove('%s.bam'%(prefix))


def cufflinks_fpkm(genome,gtf,bam,output,prefix):
    # check if fpkm files is already exsit, if not cufflinks with bam
    files = pathlib.Path('%s.fpkm.txt' % (prefix))
    if files.exists() and os.path.getsize(files) > 0:
        print(" fpkms file for " + prefix + " is already existed,please check!Skipping!")
    # run cufflinks
    else:
        cmd1 = '%s --no-update-check -o %s -u  --max-bundle-frags 2000000 -p 4   --library-type fr-firststrand -b %s -G %s %s' % (
        cufflinks, output, genome, gtf, bam)
        print(cmd1)
        subprocess.run(cmd1, shell=True, check=True)
        os.rename("genes.fpkm_tracking", prefix + ".fpkm.txt")
        os.remove("isoforms.fpkm_tracking")
        os.remove("transcripts.gtf")
        os.remove("skipped.gtf")

if __name__ == "__main__":
    #deal with the prefix
    Detect_path=bam+"/"+"*.bam"
    BAM_file_list=glob.glob(Detect_path)
    for i in BAM_file_list:
        sample_name_tmp=os.path.split(i)[1]
        prefix=os.path.splitext(sample_name_tmp)[0]
        htseq_counts(i, gtf,outputDirectory, prefix)
#    cufflinks_fpkm(genome,gtf,bam,outputDirectory,prefix)

