import argparse
import glob
import subprocess
import shutil
import os
import pathlib
import configparser
import sys
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool

# Read the command line arguments.
parser = argparse.ArgumentParser(description=" RSeqQC :geneBody_coverage,read_distribution,RPKM_saturation")
parser.add_argument("-b1", "--inputfile1", help="Input bam file")
parser.add_argument("-b2", "--inputfile2", help="Input genome bed file")
parser.add_argument("-p", "--output_prefix", help="output file's prefix,eg: Sample_A")
parser.add_argument("-o", "--outputDirectory", help="Output directory ")
args = parser.parse_args()

# Process the command line arguments.
bamfile = os.path.abspath(args.inputfile1)
bedfile = os.path.abspath(args.inputfile2)
prefix = args.output_prefix
outputDirectory = os.path.abspath(args.outputDirectory)

if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)


# define function for RseqQC
def geneBody_coverage(bed, bam, output, prefix):
    files = pathlib.Path(os.path.join(output, prefix + ".coverage.geneBodyCoverage.txt"))
    if files.exists() and len(files.read_text()) > 0:
        print(output + "/" + prefix + ".coverage.geneBodyCoverage.txt" + "is exsited,skip!")
    else:
        command1 = "geneBody_coverage.py -r {}  -i {} -o {}/{}.coverage".format(bed, bam, output, prefix)
        print(command1)
        subprocess.run(command1, shell=True, check=True)


def read_distribution(bed, bam, output, prefix):
    files = pathlib.Path(os.path.join(output, prefix + ".read_distribution.txt"))
    if files.exists() and len(files.read_text()) > 0:
        print(output + "/" + prefix + ".read_distribution.txt" + "is exsited,skip!")
    else:
        command2 = "read_distribution.py -r {}  -i {} > {}/{}.read_distribution.txt".format(bed, bam, output, prefix)
        print(command2)
        subprocess.run(command2, shell=True, check=True)
        odata = pd.read_table("./" + prefix + ".read_distribution.txt", header=None, sep='\s+', skiprows=5, nrows=10,
                              names=["Total_bases", "Tag_count", "Tags/Kb"])
        df = DataFrame(odata[["Tag_count"]])
        plt.switch_backend('agg')
        df.plot(kind='bar',
                color=['#9BCD9B', '#9BCD9B', '#9BCD9B', '#C1FFC1', 'lightskyblue', 'lightskyblue', 'lightskyblue',
                       'lightcoral', 'lightcoral', 'lightcoral'], legend=None, figsize=(8, 8), fontsize=8, width=0.7)
        plt.ylabel('Tag_count', fontsize=10)
        plt.xlabel('Group', fontsize=10)
        plt.title('read_distribution(' + prefix + ")", fontsize=12)
        plt.show()
        plt.savefig(prefix + ".read_distribution.pdf", bbox_inches='tight')
        plt.savefig(prefix + ".read_distribution.png", dpi=600)


def RPKM_saturation(bed, bam, output, prefix):
    files = pathlib.Path(os.path.join(output, prefix + ".saturation.eRPKM.xls"))
    if files.exists() and len(files.read_text()) > 0:
        print(output + "/" + prefix + ".saturation.pdf" + "is exsited,skip!")
    else:
        command3 = "RPKM_saturation.py   -r {}  -d '1++,1--,2+-,2-+'  -i {} -o {}/{}".format(bed, bam,
                                                                                             output, prefix)
        print(command3)
        subprocess.run(command3, shell=True,check=True)
        command4 = "{} {}/run_Saturation.r -i {}/{}.eRPKM.xls -o {} -s {}".format(Rscript, sys.path[0], output, prefix,
                                                                                  output, prefix)
        print(command4)
        subprocess.run(command4, shell=True,check=True)
        os.remove("{}/{}.saturation.pdf ".format(output, prefix))
        os.remove("{}/{}.saturation.r ".format(output, prefix))


def bam_stat(bam, output, prefix):
    files = pathlib.Path(os.path.join(output, prefix + ".bam_stat.out"))
    if files.exists() and len(files.read_text()) > 0:
        print(output + "/" + prefix + ".bam_stat.out" + " is exsited,skip!")
    else:
        command4 = "bam_stat.py -i {}  &> {}/{}.bam_stat.out ".format(bam, output, prefix)
        print(command4)
        subprocess.run(command4, shell=True, check=True)


# main functions
def main():
    pool = ThreadPool(4)
    pool.apply_async(geneBody_coverage, (bedfile, bamfile, outputDirectory, prefix))
    pool.apply_async(read_distribution, (bedfile, bamfile, outputDirectory, prefix))
    pool.apply_async(RPKM_saturation, (bedfile, bamfile, outputDirectory, prefix))
    pool.apply_async(bam_stat, (bamfile, outputDirectory, prefix))
    pool.close()
    pool.join()


# exe
if __name__ == "__main__":
    main()
