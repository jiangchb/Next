# Version 1.1
# Author: xiufeng yang
# Date: 10/04/2017
# Email: yxfhenu@163.com
# add gene by wsb 2019.08.01

import argparse
import os
import shlex
import subprocess
import shutil
import sys
import glob

# Define software path 

ASprofile_path = "/home/fanyucai/software/ASprofile/ASprofile.b-1.0.4"
bedtools = "/home/fanyucai/software/bedtools/bedtools2/bin/bedtools"
#bedtools = "/home/xtsheng/miniconda2/bin/bedtools"
Rscript = "/home/fanyucai/software/R/R-v3.4.0/bin/Rscript"

# Read the command line arguments.
parser = argparse.ArgumentParser(description="run Asprofile for each sample.")
parser.add_argument("-i", "--inputDirectory",
                    help="Input directory with reconstruction gtf. (default: ./reconstruction)",
                    default="./reconstruction")
parser.add_argument("-g", "--refgenome", help="Input genome.fa (force).")
parser.add_argument("-b", "--refbed", help="Input exon.bed (force).")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with Advanced analysis results.[default:./Advanced_analysis]",
                    default="./Advanced_analysis")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
genome = os.path.abspath(args.refgenome)
bed = os.path.abspath(args.refbed)
outputDirectory = os.path.abspath(args.outputDirectory)

# Creat outputDirectory if it is not exist
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)


# define function for CIRI2

def Asprofile_run(input, genome, bed):
    if not os.path.exists("./ASprofile"):
        os.makedirs("./ASprofile")
    os.chdir("./ASprofile")
    files = glob.glob(input + "/*/*.gtf")
    files.sort()
    for i in range(0, len(files), 1):
        filename = os.path.basename(os.path.dirname(files[i]))
        cmd0 = "export LD_LIBRARY_PATH=/home/fanyucai/software/zlib/zlib-v1.2.11/lib:$LD_LIBRARY_PATH"

        cmd1 = '%s/extract-as  %s  %s | cut -f 1,2,3,4,5,6,8,9,15,16 > %s_AS_result_tmp.txt' % (
            ASprofile_path, files[i], genome, filename)
        cmd2 = 'awk -F \"\\t\"  -v OFS=\"\\t\"  \'{if(NR>1) print $4,$5-1,$6, $1,"0",$7}\'  %s_AS_result_tmp.txt > tmp.bed' % (
            filename)
        cmd3 = '%s intersect -wa -wb -a tmp.bed -b %s | cut -f 4,10 | sort -u |awk \'BEGIN{FS=OFS=\"\\t\"}{a[$1]=a[$1]\";\"$2} END{for(j in a) print j,a[j] }\' |sed \'s:\\t;:\\t:\' > intersect.tmp' % (
            bedtools, bed)
        cmd4 = 'awk -F \"\\t\"  -v OFS=\"\\t\"  \'NR==FNR{a[$1]=$2; next} {print $0,a[$1]}\'  intersect.tmp %s_AS_result_tmp.txt | awk -F \"\\t\" \'{if(NR==1) $11=\"ref_id\"; print $0}\' OFS=\"\\t\" |cut -f1-8,10-11 > %s_AS_result.xls' % (
            filename, filename)
        cmd5 = 'rm  %s_AS_result_tmp.txt tmp.bed  intersect.tmp' % (filename)
        print(cmd0)
        subprocess.run(cmd0, shell=True, check=True)
        print(cmd1)
        subprocess.run(cmd1, shell=True, check=True)
        print(cmd2)
        subprocess.run(cmd2, shell=True, check=True)
        print(cmd3)
        subprocess.run(cmd3, shell=True, check=True)
        print(cmd4)
        subprocess.run(cmd4, shell=True, check=True)
        print(cmd5)
        subprocess.run(cmd5, shell=True, check=True)
    for i in range(0, len(files), 1):
        filename = os.path.basename(os.path.dirname(files[i]))
        cmd6 = '%s  %s/ASprofile_plot.r -i %s_AS_result.xls -o . -s %s ' % (Rscript, sys.path[0], filename, filename)
        print(cmd6)
        subprocess.run(cmd6, shell=True, check=True)
    os.chdir(outputDirectory)

def Asprofile_add_gene():
    mRNA_2_gene = os.path.join(os.path.dirname(genome), "mRNA_2_gene.xls")
    if mRNA_2_gene :
        cmd7 = 'sh %s/add_enrichment.sh %s %s %s ' % (sys.path[0], sys.path[0], os.path.join(outputDirectory, "ASprofile"), mRNA_2_gene)
        subprocess.run(cmd7, shell=True, check=True)
    else:
        print ("There is no mRNA_2_gene.xls, please Check out !")



if __name__ == "__main__":
    Asprofile_run(inputDirectory, genome, bed)
    Asprofile_add_gene()