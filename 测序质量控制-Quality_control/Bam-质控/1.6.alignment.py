#!/usr/bin/env python3

# Version 1.1
# Author Xiufeng Yang
# Date: 08/18/2017

import argparse
import glob
import os
import subprocess
import shutil
import sys
import configparser

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates alignment picture draw scripts.")
parser.add_argument("-i", "--inputDirectory",
                    help="Input directory with alignment files(bam from hisat2 or tophat2). (default:./alignment) ",
                    default="./alignment")
parser.add_argument("-b", "--inputbed", help="Input bed files. (default:./genome/exon.bed)",
                    default="./genome/exon.bed")
parser.add_argument("-f", "--inputgff", help="Input gff files. (default:./genome/transcriptome_exon.gff)",
                    default="./genome/transcriptome_exon.gff")
parser.add_argument("-g", "--inputgenome", help="Input genome files. (default:./genome/genome.fa)",
                    default="./genome/genome.fa")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with RseqQC  results. (default:./Quality_control/genome_gene_coverage_plot) ",
                    default="./Quality_control/genome_gene_coverage_plot")
parser.add_argument("-c", "--configfile", help="The config file.  (default:./config_files/config.ini)",
                    default="./config_files/config.ini")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all) ",
                    choices=["all", "fat", "big"], default="all")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
bed = os.path.abspath(args.inputbed)
gff = os.path.abspath(args.inputgff)
genome = os.path.abspath(args.inputgenome)
outputDirectory = os.path.abspath(args.outputDirectory)
configfile = args.configfile
qsubquene = args.qsub_queue_server

# Get the software path and Parameter
class myconf(configparser.ConfigParser):
    def __init__(self, defaults=None):
        configparser.ConfigParser.__init__(self, defaults=defaults)
    def optionxform(self, optionstr):
        return optionstr
config = myconf()
config.read(configfile, encoding="utf-8")
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
python3 = config.get("software", "python3")
bedtools = config.get("software", "bedtools")
thread_number = config.get("parameter", "thread_number")
program = config.get("wkdir", "program")
sample_list = config.items("raw_data")
R_path = config.get("software", "R")

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Write the script(s)
if not os.path.exists(program + "/script/genome_gene_coverage_plot"):
    os.makedirs(program + "/script/genome_gene_coverage_plot")
os.chdir(program + "/script/genome_gene_coverage_plot")

if os.path.isfile("a.genome_gene_coverage_plot.sh"):
    os.remove("a.genome_gene_coverage_plot.sh")
# Cycle through all the samples, 1 by 1.

for i in range(0, len(sample_list), 1):
    sample_name = sample_list[i][0]
    if not os.path.exists(os.path.join(outputDirectory, sample_name )):
        os.makedirs(os.path.join(outputDirectory, sample_name ))
    # Create script file.
    scriptName = 'a.genome_gene_coverage_plot' + '.sh'
    script = open(scriptName, 'a')
    script.write("export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:$LD_LIBRARY_PATH" + " && ")    
    script.write("export LD_LIBRARY_PATH=/home/fanyucai/software/zlib/zlib-v1.2.11/lib:$LD_LIBRARY_PATH" + "&&")
    script.write("export  PATH=/home/fanyucai/software/bedtools/bedtools2/bin:$PATH" + "&&")
    script.write("export PATH=" + R_path + "/..:$PATH " + " && ")
    script.write(python3 + " " + sys.path[0] + "/" + "alignment_coverage_plot.py ")
    script.write("-b1 " + os.path.join(inputDirectory, sample_name , (sample_name  + ".bam")) + " ")
    script.write("-b2 " + bed + " ")
    script.write("-b3 " + gff + " ")
    script.write("-b4 " + genome + " ")
    script.write("-p " + sample_name + " ")
    script.write("-o " + os.path.join(outputDirectory, sample_name ) + "\n")
    script.close()

scripts = glob.glob("*.sh")
scripts.sort()
for i in range(0, len(scripts), 1):
    command = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + " --maxproc " + thread_number + " " + \
              scripts[i]
    subprocess.run(command, shell=True,check=True)
