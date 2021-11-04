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
parser = argparse.ArgumentParser(description="Generates hisat2 scripts.")
parser.add_argument("-i", "--inputDirectory", help="Input directory with clean data files.(default:./clean_data)",
                    default="./clean_data")
parser.add_argument("-g", "--genomefasta", help="Input genome fasta for hisat2.(default:./genome/genome.fa)",
                    default="./genome/genome.fa")
parser.add_argument("-o", "--outputDirectory", help="Output directory with hisat2 results.(default:./alignment)",
                    default="./alignment")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config_files/config.ini)",
                    default="./config_files/config.ini")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all) ",
                    choices=["all", "fat", "big"], default="all")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
genome = os.path.abspath(args.genomefasta)
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
python3 = config.get("software", "python3")
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
thread_number = config.get("parameter", "thread_number")
program = config.get("wkdir", "program")
sample_list = config.items("raw_data")

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Write the script(s)
if not os.path.exists(program + "/script/hisat2_alignment"):
    os.makedirs(program + "/script/hisat2_alignment")
os.chdir(program + "/script/hisat2_alignment")

if os.path.isfile("a.hisat2_index.sh"):
    os.remove("a.hisat2_index.sh")
if os.path.isfile("b.hisat2_align.sh"):
    os.remove("b.hisat2_align.sh")

# Index for genome
scriptName = 'a.hisat2_index' + '.sh'
script = open(scriptName, 'a')
script.write(python3 + " " + sys.path[0] + "/hisat2_index.py" + " ")
script.write("-g " + genome + " ")
script.write("-n 8 " + "\n")
script.close()

# Cycle through all the samples, 2 by 2.

for i in range(0, len(sample_list), 1):
    sample_name = sample_list[i][0]
    # Create script file.
    scriptName = 'b.hisat2_align' + '.sh'
    script = open(scriptName, 'a')
    script.write(python3 + " " + sys.path[0] + "/hisat2_align.py" + " ")
    script.write("-pe1 " + os.path.join(inputDirectory, sample_name+".R1.fq.gz") + " ")
    script.write("-pe2 " + os.path.join(inputDirectory, sample_name+".R2.fq.gz") + " ")
    script.write("-g " + genome + " ")
    script.write("-n 8 " + " ")
    script.write("-l " + " RF ")
    script.write("-o " + os.path.join(outputDirectory, sample_name ) + " ")
    script.write("-p " + sample_name + " \n")
    script.close()

scripts = glob.glob("*.sh")
scripts.sort()
for i in range(0, len(scripts), 1):
    command = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + " --maxproc " + thread_number + " " + \
              scripts[i]
    subprocess.run(command, shell=True,check=True)
