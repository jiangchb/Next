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
parser = argparse.ArgumentParser(description="Generates RseqQC scripts.")
parser.add_argument("-i", "--inputDirectory",
                    help="Input directory with alignment files(bam from hisat2 or tophat2). (default:./alignment) ",
                    default="./alignment")
parser.add_argument("-b", "--inputbed", help="Input bed files. (default:./genome/exon.bed)",
                    default="./genome/exon.bed")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with RseqQC  results. (default:./Quality_control/RseqQC) ",
                    default="./Quality_control/RseqQC")
parser.add_argument("-c", "--configfile", help="The config file.  (default:./config_files/config.ini)",
                    default="./config_files/config.ini")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all) ",
                    choices=["all", "fat", "big"], default="all")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
genebed = os.path.abspath(args.inputbed)
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
thread_number = config.get("parameter", "thread_number")
python2_path = config.get("software", "python2_path")
R_path = config.get("software", "R")
cuffcompare = config.get("software", "cuffcompare")
sample_num = config.get("par", "sample_num")
email = config.get("parameter", "email")
gemome_alignment_thresh = config.get("parameter", "gemome_alignment_thresh")
program = config.get("wkdir", "program")
sample_list = config.items("raw_data")

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
# Write the script(s)
if not os.path.exists(program + "/script/bam_RseqQC"):
    os.makedirs(program + "/script/bam_RseqQC")
os.chdir(program + "/script/bam_RseqQC")
if os.path.isfile("a.RseqQC.sh"):
    os.remove("a.RseqQC.sh")
if os.path.isfile("b.merged.bam_out.sh"):
    os.remove("b.merged.bam_out.sh")
if os.path.isfile("c.genome_alignment_check.cmd"):
    os.remove("c.genome_alignment_check.cmd")

# Cycle through all the samples, 1 by 1.

for i in range(0, len(sample_list), 1):
    sample_name = sample_list[i][0]
    # Create script file.
    scriptName = 'a.RseqQC' + '.sh'
    script = open(scriptName, 'a')
    script.write("export PATH=" + os.path.dirname(R_path) + ":$PATH " + " && ")
    script.write("export PATH=" + python2_path + ":$PATH " + " && ")
    script.write("export LD_LIBRARY_PATH=/home/fanyucai/software/python/Python-v2.7.9/lib:$LD_LIBRARY_PATH" + " && ")
    script.write(python3 + " " + sys.path[0] + "/" + "bam_RseqQC.py ")
    script.write("-b1 " + os.path.join(inputDirectory, sample_name, (sample_name + ".bam")) + " ")
    script.write("-b2 " + genebed + " ")
    script.write("-p " + sample_name + " ")
    script.write("-o " + os.path.join(outputDirectory, sample_name) + "\n")
    script.close()

scriptName = 'b.merged.bam_out' + '.sh'
script = open(scriptName, 'w')
script.write(python3 + " " + sys.path[0] + "/" + "bam_stat_merge.py ")
script.write("-i " + outputDirectory + " ")
script.write("-n " + os.path.join(outputDirectory, "genome_mapping_stat.xls") + "\n")
script.close()

# email the check results
scriptName = 'c.genome_alignment_check' + '.cmd'
script = open(scriptName, 'a')
script.write(python3 + " " + sys.path[0] + "/" + "genome_alignment_check.py " + " ")
script.write("-f " + os.path.join(outputDirectory, "genome_mapping_stat.xls") + " ")
script.write("-t " + gemome_alignment_thresh + " ")
script.write("-e " + email + " ")
script.write("-o " + outputDirectory + "/../../log/email_genome_alignment_infor.txt  ")
script.write("-s " + sample_num + "\n")
script.close()

scripts = glob.glob("*.sh")
scripts.sort()
for i in range(0, len(scripts), 1):
    command = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + " --maxproc " + thread_number + " " + \
              scripts[i]
    subprocess.run(command, shell=True,check=True)
check_cmd = "sh c.genome_alignment_check.cmd  >email_genome_alignment_check.log "
subprocess.run(check_cmd, shell=True,check=True)
