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
parser = argparse.ArgumentParser(description="Generates trim_galore scripts.")
parser.add_argument("-i", "--inputDirectory", help="Input directory with FASTQ files.(force)")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with trimmed FASTQ files. (default:./clean_data) ", default="./clean_data")
parser.add_argument("-c", "--configfile", help="The config file. (default:./config_files/config.ini) ",
                    default="./config_files/config.ini")
parser.add_argument("-f", "--fastqcoutputDirectory",
                    help="Output directory with FASTQC report. (default:./Quality_control/cleandata_fastqc)",
                    default="./Quality_control/cleandata_fastqc")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default: all) ",
                    choices=["all", "fat", "big"], default="all")
args = parser.parse_args()

# Process the command line arguments.

inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
fastqcoutputDirectory = os.path.abspath(args.fastqcoutputDirectory)
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
seqkit = config.get("software", "seqkit")
python3 = config.get("software", "python3")
itools = config.get("software", "itools")
java = config.get("software", "java")
Trimmomatic = config.get("software", "Trimmomatic")
Trimmomatic_adapter = config.get("software", "Trimmomatic_adapter")
fastqc = config.get("software", "fastqc")
sample_num = config.get("par", "sample_num")
clean_base = config.get("parameter", "clean_base")
email = config.get("parameter", "email")
thread_number = config.get("parameter", "thread_number")
read_length = config.get("parameter", "read_length")
# min_length = config.get("parameter", "min_length")
program = config.get("wkdir", "program")
sample_list = config.items("raw_data")

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
if not os.path.exists(fastqcoutputDirectory):
    os.makedirs(fastqcoutputDirectory)

# Write the script(s)

if not os.path.exists(program + "/script/clean_data_trimmomatic"):
    os.makedirs(program + "/script/clean_data_trimmomatic")
os.chdir(program + "/script/clean_data_trimmomatic")

if os.path.isfile("a.trimmomatic.sh"):
    os.remove("a.trimmomatic.sh")
if os.path.isfile("b.data_stat.sh"):
    os.remove("b.data_stat.sh")
if os.path.isfile("c.filter_ratio_stat.sh"):
    os.remove("c.filter_ratio_stat.sh")
if os.path.isfile("d.md5sum.sh"):
    os.remove("d.md5sum.sh")
# Cycle through all the samples, 2 by 2.

for i in range(0, len(sample_list), 1):
    sample_name = sample_list[i][0]
    # trimmomatic
    scriptName = 'a.trimmomatic' + '.sh'
    script = open(scriptName, 'a')
    script.write(java + " -jar  -Xmx40g " + " ")
    script.write(Trimmomatic + " ")
    script.write("PE " + " ")
    script.write("-threads 10 " + " ")
    script.write(os.path.join(inputDirectory, sample_name, sample_name + ".R1.fq.gz") + " ")
    script.write(os.path.join(inputDirectory, sample_name, sample_name + ".R2.fq.gz") + " ")
    script.write(os.path.join(outputDirectory, sample_name + ".R1.fq.gz") + " ")
    script.write(os.path.join(outputDirectory, sample_name + ".R1_un.fq.gz") + " ")
    script.write(os.path.join(outputDirectory, sample_name + ".R2.fq.gz") + "  ")
    script.write(os.path.join(outputDirectory, sample_name + ".R2_un.fq.gz") + "  ")
    script.write("CROP:" + read_length + " ")
#    script.write("ILLUMINACLIP:" + Trimmomatic_adapter + ":2:30:10:1:true" + " ")
    script.write("ILLUMINACLIP:" + Trimmomatic_adapter + ":2:30:10:8:true" + " ")
    script.write("LEADING:3" + " ")
    script.write("TRAILING:3 " + " ")
    script.write("SLIDINGWINDOW:4:15" + " ")
    script.write("MINLEN:" + "50" + " ")
    script.write("&& rm " + os.path.join(outputDirectory, sample_name + ".R1_un.fq.gz") + " ")
    script.write("&& rm " + os.path.join(outputDirectory, sample_name + ".R2_un.fq.gz") + "\n")

    # clean_data_stat
    scriptName = 'b.data_stat' + '.sh'
    script = open(scriptName, 'a')
    script.write(perl_path + " " + sys.path[0] + "/fastp2summary.pl " + " ")
    script.write(" -r1 " + os.path.join(inputDirectory, sample_name, sample_name + ".R1.fq.gz") + " ")
    script.write(" -r2 " + os.path.join(inputDirectory, sample_name, sample_name + ".R2.fq.gz") + " ")
    script.write(" -c1 " + os.path.join(outputDirectory, sample_name + ".R1.fq.gz") + " ")
    script.write(" -c2 " + os.path.join(outputDirectory, sample_name + ".R2.fq.gz") + " ")
    script.write(" -s " + sample_name + " ")
    script.write(" -o " + fastqcoutputDirectory + "\n")
    script.close()

##stat filter ratio 
scriptName = 'c.filter_ratio_stat' + '.sh'
script = open(scriptName, 'a')
script.write(perl_path + " " + sys.path[0] + "/" + "merge_fastp_summary.pl" + " ")
script.write("-i " + fastqcoutputDirectory + " ")
script.write("-o " + fastqcoutputDirectory + "/.. " + "\n")

##md5sum
scriptName = 'd.md5sum' + '.sh'
script.write("cd " + outputDirectory + " && ")
script.write("md5sum  " + " *.fq.gz  > " + outputDirectory + "/clean_data_md5.xls " + " && ")
script.write("cp " + sys.path[0] + "/md5* " + outputDirectory + "\n")
script.close()

scripts = glob.glob("*.sh")
scripts.sort()
for i in range(0, len(scripts), 1):
    command = "sh " + scripts[i]
    subprocess.run(command, shell=True,check=True)

# email the check results
scriptName = 'e.clean_data_check_email' + '.cmd'
script = open(scriptName, 'w')
script.write(python3 + " " + sys.path[0] + "/" + "cleandata_check.py " + " ")
script.write("-f " + fastqcoutputDirectory + "/../report_Filtered_data_stat.xls ")
script.write("-c " + outputDirectory + " ")
script.write("-t " + clean_base + "G" + " ")
script.write("-e " + email + " ")
script.write("-o " + outputDirectory + "/../log/email_clean_data_infor.txt  ")
script.write("-s " + sample_num + "\n")
script.close()
filter_check = "sh e.clean_data_check_email.cmd >email_filter.log "
subprocess.run(filter_check, shell=True,check=True)
