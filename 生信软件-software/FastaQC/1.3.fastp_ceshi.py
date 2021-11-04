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
parser.add_argument("-i", "--inputDirectory",
                    help="Input directory with FASTQ files. (default:./raw_data_after_rRNA_remove)",
                    default="./raw_data_after_rRNA_remove")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with trimmed FASTQ files. (default:./clean_data) ", default="./clean_data")
parser.add_argument("-c", "--configfile", help="The config file. (default:./rRNA_remove_RNAseq.conf) ",
                    default="./rRNA_remove_RNAseq.conf")
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
config = configparser.ConfigParser()
config.read(configfile)
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
seqkit = config.get("software", "seqkit")
python3 = config.get("software", "python3")
itools = config.get("software", "itools")
java = config.get("software", "java")
Trimmomatic = config.get("software", "Trimmomatic")
Trimmomatic_adapter = config.get("software", "Trimmomatic_adapter")
fastqc = config.get("software", "fastqc")
program = config.get("wkdir", "program")
read_length = config.get("parameter", "read_length")
# min_length = config.get("parameter", "min_length")
fastp = "/home/fanyucai/software/fastp/fastp-master/fastp"
# Check if the inputDirectory exists, and is a directory.
if not os.path.exists(args.inputDirectory):
    if args.inputDirectory == "./raw_data_after_rRNA_remove":
        sys.stderr.write("Error!\n")
        sys.stderr.write("The input directory (default) does not exist.\n")
        sys.stderr.write("inputDirectory (default): " + args.inputDirectory + "\n")
    else:
        sys.stderr.write("Error!\n")
        sys.stderr.write("The input directory does not exist.\n")
        sys.stderr.write("inputDirectory: " + args.inputDirectory + "\n")
    sys.exit(1)

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
if not os.path.exists(fastqcoutputDirectory):
    os.makedirs(fastqcoutputDirectory)

# Store the list of files with the extensions fastq or fastq.gz
files = glob.glob(inputDirectory + "/*/*.fq.gz") + glob.glob(inputDirectory + "/*/*.fastq.gz") + glob.glob(
    inputDirectory + "/*/*.fq") + glob.glob(inputDirectory + "/*/*.fastq")
files.sort()

# Write the script(s)
if not os.path.exists(program + "/script/fastp"):
    os.makedirs(program + "/script/fastp")
os.chdir(program + "/script/fastp")

if os.path.isfile("a.rawdata_stat.sh"):
    os.remove("a.rawdata_stat.sh")
if os.path.isfile("b.fastp.sh"):
    os.remove("b.fastp.sh")
if os.path.isfile("c.cleandata_stat.sh"):
    os.remove("c.cleandata_stat.sh")
if os.path.isfile("d.filter_ratio_stat.sh"):
    os.remove("d.filter_ratio_stat.sh")

# Cycle through all the samples, 2 by 2.

for i in range(0, len(files), 2):
    filename = os.path.basename(os.path.dirname(files[i]))
    fileR1 = os.path.basename(files[i])
    fileR2 = os.path.basename(files[i + 1])

    # Create script file.
    scriptName = 'a.rawdata_stat' + '.sh'
    script = open(scriptName, 'a')

    # raw_data_stat
    script.write(itools + " Fqtools stat" + " ")
    script.write(" -InFq " + os.path.join(inputDirectory, filename, fileR1) + " ")
    script.write(" -InFq " + os.path.join(inputDirectory, filename, fileR2) + " ")
    script.write(" -OutStat " + outputDirectory + "/" + filename + ".raw_data.stat.txt")
    script.write(" -CPU " + "2" + "\n")
    script.close()

    # trimmomatic
    scriptName = 'b.fastp' + '.sh'
    script = open(scriptName, 'a')
    script.write("export LD_LIBRARY_PATH=/home/fanyucai/software/zlib/zlib-v1.2.11/lib:$LD_LIBRARY_PATH" + "&&")
    script.write(fastp + " ")
    script.write("-l " + " 50  ")
    script.write("-w " + " 10 ")
    script.write("-i " + os.path.join(inputDirectory, filename, fileR1) + " ")
    script.write("-I " + os.path.join(inputDirectory, filename, fileR2) + " ")
    script.write("-o " + os.path.join(outputDirectory, filename) + ".R1.fq.gz  ")
    script.write("-O " + os.path.join(outputDirectory, filename) + ".R2.fq.gz  ")
    script.write("-j " + os.path.join(outputDirectory, filename) + ".json  ")
    script.write("-h " + os.path.join(outputDirectory, filename) + ".html ")

    # fastqc
    script.write("&& " + fastqc + " -o " + fastqcoutputDirectory + " ")
    script.write("--extract  " + " -f fastq " + "-t 5 " + " ")
    script.write(os.path.join(outputDirectory, filename) + ".R1.fq.gz" + " ")
    script.write(os.path.join(outputDirectory, filename) + ".R2.fq.gz" + "\n")
    script.close()

    # clean_data_stat
    scriptName = 'c.cleandata_stat' + '.sh'
    script = open(scriptName, 'a')
    script.write(itools + " Fqtools stat" + " ")
    script.write(" -InFq " + os.path.join(outputDirectory, filename + ".R1.fq.gz") + " ")
    script.write(" -InFq " + os.path.join(outputDirectory, filename + ".R2.fq.gz") + " ")
    script.write(" -OutStat " + outputDirectory + "/" + filename + ".clean_data.stat.txt")
    script.write(" -CPU " + " 2 " + "\n")
    script.close()

##stat filter ratio 
scriptName = 'd.filter_ratio_stat' + '.sh'
script = open(scriptName, 'a')
script.write(python3 + " " + sys.path[0] + "/" + "stat_filter_ratio.py " + " ")
script.write("-i " + outputDirectory + " ")
script.write("-o " + fastqcoutputDirectory + "/.. " + " ")
script.write("-n " + " Filtered_data_stat.txt " + " && ")
script.write("cd " + outputDirectory + " && ")
script.write("md5sum  " + " *.fq.gz > " + outputDirectory + "/clean_data_md5.txt " + " && ")
script.write("cp " + sys.path[0] + "/md5* " + outputDirectory + "\n")
script.close()

scripts = glob.glob("*.sh")
scripts.sort()
for i in range(0, len(scripts), 1):
    command = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + scripts[i]
    subprocess.run(command, shell=True,check=True)
