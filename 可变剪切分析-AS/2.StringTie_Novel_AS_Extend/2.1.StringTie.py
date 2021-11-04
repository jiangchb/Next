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
parser = argparse.ArgumentParser(description="Generates StringTie scripts.")
parser.add_argument("-i", "--inputDirectory",
                    help="Input directory with alignment files(bam from hisat2 or tophat2).(default:./alignment)",
                    default="./alignment")
parser.add_argument("-g", "--genomefasta", help="Reference genome fasta.(default:./genome/genome.fa)",
                    default="./genome/genome.fa ")
parser.add_argument("-f", "--referencegff",
                    help="Reference mRNA annotatation gff .(default:genome/transcriptome_exon.gff )",
                    default="./genome/transcriptome_exon.gff")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config_files/config.ini)",
                    default="./config_files/config.ini")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with StringTie results.(default:./reconstruction)",
                    default="./reconstruction")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all)",
                    choices=["all", "fat", "big"], default="all")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
outputDirectory = os.path.abspath(args.outputDirectory)
referencegff = os.path.abspath(args.referencegff)
genomefasta = os.path.abspath(args.genomefasta)
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
sample_list = config.items("raw_data")
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
gffread = config.get("software", "gffread")
stringtie = config.get("software", "stringtie")
cuffcompare = config.get("software", "cuffcompare")
thread_number = config.get("parameter", "thread_number")
program = config.get("wkdir", "program")

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)

# Write the script(s)
if not os.path.exists(program + "/script/StringTie"):
    os.makedirs(program + "/script/StringTie")
os.chdir(program + "/script/StringTie")
if os.path.isfile("a.stringtie.sh"):
    os.remove("a.stringtie.sh")
if os.path.isfile("b.merged.compare.sh"):
    os.remove("b.merged.compare.sh")
if os.path.isfile(os.path.join(outputDirectory,"mergelist.txt")):
    os.remove(os.path.join(outputDirectory,"mergelist.txt"))

# Cycle through all the samples, 1 by 1.
for i in range(0, len(sample_list), 1):
    sample_name = sample_list[i][0]
    if not os.path.exists(os.path.join(outputDirectory, sample_name)):
        os.makedirs(os.path.join(outputDirectory, sample_name))
    # Create script file.
    scriptName = 'a.stringtie' + '.sh'
    script = open(scriptName, 'a')
    script.write(stringtie + " ")
    script.write("-p 8 " + " ")
    script.write("--rf" + " ")
    script.write("-l " + sample_name + " ")
    script.write("-o " + os.path.join(outputDirectory, sample_name, (sample_name + ".gtf")) + "  ")
    script.write(os.path.join(inputDirectory, sample_name, (sample_name + ".bam")) + " \n")
    script.close()
    scriptName = os.path.join( outputDirectory, 'mergelist.txt')
    script = open(scriptName, 'a')
    script.write(os.path.join(outputDirectory, sample_name, (sample_name + ".gtf")) + " \n")
    script.close()

scriptName = 'b.merged.compare' + '.sh'
script = open(scriptName, 'a')
script.write(gffread + " -T -o " + os.path.join(outputDirectory, "anno.gtf ") + " " + referencegff + "&&")
script.write("cut -d\';\' -f1 " + os.path.join(outputDirectory, "anno.gtf ") + " >" + os.path.join(outputDirectory,
                                                                                                   "tmp.gtf ") + "&&")
script.write(
    "mv " + os.path.join(outputDirectory, "tmp.gtf ") + " " + os.path.join(outputDirectory, "anno.gtf ") + "&&")
script.write(stringtie + " ")
script.write("--merge " + " ")
script.write("-p 8 " + " ")
script.write("-o " + os.path.join(outputDirectory, "stringtie_merged.gtf") + "  ")
script.write(os.path.join(outputDirectory, "mergelist.txt") + " && ")
script.write("cd " + outputDirectory + " && ")
script.write(cuffcompare + " ")
script.write("-s  " + genomefasta + " ")
script.write("-r " + os.path.join(outputDirectory, "anno.gtf ") + " ")
script.write(os.path.join(outputDirectory, "stringtie_merged.gtf") + "\n")
script.close()

scripts = glob.glob("*.sh")
scripts.sort()
for i in range(0, len(scripts), 1):
    command = perl_path + " " + qsub_pbs + " --queue " + "fat" + " " + " --maxproc " + "8" + " " + scripts[i]
    subprocess.run(command, shell=True,check=True)
