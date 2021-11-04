#!/usr/bin/env python3

# Version 1.1
# Author: xiufeng yang
# Date: 10/10/2017
# Email: yxfhenu@163.com

import argparse
import glob
import os
import subprocess
import shutil
import sys
import configparser

# Read the command line arguments.
parser = argparse.ArgumentParser(description="script for cpc  analysis.")
parser.add_argument("-i", "--fasta", help="Input transcript fasta (force).")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config_files/config.ini)",
                    default="./config_files/config.ini")
parser.add_argument("-o", "--outputDirectory", help="Output directory with cpc results")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all)",
                    choices=["all", "fat", "big"], default="all")
args = parser.parse_args()

# Process the command line arguments.

fasta = os.path.abspath(args.fasta)
number = args.splitnumber
outputDirectory = os.path.abspath(args.outputDirectory)
configfile = args.configfile
qsubquene = args.qsub_queue_server

# Get the software path and Parameter
config = configparser.ConfigParser()
config.read(configfile, encoding="utf-8")
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
python2 = config.get("software", "python2")
cpc2 = config.get("software", "cpc2")

# define function for get lncRNA canditate
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)

if os.path.isfile("a.cpc_run.sh"):
    os.remove("a.cpc_run.sh")

########cpc###################
scriptName = 'a.cpc_run' + '.sh'
script = open(scriptName, 'a')
script.write('{python2}   {cpc2}  -i {fasta} -o {results} '.format(python2=python2, cpc2=cpc2, fasta=fasta,
                                                                   results=os.path.join(outputDirectory,
                                                                                        filename + ".xls")))
script.close()
######run#####################
scripts = glob.glob("*.sh")
scripts.sort()
for i in range(1, len(scripts), 1):
    command = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " --maxproc " + number + " " + scripts[i]
    subprocess.run(command, shell=True, check=True)
