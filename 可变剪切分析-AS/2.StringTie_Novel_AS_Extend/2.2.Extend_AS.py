# Version 1.1
# Author: xiufeng yang
# Date: 10/11/2017
# Email: yxfhenu@163.com

import argparse
import os
import shlex
import subprocess
import shutil
import sys
import glob
import configparser
from concurrent.futures import ThreadPoolExecutor, as_completed

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates nover transcript and extend gene candidate .")
parser.add_argument("-i", "--inputDirectory",
                    help="Input directory with reconstruction gtf. (default: ./reconstruction)",
                    default="./reconstruction")
parser.add_argument("-g", "--refgenome", help="Input genome.fa (force).")
parser.add_argument("-f", "--refgff", help="Input transcriptome_exon.gff (force).")
parser.add_argument("-b", "--refbed", help="Input exon.bed (force).")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config_files/config.ini)",
                    default="./config_files/config.ini")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all)",
                    choices=["all", "fat", "big"], default="all")
parser.add_argument("-o", "--outputDirectory",
                    help="Output directory with Advanced analysis results.[default:./Advanced_analysis]",
                    default="./Advanced_analysis")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
gff = os.path.abspath(args.refgff)
genome = os.path.abspath(args.refgenome)
bed = os.path.abspath(args.refbed)
outputDirectory = os.path.abspath(args.outputDirectory)
qsubquene = args.qsub_queue_server
configfile = args.configfile

# Get the software path and Parameter
config = configparser.ConfigParser()
config.read(configfile, encoding="utf-8")
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
python3 = config.get("software", "python3")
thread_number = config.get("parameter", "thread_number")
program = config.get("wkdir", "program")

# Check if the inputDirectory exists, and is a directory.
if not os.path.exists(args.inputDirectory):
    if args.inputDirectory == "./reconstruction":
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
if not os.path.exists(program + "/script/Extend_AS"):
    os.makedirs(program + "/script/Extend_AS")
os.chdir(program + "/script/Extend_AS")

if os.path.exists("a.get_candidate_run.sh"):
    os.remove("a.get_candidate_run.sh")
if os.path.exists("b.pblat_as.sh"):
    os.remove("b.pblat_as.sh")

##print script for get_candidate
scriptName = 'a.get_candidate_run' + '.sh'
script = open(scriptName, 'a')
script.write(python3 + " " + sys.path[0] + "/get_novel_exend_candidate.py" + " ")
script.write("-i " + inputDirectory + " ")
script.write("-f " + gff + " ")
script.write("-g " + genome + " ")
script.write("-b " + bed + " ")
script.write("-o " + outputDirectory + "\n")
script.close()

########pblat###################    
scriptName = 'b.pblat_as' + '.sh'
script = open(scriptName, 'a')
script.write(python3 + " " + os.path.join(sys.path[0], "pblat_Asprofile_run.py") + " ")
script.write(" -a " + inputDirectory + " ")
script.write(" -e " + os.path.join(outputDirectory, "Extend_gene") + " ")
script.write(" -g " + genome + " ")
script.write(" -b " + bed + " ")
script.write(" -c " + configfile + " ")
script.write(" -o " + os.path.join(outputDirectory) + " ")
script.write(" -q  " + qsubquene + "\n")
script.close()


def pblat_run(files):
    f = open(files, 'r')
    lines = f.readlines()
    cmd1 = lines[0].strip('\n')
    subprocess.run(cmd1, shell=True, check=True)


if __name__ == "__main__":
    cmd1 = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + " --maxproc " + thread_number + " " + "a.get_candidate_run.sh"
    subprocess.run(cmd1, shell=True, check=True)
    pblat_run("b.pblat_as.sh")
