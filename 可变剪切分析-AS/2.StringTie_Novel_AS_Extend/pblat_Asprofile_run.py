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

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Run Asprofile and Pblat.")
parser.add_argument("-a", "--inputDirectory", help="Input directory with reconstruction gtf. (default: ./reconstruction)", default="./reconstruction")
parser.add_argument("-e", "--extenddirectory", help="Input directory with extend gene candidate. (default: ./Advanced_analysis/Extend_gene)", default="./Advanced_analysis/Extend_gene")
parser.add_argument("-g", "--refgenome", help="Input genome.fa (force).")
parser.add_argument("-b", "--refbed", help="Input exon.bed (force).")
parser.add_argument("-c", "--configfile", help="The config file.(default:./config_files/config.ini)", default="./config_files/config.ini")
parser.add_argument("-q", "--qsub_queue_server", help="Submit jobs to queue.(default:all)", choices=["all", "fat", "big"], default="all")
parser.add_argument("-o", "--outputDirectory", help="Output directory with Advanced analysis results.[default:./Advanced_analysis]", default="./Advanced_analysis")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
extenddirectory = os.path.abspath(args.extenddirectory)
genome = os.path.abspath(args.refgenome)
bed = os.path.abspath(args.refbed)
outputDirectory = os.path.abspath(args.outputDirectory)
qsubquene = args.qsub_queue_server
configfile =args.configfile

#Get the software path and Parameter 
config = configparser.ConfigParser()
config.read(configfile,encoding="utf-8")
qsub_pbs = config.get("software", "qsub_pbs")
perl_path = config.get("software", "perl_path")
python3 = config.get("software", "python3")



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
        
        
#Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)   
if not os.path.exists( os.path.dirname(outputDirectory) + "/script/Extend_AS"):
    os.makedirs(os.path.dirname(outputDirectory)+ "/script/Extend_AS")
os.chdir(outputDirectory+ "/../script/Extend_AS")


if os.path.isfile("c.Aspofile_run.sh"):
   os.remove("c.Aspofile_run.sh")
   
if os.path.isfile("d.pblat_run.sh"):
   os.remove("d.pblat_run.sh")

#script for asprofile
scriptName = 'c.Aspofile_run'  + '.sh'
script = open(scriptName, 'a')
script.write("LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:$LD_LIBRARY_PATH"+"&&")
script.write(python3 +" " + sys.path[0] +"/Asprofile.py"+ " ")
script.write("-i "+inputDirectory+ " ")
script.write("-g "+genome+" ")
script.write("-b "+bed+" ")
script.write("-o "+outputDirectory+"\n")
script.close() 

#script for pblat
scriptName = 'd.pblat_run'  + '.sh'
script = open(scriptName, 'a')
script.write("LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:$LD_LIBRARY_PATH"+"&&")
script.write(python3 +" " + sys.path[0] +"/pblat_run.py "+ " ")
script.write("-i "+extenddirectory+" ")
script.write("-o "+extenddirectory+"\n")
script.close() 

#qsub scipt 

command1 = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + "c.Aspofile_run.sh"
subprocess.run(command1, shell=True, check=True)
command2 = perl_path + " " + qsub_pbs + " --queue " + qsubquene + " " + "d.pblat_run.sh"
subprocess.run(command2, shell=True, check=True)
