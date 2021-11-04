# coding:utf-8
#!/usr/bin/env python3
# Author Afan
# Date: 12/29/2019

import argparse
import glob
import os
import subprocess
import shutil
import sys

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Generates htseq-counts script.")
parser.add_argument("-i", "--inputDirectory", help="Input directory with DEG files or gene list file.[default : ./diff_gene] .",
                    default="./diff_gene")
parser.add_argument("-g", "--genome",
                    help="Input directory with gene_kegg.genome.xls, gene_go.genome.xls [default: ./genome]",
                    default="./genome")
parser.add_argument("-t", "--type",
                    help="Input file type,gene or mRNA [default: gene]",
                    default="./genome")
parser.add_argument("-o", "--outputDirectory", help="Output directory with GO and KEGG enrichment files [default:./Enrichment].",
                    default="./Enrichment")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
genome = os.path.abspath(args.genome)
type = args.type
outputDirectory = os.path.abspath(args.outputDirectory)
script_paths = "/public/cluster2/works/lipeng/STUDY/onekey"
# Get the software path
#software = "/media/nbfs/nbCloud/public/task/bin/RNAseq_v1.04/5.Quantification_DEG_Enrichment/enrich"

go_bg = os.path.join(genome + "/" + type + "_go.backgroud.xls")
kegg_bg = os.path.join(genome + "/" + type + "_kegg.backgroud.xls")
category = "/public/cluster2/works/lipeng/pipline/5.Quantification_DEG_Enrichment/enrich/category.xls"

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
else:
    shutil.rmtree(outputDirectory, ignore_errors=True)
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)

os.makedirs("%s/GO_enrichment" % outputDirectory)
os.makedirs("%s/KEGG_enrichment" % outputDirectory)

#os.chdir(outputDirectory)
cmd1 = 'Rscript %s/enrichment.r -i %s -j %s -k %s -c %s -o %s' % (
    script_paths, inputDirectory, go_bg, kegg_bg, category , outputDirectory)
print("step1:"+cmd1)
subprocess.run(cmd1, shell=True, check=True)

GO_result = glob.glob('%s/enrich*go*.xls' % outputDirectory)
GO_list = "".join(GO_result)
if os.path.exists(GO_list):
    GO_enrich=GO_list
else:
    GO_enrich=False
def GO_enrich(GO_list):
    subprocess.call('Rscript %s/top10X3_GO.r -i %s -m Total -o %s/GO_enrichment' % ( script_paths ,GO_list,outputDirectory), shell=True)

KEGG_result = glob.glob('%s/enrich*kegg*.xls' % outputDirectory)
KEGG_list = "".join(KEGG_result)
if os.path.exists(KEGG_list):
    KEGG_enrich=KEGG_list
else:
    KEGG_enrich=False
def KEGG_enrich(KEGG_list):
    subprocess.call('Rscript %s/top20_KEGG.r -i %s -m Total -o %s/KEGG_enrichment' % ( script_paths ,KEGG_list,outputDirectory), shell=True)

if __name__ == "__main__":
    GO_enrich(GO_list)
    KEGG_enrich(KEGG_list)

cmd2 = 'mv %s %s/GO_enrichment '% (GO_list,outputDirectory)
print("step2:"+cmd2)
subprocess.run(cmd2, shell=True, check=True)
cmd3 = 'mv %s %s/KEGG_enrichment '% (KEGG_list,outputDirectory)
print("step3:"+cmd3)
subprocess.run(cmd3, shell=True, check=True)

