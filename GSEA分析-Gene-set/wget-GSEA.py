# coding:utf-8
#!/usr/bin/env python3
# Author Congjia.chen
# Date: 12/29/2019

import argparse
import glob
import os
import subprocess
import shutil
import sys
import re
import pandas

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Download from GSEA.")
parser.add_argument("-i", "--inputDirectory", help="Input Download list file.[default : ./geneset.xls] .",
                    default="./geneset.xls")
parser.add_argument("-b", "--ANNOTATION",
                    help="If you want the download content in one file T/F[default: F]",
                    default="F")
parser.add_argument("-o", "--outputDirectory", help="Output directory with GO and KEGG enrichment files [default:./Download].",
                    default="./Download")
args = parser.parse_args()

# Process the command line arguments.
inputDirectory = os.path.abspath(args.inputDirectory)
anno = args.ANNOTATION
outputDirectory = os.path.abspath(args.outputDirectory)

# Create output directories, if they do not exist yet.
if not os.path.exists(outputDirectory):
    os.makedirs(outputDirectory)
os.chdir(outputDirectory)

#read the list
import pandas as pd
import numpy as np
import re
import os
df=pd.read_csv(inputDirectory, sep='\t')

list1=df.iloc[:,0].to_list()

##########Allin or seperate #####################
if anno == "F":
	for i in list1:
		web="http://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=%s&fileType=gmt" % (i) 
		command="wget -c %s" % (i, web)
		print( command )
#		wget.download(DATA_URL, out=out_fname)
		subprocess.run(command, shell=True, check=False)
else:
	for i in list1:
		web="http://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=%s&fileType=gmt" % (i) 
		command="wget -O %s.txt -c %s" % ("all_gene_sets", web)
		print( command )
		subprocess.run(command, shell=True, check=False)
	print ("The download content will be saved in one file")
print ("##################################################Script generate successfully##############################################")

