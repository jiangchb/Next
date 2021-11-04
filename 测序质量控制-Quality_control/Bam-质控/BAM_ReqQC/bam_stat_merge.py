# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 18:42:01 CST 2017
@author: XiufengYang
to summary genome mapping information from bam_stat.out 
"""
import os.path, time
import glob
import os
import argparse
import string 
import math
import pandas as pd

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Stat filter ratio  scripts.")
parser.add_argument("-i", "--inputDirectory", help="Input directory with bam_stat.out.", default="./Quality_control/RseqQC")
parser.add_argument("-n", "--outputfile",help="Output file with genome mapping stat file", default="./Quality_control/genome_mapping_stat.txt")
args = parser.parse_args()

# Process the command line arguments.

inputDirectory = os.path.abspath(args.inputDirectory)
outputfile = os.path.abspath(args.outputfile)

#define function 
global output

def main(inputdir,outputstat):
   files = glob.glob(inputdir+"/*/*"+"bam_stat.out")
   files = sorted(files)
   for i in range(0, len(files), 1):
      filename=files[i]
      sample_name=(os.path.basename(os.path.dirname(files[i]))).replace("Sample_","")
      file1=filename
      odata = pd.read_table(filename,header=None,sep=':',skiprows=11,nrows=10,names=["Sample",sample_name])
      Total_reads=odata.iloc[0,1]+odata.iloc[1,1]+odata.iloc[2,1]
      Total_mapped_reads=odata.iloc[1,1]+odata.iloc[2,1]
      s1="" 
      sequence=str(Total_mapped_reads)
      line0 = s1.join(sequence)
      odata.loc[odata.index.min() - 1]=["Total reads",Total_reads]
      odata.iloc[0]={"Total mapped reads":1,line0:2}
      odata=odata.sort_index()
      odata=odata.replace("mapq < mapq_cut (non-unique)","Multiple mapped")
      odata=odata.replace("mapq >= mapq_cut (unique)","Uniquely mapped")
      header=odata[["Sample"]]
      odata=odata.drop("Sample",1)
      for j in range(1,len(odata)):
          mapped_ratio=("%.2f" % (int(odata.iloc[j,0])*100/Total_reads) )
          odata.iloc[j,0]=str(odata.iloc[j,0])+"("+str(mapped_ratio)+"%)"
      if i==0: 
           output=pd.concat([odata], axis=1)
      if i>0:
           output=pd.concat([output,odata], axis=1)
   output=pd.concat([header,output], axis=1)
   #output.to_csv(outputstat,sep='\t', encoding='utf-8', index=False)
   output=output.T
   output.to_csv(outputstat,sep='\t', encoding='utf-8', index=True,header=None)

   
   
#action 
if __name__ == '__main__':
   main(inputDirectory,outputfile) 

  
