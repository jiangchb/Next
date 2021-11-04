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
import subprocess
import numpy as np
import sys 

# Read the command line arguments.
parser = argparse.ArgumentParser(description="Stat filter ratio  scripts.")
parser.add_argument("-f", "--align_stat", help="Input genome_mapping_stat.xls")
parser.add_argument("-t", "--thresh", help="Genome mapping ratio thresh.",default="0.80")
parser.add_argument("-s","--sample_num", help="Sample_number from  config")
parser.add_argument("-e", "--email", help="email address for results.")
parser.add_argument("-o", "--outputfile", help="output check result ")
args = parser.parse_args()

# Process the command line arguments.

align_stat_file = os.path.abspath(args.align_stat)
thresh = args.thresh
outputfile = os.path.abspath(args.outputfile)
email = args.email
sample_num = args.sample_num 

if os.path.exists(outputfile):
    os.remove(outputfile)
def main(align_stat_file, thresh,sample_num, outputfile, email):
    data = pd.read_table(align_stat_file, sep='\t', header=0, index_col=0)
    output = open(outputfile, 'a',newline='',encoding='utf-8')
    output.write(os.path.dirname(align_stat_file)+"：基因组比对情况反馈:"+"\n")
    output.write(u"config文件中提供样本数目："+str(sample_num)+"\n")
    output.write(u"基因组统计文件中样本数目："+str(len(data.columns))+"\n")
    if sample_num != str(int(len(data.columns))):
        output.write(u"！！！！！！报警：config文件中提供样本数目和基因组统计文件中样本数目不吻合，请仔细核查！！！！！！"+"\n")
    for j in range(0,len(data.columns)):
        sample_name = data.columns[j]
        mapped_ratio = float(data.iloc[1,j].split("(")[1].replace("%)",""))
        read_12_ratio = round(abs(float(data.iloc[4,j].split("(")[1].replace("%)",""))-float(data.iloc[5,0].split("(")[1].replace("%)",""))),2)
        thresh_ratio = float(thresh)*100
        if mapped_ratio > thresh_ratio :
            output.write(u'样本%s基因组比对率%s 超过设定阈值%s ,正常; '%(sample_name ,str(mapped_ratio)+"%",str(thresh_ratio)+"%"))
            if read_12_ratio > 5 : 
                output.write(u'！！！报警:R1,R2比对率差别较大，为%s, 超过%s ,请核查数据质量！！！\n'%(str(read_12_ratio)+"%","5%"))
            else:
                output.write(u'R1,R2比对率差别不大，为%s, 未超过%s ,正常.\n'%(str(read_12_ratio)+"%","5%"))
        else:
            output.write(u'！！！！！！报警：样本 %s 基因组比对率 %s,低于设定阈值 %s ,不正常,请仔细核查！！！！！！'%(sample_name ,str(mapped_ratio)+"%",str(thresh_ratio)+"%"))
            if read_12_ratio > 5 : 
                output.write(u'！！！报警:R1,R2比对率差别较大，为%s, 超过%s ,请核查数据质量！！！\n'%(str(read_12_ratio)+"%","5%"))
            else:
                output.write(u'R1,R2比对率差别不大，为%s, 未超过%s ,正常.\n'%(str(read_12_ratio)+"%","5%"))
    output.close()
    cmd1 = 'cat  %s %s | mailx  -s %s  -a %s %s' %(outputfile,align_stat_file,"基因组比对率反馈："+os.path.dirname(align_stat_file), align_stat_file, email)
    subprocess.run(cmd1.encode(encoding="utf-8"), shell=True, check=True)
#action
if __name__ == '__main__':
    main(align_stat_file, thresh,sample_num, outputfile, email)
