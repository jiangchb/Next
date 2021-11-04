# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 11:14:52 2021

@author: Administrator
"""

path1="Unigene.KEGG.gene.anno.xls"
Outputfile1="mRNA_kegg.backgroud.xls"
Outputfile2="mRNA_anno-kegg.backgroud.xls"

import pandas as pd
import numpy as np

df=pd.read_csv(path1, sep='\t')

#生成mRNA_anno-kegg.backgroud.xls

df2=df[["#Gene_id","KO","Pathway"]]

#替换 |        
df3 = df2.drop('Pathway',axis =1).join(df2["Pathway"].str.replace("|",","))

#去除包含--的列
方法1
df3_final= df3[df3['Pathway'].str.contains('--') == False]

方法2
df3_final=df3[~df3['Pathway'].str.contains('--')]


df3_final.to_csv(Outputfile2,sep='\t',index=False,header=False)


#生成mRNA_kegg.backgroud.xls

df4=df[["#Gene_id","Pathway","Pathway_definition"]]

df5 = df4.drop('Pathway',axis =1).join(df4["Pathway"].str.replace("|",","))

df5=df5[["#Gene_id","Pathway","Pathway_definition"]]
#去除包含--的列
df5_final=df5[~df5['Pathway'].str.contains('--')]

df5_final.to_csv(Outputfile1,sep='\t',index=False,header=False)
  
        