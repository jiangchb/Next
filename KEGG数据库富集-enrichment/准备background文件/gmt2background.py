# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 09:48:55 2021

@author: Administrator
"""

import pandas as pd
import numpy as np
import re
import os

path1 = "Homo_sapiens.path"
Outputfile="result.xls"
#df=pd.read_csv(path, sep='\t')
df=pd.read_csv(path1, sep='\t')

#清洗数据
list=["--"]
df=df[~(df["Gene(Unigene)"].isin(list))]

#筛选
df2=df[["Proteins","#Pathway","Pathway ID"]]

#拆分基因
df3=df2.drop("Proteins",axis=1).join(df["Proteins"].str.split(";",expand = True).stack().reset_index(level=1,drop = True).rename("Proteins"))
df4=df3.reset_index(drop= True)

#定义合并函数
def concat_func(s):
    return pd.Series({
        'Pathway_ID':','.join(s['Pathway ID'].values),
        'Pathway_Term':"|".join(s['#Pathway'].values)
    }
    )

result=df4.groupby(['Proteins']).apply(concat_func).reset_index()

result.to_csv(Outputfile,sep='\t',index=False,header=False)
