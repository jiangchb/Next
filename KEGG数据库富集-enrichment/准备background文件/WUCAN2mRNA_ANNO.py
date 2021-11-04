# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 10:53:23 2021

@author: Administrator
"""

path1="Unigene.GO.gene.anno.xls"
Outputfile="mRNA_go.backgroud.xls"

import pandas as pd
import numpy as np

df=pd.read_csv(path1, sep='\t')

#提取无参的表格

df5=df.drop("GO_Anno",axis=1).join(df["GO_Anno"].str.split("|",expand = True).stack().reset_index(level=1,drop = True).rename("Anno_tmp"))
df5=df5.reset_index(drop= True)

#提取GO NUMBER
df6=df5.join(df5["Anno_tmp"].str.extract(r"(GO:\d+)"))
df6.columns = ['gene_id','Num','GO_Anno',"GO_id"]  

#提取term
df6_new=df6.drop('GO_Anno',axis =1).join(df6["GO_Anno"].str.extract(r":\s(.+)\s"))

#\s表示匹配任意空白字符 
#. 表示匹配除 "\n" 之外的任何单个字符。要匹配包括 '\n' 在内的任何字符，请使用象 '[.\n]' 的模式。
#+ 表示匹配前面的子表达式一次或多次，不包括0次
#()表示需要这里的内容

df6_new.columns = ['gene_id','Num','GO_id',"GO_term"]  

#合并
def concat_func(s):
    return pd.Series({
        'GO_id':','.join(s['GO_id'].values),
        'GO_term':"|".join(s['GO_term'].values)
    }
    )

result=df6_new.groupby(['gene_id']).apply(concat_func).reset_index()

result.to_csv(Outputfile,sep='\t',index=False,header=False)

