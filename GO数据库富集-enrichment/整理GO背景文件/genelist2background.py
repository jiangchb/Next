# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 09:40:23 2021

@author: Administrator
"""

import pandas as pd
import numpy as np
import re
import os
#path = "annotation2.xls"
path1 = "annotation.xls"
Outputfile="gene_go_backgroud.xls"
#df=pd.read_csv(path, sep='\t')
df=pd.read_csv(path1, sep='\t')
#清洗数据
list=["--"]
df=df[~(df["Unigene"].isin(list))]

#筛选
df2=df[["GO_ID","GO_classify2","Gene"]]

#拆分基因
df3=df2.drop("Gene",axis=1).join(df["Gene"].str.split("; ",expand = True).stack().reset_index(level=1,drop = True).rename("gene_id"))
df4=df3.reset_index(drop= True)

#定义合并函数
def concat_func(s):
    return pd.Series({
        'GO_id':','.join(s['GO_ID'].values),
        'GO_term':"|".join(s['GO_classify2'].values)
    }
    )

result=df4.groupby(['gene_id']).apply(concat_func).reset_index()

result.to_csv(Outputfile,sep='\t',index=False,header=False)
