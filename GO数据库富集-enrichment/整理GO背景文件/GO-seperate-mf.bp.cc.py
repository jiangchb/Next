# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 08:53:40 2021

@author: Administrator
"""

import pandas as pd
import numpy as np
import re
import os

path1 = "Homo_sapiens_GO_detail.xls"
Outputfile="go_backgroud.xls"
#df=pd.read_csv(path, sep='\t')
df=pd.read_csv(path1, sep='\t')

df_BP=df[["Protein_ID","GO Biological_Process"]]
df_MF=df[["Protein_ID","GO Molecular_Function"]]
df_CC=df[["Protein_ID","GO Cellular_Component"]]

def splitGOdf(df2):
    df2=renameGO(df2)
    df3=df2.drop("GO_TERM",axis=1).join(df2["GO_TERM"].str.split(";",expand = True).stack().reset_index(level=1,drop = True).rename("GO"))
    df3=df2.drop("GO_TERM",axis=1).join(df2["GO_TERM"].str.split(";",expand = True).stack().reset_index(level=1,drop = True).rename("GO"))
    df4=df3.drop("GO",axis=1).join(df3["GO"].str.split("//",expand = True))
    df4.columns = ['Protein_ID','GO_ID',"GO_classify2"]
    df4=df4.dropna(axis=0,how="any",inplace=False)
    result=df4.groupby(['Protein_ID']).apply(concat_func).reset_index()
#    list=["nan"]
#    result=result[~(result["GO_id"].isin(list))]
    return result

def renameGO(df):
    return df.rename(columns={ df.columns[1]: "GO_TERM" }, inplace=False)

def concat_func(s):
    return pd.Series({
        'GO_id':','.join(s['GO_ID'].map(str).values),
        'GO_term':"|".join(s['GO_classify2'].map(str).values)
    }
    )

result_BP=splitGOdf(df_BP)
result_BP2=splitGOdf(df_BP)
result_MF=splitGOdf(df_MF)
result_CC=splitGOdf(df_CC)
    

    
result_BP.to_csv("GO.BP.background.xls",sep='\t',index=False,header=False)
result_MF.to_csv("GO.MF.background.xls",sep='\t',index=False,header=False)
result_CC.to_csv("GO.CC.background.xls",sep='\t',index=False,header=False)


def splitGOdf2(df2):
    df2=renameGO(df2)
    df3=df2.drop("GO_TERM",axis=1).join(df2["GO_TERM"].str.split(";",expand = True).stack().reset_index(level=1,drop = True).rename("GO"))
    df3=df2.drop("GO_TERM",axis=1).join(df2["GO_TERM"].str.split(";",expand = True).stack().reset_index(level=1,drop = True).rename("GO"))
    df4=df3.drop("GO",axis=1).join(df3["GO"].str.split("//",expand = True))
    df4.columns = ['Protein_ID','GO_ID',"GO_classify2"]  
    return df4

result_BP=splitGOdf2(df_BP)
result_MF=splitGOdf2(df_MF)
result_CC=splitGOdf2(df_CC)

list=[result_BP,result_MF,result_CC]

result_tmp=pd.concat(list, axis=0, join='outer',ignore_index = True)

#去重
result_tmp2 = result_tmp.drop_duplicates(subset=['GO_ID',"Protein_ID","GO_classify2"])

#list=["nan"]
#result_tmp=result_tmp[~(result_tmp["GO_ID"].isin(list))]  
result_tmp2=result_tmp2.dropna(axis=0,how="any",inplace=False)
result=result_tmp2.groupby(['Protein_ID']).apply(concat_func).reset_index()

result.to_csv("GO.ALL.background.xls",sep='\t',index=False,header=False)





