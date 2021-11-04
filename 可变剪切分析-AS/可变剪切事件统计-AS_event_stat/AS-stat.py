# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 14:01:36 2021

@author: Administrator
"""

import pandas as pd
import numpy as np
import re
import os
#path = "annotation2.xls"
path1 = "LINC00278_AS_result.xls"
path2= "PCDNA3_1_AS_result.xls"
Outputfile="result.xls"
#df=pd.read_csv(path, sep='\t')
df=pd.read_csv(path1, sep='\t')
df2=pd.read_csv(path2, sep="\t")

def create_supportcol(df):
    df["combination"]=df["event_type"].map(str)+":"+df["chrom"].map(str)+","+df["event_start"].map(str)+df["strand"].map(str)+df["event_end"].map(str)+","+df["txpt_len"].map(str)
    return df

df_2=create_supportcol(df)
df2_2=create_supportcol(df2)

def split_table(df):#统计单个基因的可变剪切事件不需要去冗余，因为可能与其他基因共享
    df_result=df[["ref_id","combination","event_type"]]#.drop_duplicates(subset=('combination'))
    df_result2=df_result.drop("ref_id",axis=1).join(df_result["ref_id"].str.split(";",expand = True).stack().reset_index(level=1,drop = True).rename("gene_id"))
    df_result2.dropna(axis=0, how='any', inplace=True)
    return df_result2

df_3=split_table(df_2)
df2_3=split_table(df2_2)

def summary_(df):
    result=df.groupby('gene_id')['event_type'].value_counts().unstack()#得到同一基因不同属性的数量
    result=result.fillna(value=0) #空值填入
    result2=result.apply(lambda x:x.sum(),axis=1).rename("Sum")#求和
    result3=result.merge(result2,on="gene_id")#将求完的和加入进去
    
    return result3

df_result=summary_(df_3)
df2_result=summary_(df2_3)
df_result.to_csv("LINC00278_AS_stat.xls",sep='\t',index=True)
df2_result.to_csv("PCDNA3_1_AS_stat.xls",sep='\t',index=True)

def sum_(df,name):
    result2=df.apply(lambda x:x.sum()).rename(name)
    return result2

df_result2=sum_(df_result,"LINC00278")
df2_result2=sum_(df2_result,"PCDNA3.1")
df_result2.to_csv("LINC00278_AS_sum.xls",sep='\t',index=True)
df2_result2.to_csv("PCDNA3_1_AS_sum.xls",sep='\t',index=True)

def calcuIR(df):
    df["IR_Fraction"]=(df["IR_OFF"]+df["IR_ON"])/df["Sum"]
    return df
#################统计可变剪切的事件，可能在两个不同基因中###################
每个单独event的各个事件类型的数量。
df_result=summary_(df)
df2_result=summary_(df2)

################查找重复的可变剪切#####################################
def find_replicate(df):
    a = df.groupby('combination').count() >1
    #print (df.duplicated(keep = False))
    price = a[a['event_id'] == True].index.to_list() #输出重复的行
    #print (price)
    repeat_df = df[df['combination'].isin(price)] #根据重复的行 输出数据框
    return repeat_df

df_test2=find_replicate(df_test)

#################统计单个基因的可变剪切事件数量，并与表达量合并##############
因为不同的基因可能和同一个可变剪切事件有关系，所以不能随便乱去重复

#拆分表格
df_2=create_supportcol(df)
df2_2=create_supportcol(df2)

df_3=split_table(df_2)
df2_3=split_table(df2_2)

#统计数量
df_result=summary_(df_3)

df2_result=summary_(df2_3)

#计算IR fraction
df_result2= calcuIR(df_result)
df2_result2= calcuIR(df2_result)

#读取表达量
path3= "gene_fpkm_LINC.xls"
path4= "gene_fpkm_PCDNA.xls"
df_fpkm1=pd.read_csv(path3, sep='\t')
df_fpkm2=pd.read_csv(path4, sep='\t')

#合并
df_result3=df_result2.merge(df_fpkm1,on="gene_id")
df2_result3=df2_result2.merge(df_fpkm2,on="gene_id")

#筛选
df_tmp=df_result3[["gene_id","IR_Fraction","LINC00278"]]
df_final=df_tmp[(df_tmp.IR_Fraction!=0)]
df2_tmp=df2_result3[["gene_id","IR_Fraction","PCDNA3_1"]]
df2_final=df2_tmp[(df2_tmp.IR_Fraction!=0)]

df_final.to_csv("LINC00278_AS_stat_forplot.xls",sep='\t',index=False)
df2_final.to_csv("PCDNA3_1_AS_stat_forplot.xls",sep='\t',index=False)

