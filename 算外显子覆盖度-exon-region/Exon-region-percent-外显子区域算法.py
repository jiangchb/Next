# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 14:24:37 2021

@author: nedch
"""

import pandas as pd
import numpy as np

class Solution2(object):
    def read_information(self, path1,path2):
        df=pd.read_csv(path1, sep='\t',header=None)
        df2=pd.read_csv(path2, sep='\t',header=None)
        df_read=df.loc[:,[2,9,10]]
        df_read.columns=["Chrom_id","sp","ep"]
        df2.columns=["Chrom_id","Chrom_szie"]
        return df_read,df2
    def split_col(self, df_read):
        '''
        

        Parameters
        ----------
        df_read : TYPE dataframe
            DESCRIPTION.

        Returns
        -------
        df_result : TYPE dataframe with split col
            DESCRIPTION.

        '''
        dfa=df_read.drop("ep",axis=1).join(df_read["sp"].str.split(",",expand = True).stack().reset_index(level=1,drop = True).rename("start point"))
        dfb=df_read.drop("Chrom_id",axis=1).join(df_read["ep"].str.split(",",expand = True).stack().reset_index(level=1,drop = True).rename("end point"))
        df_ab= pd.concat([dfa, dfb], axis=1).reset_index()
        df_ab.replace('', np.nan, inplace=True)
        df_result=df_ab.dropna(axis=0)[["start point","end point","Chrom_id"]]
        #df_result["Comb_set"]=df_result["start point"]+","+df_result["end point"]
        return df_result
    def test(self,x):
        '''

        Parameters
        ----------
        x : TYPE    Groupby object
            DESCRIPTION.     for the combination of Groupby and apply

        Returns
        -------
        TYPE     Series
            DESCRIPTION.    

        '''
        #x["starp"]=x["Comb_set"].map(lambda x:x.split(",")[0])
        #x["endp"]=x["Comb_set"].map(lambda x:x.split(",")[1])
        #x=x[["starp","endp"]]
        x=x[["start point","end point"]]
        intervals = x.astype(int).values.tolist()
        intervals = sorted(intervals, key=lambda x:x[0]) #用列表中的第一个元素来进行排序
        res = []
        exon_length=[]
        start, end = intervals[0][0], intervals[0][1]
        for i in range(len(intervals)):
            s, e = intervals[i][0], intervals[i][1]
            if s <= end: # overlap
                end = max(end, e) #当发现重叠时，最大值为两个end的最大。
                    #if 的过程是不断累积的过程，当没发现有overlap之后，就归入else，将最终答案写进结果。然后start,end又开始一个新的轮回。
            else:
                length=end-start+1
                exon_length.append(length)
                res.append([start,end,length])
                start, end = s, e 
        res.append([start, end])
        length=end-start+1
        exon_length.append(length)#这个是为了避免最后一个轮回的时候，会有缺失。
        length_summary=sum(exon_length)
            #这个算法值得学习
        return pd.Series({
        "exon_length": length_summary
        }
        )
    def calculperexon(self,exon,chrom):
        df_tmp2=exon.merge(chrom,on="Chrom_id")
        df_tmp2["exon_percent"]=df_tmp2["exon_length"]/df_tmp2["Chrom_szie"]
        return df_tmp2


merge_df=Solution2()
path1="ncbiRefSeqCurated.txt"
path2="hg38.chrom.sizes.txt"
df_test,df_chrinfo=merge_df.read_information(path1, path2)
df_tmp1=merge_df.split_col(df_test)
result=df_tmp1.groupby(df_tmp1['Chrom_id']).apply(merge_df.test).reset_index()
exonper=merge_df.calculperexon(result,df_chrinfo)
