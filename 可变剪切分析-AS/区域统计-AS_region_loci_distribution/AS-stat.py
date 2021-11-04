# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 09:21:27 2021

@author: Administrator
"""
import numpy as np
import pandas as pd
path1 = "intersect_AS_REGION.xls"
Outputfile="result.xls"
#df=pd.read_csv(path, sep='\t')
df=pd.read_csv(path1, sep='\t',header=None)
df=pd.read_excel(path1, sep="\t")

df.columns = ['col1','col2',"col3","col4","col5","col6","col7","col8","col9","col10","col11","col12"] 
df2=df[["col1","col2","col3","col6"]]

df2["combination"]=df["col1"].map(str)+":"+df["col2"].map(str)+"-"+df["col3"].map(str)

df_uniq = df_uniq["combination","col6"]

def concat_func(s):
    return pd.Series({
        #'GO_id':','.join(s['GO_id'].values),
        #'GO_term':"|".join(s['GO_term'].values)
        #'Gene':','.join(s['gene_id'].values)
        "Label":'-'.join(s['col6'].values)
    }
    )
#聚合
result_uniq=df_uniq.groupby(df_uniq['combination']).apply(concat_func)

mapping={"three_prime_UTR-CDS":"CDS-3'UTR","three_prime_UTR":"3'UTR","five_prime_UTR":"5'UTR","five_prime_UTR-CDS":"5'UTR-CDS"}

list_i_want=["CDS-3'UTR","3'UTR","5'UTR","5'UTR-CDS","CDS"]

result_tmp=result_uniq.replace(mapping)

result_final=result_tmp[(result_tmp["Label"].isin(list_i_want))].reset_index()

a=result_final.groupby(result_final["Label"]).size().rename("counts")

a.to_csv(Outputfile,sep='\t',index=True,header=True)

