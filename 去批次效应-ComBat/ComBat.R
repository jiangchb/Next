#Batch-effect

# http://www.bioconductor.org/packages/release/bioc/html/sva.html
# by wsb 20200107
# ComBat-Seq 消除批次效应
# 安装sva包
#    if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 安装 ComBat_seq
#    devtools::install_github("zhangyuqing/sva-devel")

library("sva")

'''
遇到出错
I encountered the exact same error when I try to remove known batches using combat function, I found 2 problems matters:
  
1.the type of 'dat', it should be matrix instead of data.frame;
2.variance of variables in 'dat' should not equal zero.
Once these 2 conditions satisfied, you can run combat successfully.
'''

path1="gene_counts_new.xls"

cdata <- read.delim("gene_counts_new.xls", header = T, sep = "\t", row.names = 1,comment.char = '')

cdata[is.na(cdata)] <- 0
cdata <- as.matrix(cdata)
"""
#Force variance not to be zero！
library("dplyr")
cdata["sd"] <- apply(cdata, 1, sd)
cdata2 <- filter(cdata , sd != 0)
cdata2 <- as.matrix(cdata2)
#
"""
#转换为矩阵
csif <- read.table("batch_infor.txt", header = T, sep = "\t", row.names = 1)

modcombat = model.matrix(~1, data = csif)

batch = csif$batch

combat_edata = ComBat_seq(cdata, batch=batch,group = NULL)
#,mod=modcombat
#去除批次效应
write.table(combat_edata, "ComBat_data(all)(counts).xls", sep = "\t", quote = F)