#!/usr/bin/env Rscript
data <- read.table("./VennData.xls", sep="\t",header=T, quote="", check.names=FALSE, row.names=1)
library("venn")
library("ggplot2")
pdf("./VennGraph.pdf", width=10, height=10, onefile=F)
venn(data,zcolor="red,yellow,blue,green,orange,purple")
dev.off()

png("./VennGraph.png", width=3000, height=3000, res=300)
venn(data,zcolor="red,yellow,blue,green,orange,purple")
dev.off()