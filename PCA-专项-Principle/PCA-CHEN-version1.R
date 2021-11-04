if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#使用counts定量文件绘制PCA情况

library("DESeq2")
counts="gene_counts.xls"
group="sample_group.txt"

counts <- round(read.delim(counts, row.names = 1, check.names = F,quote=""))
phenodata <- read.table(group, row.names = 1, header = T, sep = "\t", check.names = F)

###################################################################################################
#y
counts_tmp<-subset(counts, select = c(rownames(phenodata))) 
counts<-counts_tmp
###################################################################################################
samples <- colnames(counts)
index_phenodata <- phenodata
for (i in 1 : length(samples))
{
  index_phenodata[i,] <- phenodata[samples[i],]
}
print(index_phenodata)
#fixed the group name showing error by define the  level of index_phenodata  
rownames(index_phenodata)<-c(as.character(samples))
index_phenodata$Group<-factor(index_phenodata$Group,level=unique(index_phenodata$Group))
groups <- as.character(index_phenodata$Group) 

#点的颜色完全是根据pcaData的顺序加的
##modified 20170911 ,show the group name  as same as  plot-3d  by reindex coldata and counts 
tmp_coldata <- data.frame(row.names = samples, Group = groups)
colData<-tmp_coldata[order(tmp_coldata$Group), , drop = FALSE]
counts<-counts[ ,c(rownames(colData))]
groups <-c(as.character(colData$Group))
samples <- colnames(counts)

#OE有参常规脚本
#创建一个DEseq的对象，指定矩阵，实验设计以及分类

dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Group) #画PCA图需要Deseq2 

##样本数量<30,使用rlog方法进行归一化

if (length(samples) < 3)
{cat("Samples number is less than 3,stop!!", call. = FALSE)}

#rld <- rlog(dds, blind = T)
if (length(samples) < 30){rld <- rlog(dds, blind = T)#Deseq2包里面的标准化的函数 负二项式分布
}else{##样本数量>30,使用vsd归一化，速度快，质量好
rld <- vst(dds, blind = FALSE)}

##PCA 2D
#至此得到标准化后的表达量表格
pcaData <- plotPCA(rld, intgroup = c("Group"), returnData = TRUE) #plotPCA
#intagroup实际上这个分组最后影响的是PCA图中的颜色，但是并不影响PCA图中各个样本的位置

#上述函数其实可以直接生成PCA的图，但这里是用来生成绘图用的数据
min <- round(min(pcaData$PC1, pcaData$PC2))
max <- ceiling(max(pcaData$PC1, pcaData$PC2))


print(attr(t(pcaData), "percentVar"))
percentVar <-  round(100 *attr(pcaData, "percentVar"),2)


#################################################

min <- round(min(pcaData$PC1, pcaData$PC2))
max <- ceiling(max(pcaData$PC1, pcaData$PC2))

#################################################

###PCA 3D
groups <- as.character(colnames(counts))
#vsd <- varianceStabilizingTransformation(dds, blind = T)
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = T)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
max_x<-round(percentVar[1]*100,2)
max_y<-round(percentVar[2]*100,2)
max_z<-round(percentVar[3]*100,2)


##################################################
ggplot(pcaData, aes(PC1, PC2, color = Group)) +
  geom_point() +
  xlab(paste0("PC1: ", round(percentVar[1],2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2],2), "% variance")) +
  guides(col = guide_legend(nrow = 15)) +
  scale_color_manual(values = oe_col_qua(1 : length(seq_along(levels(pcaData$Group)))),
                     labels = c(as.character(unique(pcaData$Group))),
                     name = "Group") +
  scale_fill_manual(values = oe_col_qua(1 : length(seq_along(levels(pcaData$Group)))),
                    labels = c(as.character(unique(pcaData$Group))),
                    name = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
  theme(text = element_text(size = 15, family = "ArialMT")) +
  scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
  scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
  geom_point(size = 4) +
  geom_vline(xintercept = 0, linetype = 4, color = "grey") +
  geom_hline(yintercept = 0, linetype = 4, color = "grey") +
  stat_ellipse(level = 0.95, show.legend = F)+
  coord_equal()          
#	ggsave("PCA_deseq2_replicate.pdf",height=4*0.05*percentVar[2],width=4*0.05*percentVar[1])
ggsave("PCA_1.pdf", height = 8, width = 8)
ggsave("PCA_1.png",type="cairo-png",height=8,width=8)