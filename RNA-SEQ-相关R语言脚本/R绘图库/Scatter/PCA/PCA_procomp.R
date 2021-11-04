library(ggplot2)
library(ggrepel)
suppressPackageStartupMessages(library(ggrepel))
#suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(RColorBrewer))
#suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(grDevices))
#suppressPackageStartupMessages(library(genefilter))
#suppressPackageStartupMessages(library(pca3d))
#suppressPackageStartupMessages(library(maptools))
#suppressPackageStartupMessages(library(oebio))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(grid))

path1 <- "protein_quan24.xls"
path2 <- "sample_group24.txt"

#=====================读取表达量矩阵==========================
counts <- read.delim(path1, sep="\t", header=T, quote="",row.names = 1)
#标准化
#data2=scale(d1, center=F, scale=T)
#=====================读取组名==========================
phenodata <- read.table(path2, row.names = 1, header = T, sep = "\t", check.names = F)
###################################################################################################
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
##modified 20170911 ,fixed the group name showing error by define the  level of index_phenodata  
rownames(index_phenodata)<-c(as.character(samples))
index_phenodata$Group<-factor(index_phenodata$Group,level=unique(index_phenodata$Group))
groups <- as.character(index_phenodata$Group) 

##modified 20170911 ,show the group name  as same as  plot-3d  by reindex coldata and counts 
tmp_coldata <- data.frame(row.names = samples, Group = groups)
colData<-tmp_coldata[order(tmp_coldata$Group), , drop = FALSE]
counts<-counts[ ,c(rownames(colData))]
groups <-c(as.character(colData$Group))
samples <- colnames(counts)

#=====================创建PCA===========================
#PCA加上标准化
pca <- prcomp(t(counts), center=T, scale.=T)
#创建PCA的数据表
pcaData <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],  group=groups, colData, name=samples)

# pca data
pca_data = data.frame(rownames(pca$x), pca$x, check.names = F)
colnames(pca_data)[1] = "IDs"
write.table(pca_data, file = "PCA_value.xls", row.names = FALSE, quote = FALSE, sep='\t')

percentVar <- pca$sdev^2/sum(pca$sdev^2)

min <- round(min(pcaData$PC1, pcaData$PC2))
max <- ceiling(max(pcaData$PC1, pcaData$PC2))

x <- round(percentVar[1]*100,2)
y <- round(percentVar[2]*100,2)
z <- round(percentVar[3]*100,2)

# PCA 2D-1
pca1 = ggplot(pcaData, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", x, "%)")) +
  ylab(paste0("PC2 (", y, "%)")) +
  guides(col = guide_legend(nrow = 15)) +
  scale_color_manual(values = c("#006400","#8B0000","#FF8C00"),
                     labels = c(as.character(unique(factor(pcaData$Group)))), name = "Group") +
  scale_fill_manual(values = c("#006400","#8B0000","#FF8C00"),
                    labels = c(as.character(unique(factor(pcaData$Group)))), name = "Group") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
  theme(text = element_text(size = 15, family = "ArialMT")) +
  scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
  scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
  geom_vline(xintercept = 0, linetype = 4, color = "grey") +
  geom_hline(yintercept = 0, linetype = 4, color = "grey") +
#  stat_ellipse(aes(fill = Group, group = Group),geom = "polygon", level = 0.95, alpha = 0.5) +#fill = Group,
  coord_equal()
ggsave("PCA_1.pdf", height = 8, width = 8,pca1)
ggsave("PCA_1.png",type="cairo-png",height=8,width=8,pca1)


pca2 = pca1 + geom_text_repel(aes(x = PC1, y = PC2, label = rownames(pcaData)), size = 3) 

ggsave("PCA_2.pdf", height = 8, width = 8,pca2)
ggsave("PCA_2.png",type="cairo-png",height=8,width=8,pca2)

