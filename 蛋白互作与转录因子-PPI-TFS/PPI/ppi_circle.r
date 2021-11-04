#!/usr/bin/env Rscript
# by The Coder, 20160527
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="top300 ppi", metavar="character"),
	make_option(c("-d", "--diff"), type="character", default=NULL, help="diff file", metavar="character"),
	make_option(c("-n", "--number"), type="character", default="30", help="show number", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript ppi_circle.r -i top_300_diff-A-vs-B_gene2gene_network.txt -o outdir/");
opt = parse_args(opt_parser);
if(is.null(opt$input) | is.null(opt$diff)){
	print_help(opt_parser)
	stop("--input --diff must be supplied", call.=FALSE)
}
if(!file.exists(opt$input)){
	print("input is not exists")
	q()
}

library(tidyverse)
library(igraph)
library(dplyr)

df_ppi <- read.table(opt$input,sep="\t",header = T)
df_fc <- read.delim(opt$diff,sep="\t",header = T)

#
if(dim(df_ppi)[1] <1){
	stop("ppi fill is empty")
}


#get links
top20 = head(df_ppi,n=opt$number)
diffname <- gsub("top_300_diff-","",gsub("_gene2gene_network.txt", "", basename(opt$input)))

#cal gene nodes and weight,add FC

list1 = as.character(top20$gene1)
list2 = as.character(top20$gene2)
list = as.character(cbind(list1,list2))
df_num = as.data.frame(table(list))
names(df_num) = c("gene_id","num")
gene2Num = df_num[order(df_num[,2],decreasing=T),]
gene2FC = df_fc[c('gene_id','FoldChange')]
mergr_nodes = merge(gene2Num,gene2FC,by.x="gene_id",by.y="gene_id",all.x=T)
#get color 

df_Up = mergr_nodes[which(mergr_nodes[,3]>= 2 ),]
df_Down = mergr_nodes[which(mergr_nodes[,3]<=0.5 ),]
df_Up = df_Up[order(df_Up[,2],decreasing=F),]
df_Down = df_Down[order(df_Down[,2],decreasing=F),]
palette_Up <- colorRampPalette(c("grey","red"))(n=length(rownames(df_Up)))
palette_Down <- colorRampPalette(c("steelblue","grey"))(n=length(rownames(df_Down)))
df_Up$color = palette_Up
df_Down$color = palette_Down
my_nodes= rbind(df_Up,df_Down)

my_link = top20[c(1,3,5)]
#rm NA if less than top num
my_link <-my_link %>% drop_na()
colnames(my_link)[3] = "weight"

color_num = length(rownames(df_num)):1
#get network


net<-graph.data.frame(my_link,my_nodes,directed = F)
#nodes size &color
deg <- igraph::degree(net, mode="all")
V(net)$size <- log2(deg*50 + 50) 

list_sort = as.data.frame(V(net)$name)
names(list_sort) = "gene_id"
re_col = left_join(list_sort,my_nodes[,c(1,4)],by="gene_id")
V(net)$color = re_col$color

#E(net)$label<-E(net)$weight/1000 #edges label

#move node_lable location
#radian.rescale <- function(x, start=0, direction=1) {
#  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
#  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
#}
#n = length(my_nodes)
#lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)
pdf(paste0(opt$outpath, "/top_",opt$number,"_diff-", diffname,"_network.pdf"), height=10, width=10)

plot(net,layout=layout.circle, 
#     vertex.color=palette,
     vertex.label.dist = 1,
#     vertex.label.degree=lab.locs,
     edge.arrow.size=0, #设置箭头大小
     vertex.frame.color="transparent",  #节点边框透明
     vertex.label.color="black",        #节点标签黑色
     vertex.label.cex=0.65,             #节点标签大小
     edge.curved=0.3,                   #边是否弯曲，取值0-1，0为不弯曲
     edge.color="grey"  ,main = paste0("top ",opt$number,"_diff-",diffname,"_network")             #边的颜色
)
dev.off()