#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list = list(
        make_option(c("-l", "--links"), type="character", default=NULL, help="top 300 ppi network results file name", metavar="character"),
	make_option(c("-c", "--config"), type="character", default=NULL, help="top 300 ppi network configuration file name", metavar="character"),
        make_option(c("-o", "--picname"), type="character", default=NULL, help="output picture name", metavar="character")
       ); 
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript Network.r -l top_300_diff-Group1-vs-Group2_gene2gene_network.xls -c top_300_diff-Group1-vs-Group2_configuration.xls -o top_300_diff-Group1-vs-Group2_gene2gene_network.pdf");
opt = parse_args(opt_parser);
if (is.null(opt$links) | is.null(opt$config) | is.null(opt$picname)){
	print_help(opt_parser)
	stop("--links --config --picname must be supplied", call.=FALSE)
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(network))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(sna))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(RColorBrewer))

correlation<-read.table(opt$links,header=T,sep="\t",check.names=F,quote="")
cor<-cbind(as.character(correlation[,1]),as.character(correlation[,3]),"1")
colnames(cor)<-c("source","target","Correlation")
nodes<-read.table(opt$config,header=T,sep="\t",check.names=F,quote="")
nodes<-as.matrix(nodes)

color<-as.character(nodes[,2])
type<-as.character(nodes[,3])
degree<-nodes[,4]
names(color)=nodes[,1]
names(type)=nodes[,1]
names(degree)=nodes[,1] #点大小


edges<-cor
em.net <- edges[,c("source","target")]
em.net <- network::network(em.net, directed = F)
network::set.vertex.attribute(em.net,"color",color[network.vertex.names(em.net)])
network::set.vertex.attribute(em.net,"type",type[network.vertex.names(em.net)])
#network::set.vertex.attribute(em.net,"degree",degree[network.vertex.names(em.net)])
network::set.vertex.attribute(em.net,"degree",sqrt(0.3*(degree(em.net)+2)))
network::set.vertex.attribute(em.net,"size",sqrt(0.3*(sna::degree(em.net)+1)))
color<-factor(get.vertex.attribute(em.net,"color"),levels=c("red", "green"))
#network::set.edge.attribute(em.net,"n",abs()/2)

pdf(opt$picname)
ggnet2(em.net, node.color="color",shape="type", 
       #layout.exp = 0.1,
       #layout.par = list(niter=1000,cell.jitter =1),
       size = "degree",label=TRUE,label.size="size",label.color="#515151",vjust=0,
       edge.alpha = 1,edge.size=0.2,edge.color="#63B8FF",
       #color.legend = "color" 
       alpha=0.9,
       mode= "fruchtermanreingold",
       legend.position = "none",legend.size="size")

dev.off()

