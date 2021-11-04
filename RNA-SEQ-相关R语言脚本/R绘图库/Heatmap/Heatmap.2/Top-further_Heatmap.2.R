#!/usr/bin/env Rscript
#############################################
#Author: congjia chen
#Creat Time: 2021-5-19
Heatmap.2
#############################################
usage = "\
usage:
    *X基因在所有样本中数值相等将报错 , 需调整输入文件！！*
    /home/hanmin/anaconda3/envs/mro3.5/bin/Rscript further_heatmap.R  -e [fpkm or data matrix file]  -d [phenodata matrix file] -o ./Heatmap
input grouping file format:
    *sample_group.xls*:
    Sample      Group
    Sample_A    group2
    Sample_B    group2
input phenotype file format:
    *anno.xls*:
    gene_id    Celltype
    gene1      Endothelial
    gene2      Fibroblasts
"
cat(usage)


#==========parameter import==========
suppressPackageStartupMessages(library(optparse))
option_list = list(
    make_option(c("-e", "--expression"), type = "character", default = NULL,
        help = "Expression matrix file name(force). "),
    make_option(c("-d", "--group"), type = "character", default = NULL,
        help = "Group information of samplenames.  e.g. sample_group.xls. " ),
    make_option(c("-t", "--title"), type = "character", default = NULL,
        help = "Graphic title information: Group_A-vs-Group_B or \"Group A vs Group B\". "),
    make_option(c("-s", "--showgenes"), type = "logical", default = "F",
        help = "Whether to show gene id, [default: F] . " ),
    make_option(c("-r", "--rowcluster"), type = "logical", default = "T",
        help = "Whether rows (genes) are clustered : T or F , [default: T] . "),
    make_option(c("-l", "--colcluster"), type = "logical", default = "T",
        help = "Whether the columns (samples) are clustered : T or F , [default: T] . "),
    make_option(c("-a", "--angle"), type = "double",  default = 90,
        help = "angle of the column labels, right now one can choose only from few predefined : 0, 45, 90, 270 and 315, [default 90]", metavar = "double"),
    make_option(c("-f", "--fontsize"), type = "double",default = NULL, 
        help = "Base fontsize for the plot, [default %default]", metavar = "double"),
    make_option(c("-p", "--phenotype"), type = "character", default = NULL,
        help = "Phenotype information of genes.  e.g. phenotype.xls. " ),
    make_option(c("-c", "--colors"), type = "character", default = "redwhiteblue",
        help = "colors choise for Heatmap picture : redwhiteblue, redblackgreen ,yellowblackblue , [default: redwhiteblue]. "),
    make_option(c("-g", "--height"), type = "double",  default = NULL,
        help = "Height limit [default %default]", metavar = "double"),
    make_option(c("-k", "--width"), type = "double",  default = NULL,
        help = "Width limit [default %default]", metavar = "double"),
    make_option(c("-o", "--outputdir"), type = "character", default = "./Heatmap",
        help = "output directory for Heatmap results ,[ default: ./Heatmap] . ")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#col=redgreen

#==========import library==========
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(grid))


#==========parameter check==========
if (! is.null(opt$expression)) {
    expression <- read.table(normalizePath(opt$expression), header=T, sep="\t", row.names=1, check.names=F, quote="")
}else {
    print_help(opt_parser)
    stop("expression matrix file is not supplied, heatmap for replicate groups  will be skipped\n", call. = FALSE)
}
if ( !is.null(opt$group) ){
    annotation_col <- read.delim(normalizePath(opt$group), header=T, sep="\t", row.names=1,check.names=F, quote="")
}else {
    annotation_col = NA
}
if ( !is.null(opt$phenotype) ){
    annotation_row <- read.delim(normalizePath(opt$phenotype), header=T, sep="\t", row.names=1,check.names=F, quote="")
}else {
    annotation_row = NA
}
if ( !is.null(opt$title) ){
    title = opt$title
}else {
    title = "Heatmap"
}
angle_list = unlist( strsplit( "0,45,90,270,315", ",", perl = T) )
if ( !opt$angle %in% angle_list ){stop("You can choose only from few predefined : 0, 45, 90, 270 and 315")}
if ( !is.null(opt$angle) ){
    angle = opt$angle
}else {
    angle = 90
}
if ( !is.null(opt$colors) ){
    if ( opt$colors == "redwhiteblue" ){
        palette <- colorRampPalette(c("RoyalBlue2", "White", "Red2"))(n=256)
    }else if ( opt$colors == "redblackgreen" ){
        palette <- colorRampPalette(c("Green", "Black", "Red"))(n=256)
    }else if ( opt$colors == "yellowblackblue" ){
        palette <- colorRampPalette(c("Blue", "Black", "Yellow"))(n=256)
    }
}else {
    palette <- colorRampPalette(c("RoyalBlue2", "White", "Red2"))(n=256)
}
if ( is.null(opt$outputdir) ){
    output_dir = "Heatmap"
}else{
    if ( file.exists(opt$outputdir) ){
        output_dir = opt$outputdir
    }else{
        output_dir = opt$outputdir
        dir.create(output_dir)
    }
}


#==========Limitations==========
#Force gene numbersx 
x <- as.matrix(t(expression))

pdf("output.pdf",width=30,height=10)
p = heatmap.2(x,dendrogram="row",col=bluered,scale="col",trace="none", Rowv = T,
          key = T,
          density.info=c("none"),
          symbreaks = TRUE,
          symkey = T,
          labCol=FALSE
          )  
dev.off()

plot_color = c('orange','green')[treatment]
ColSideColors = plot_color

#
#legend("topright",      
#       legend = unique(x),
#       col = unique(as.numeric(x)), 
#       lty= 1,             
#       lwd = 5,           
#       cex=.7
#)


# heatmap.2(x, dendrogram="col") #只显示列向量的聚类情况
#ggsave(paste0(opt$outpath,"test.png"), height=10, width=10, plot=p)
