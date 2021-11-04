#!/usr/bin/env Rscript
#############################################
#Author: Afan 
#Creat Time: 2019-11-15
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


#==========import library==========
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(oebio))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))


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
#Force gene numbers
genes <- length(rownames(expression))
if (( opt$showgenes == F ) & ( !is.null(opt$fontsize))) {
    show_rownames = F
    fontsize = opt$fontsize
}else if (( opt$showgenes == F ) & ( is.null(opt$fontsize))){
    show_rownames = F
    fontsize = 12
}else if (( opt$showgenes == T ) & ( !is.null(opt$fontsize))){
    if ( genes > 150 ){
        show_rownames = F
        fontsize = opt$fontsize
    }else{
        show_rownames = T
        fontsize = opt$fontsize
    }
}else if (( opt$showgenes == T ) & ( is.null(opt$fontsize))){
    if ( genes <= 30 ){
        show_rownames = T
        fontsize = 14
    }else if ( (genes > 30) & (genes <= 100) ){
        show_rownames = T
        fontsize = 10
    }else if ( (genes > 100) & (genes <= 150) ){
        show_rownames = T
        fontsize = 6
    }else{
        show_rownames = F
        fontsize = 12
    }
}else{
    print ("Done")
}

#Force variance not to be zero！
ind <- apply(expression, 1, mean) > 0
expression <- expression[ind, ]
sd <- apply(expression, 1, sd)
if (length(sd[sd == 0])>0) {
    stop("sd = 0, program exit!", call. = FALSE)
}
#Limit gene length to 30！
limit30 <- function(k) {
    str_length(k)
}
character = as.data.frame(rownames(expression))
clength <- apply(character , 1, limit30)
if ( any(clength > 30)) {
    show_rownames = F
}


#==========pheatmap function==========
express_file_abpath <- normalizePath(opt$expression)
path <- as.character(unlist(strsplit(express_file_abpath, split = "/")))
picname <- sub(".txt|.xls$", "", path[length(path)] )

samples <- colnames(expression)
border = F

if (length(samples) < 2) {
        cat("less than 2 diff gene/transcript,heatmap cluster is skipped!!")
}else if (length(samples) == 2) {
    if (opt$rowcluster == T ) {
        p = pheatmap(log2(expression+0.0001),
        main = title,
        treeheight_row=30, lwd=1,
        color = palette , show_rownames = show_rownames, fontsize = fontsize, border = border,
        cluster_rows = T , cluster_cols = F , angle_col = angle )
        
        data = pheatmap(log2(expression+0.0001))
        order_row = data$tree_row$order
        datat = data.frame(expression[order_row,])
        datat = data.frame(rownames(datat),datat,check.names =F)
        colnames(datat)[1] = "gene_id"
        write.table(datat,file=paste0(output_dir ,"/", picname,".genes_reorder.cluster_result.xls"),row.names=FALSE,quote = FALSE,sep='\t')
    }else {
        p = pheatmap(log2(expression+0.0001),
        main = title,
        lwd=1, color = palette , show_rownames = show_rownames, fontsize = fontsize , border = border,
        cluster_rows = F , cluster_cols = F , angle_col = angle )
    }
}else {
    p = pheatmap(log2(expression+0.0001),
    main = title, annotation_col = annotation_col, annotation_row = annotation_row, scale="row",
    treeheight_row=30, treeheight_col=30, lwd=1,
    color = palette , show_rownames = show_rownames, fontsize=fontsize, border=border,
    cluster_rows = opt$rowcluster , cluster_cols = opt$colcluster , angle_col = angle )

}


if ((!is.null(opt$height)) & (!is.null(opt$width))){
    height = opt$height
    width = opt$width
}else if ((is.null(opt$height)) & (is.null(opt$width))){
    if ( opt$showgenes == F ){
        height = 7
        width = 7
    }else{
        height = 12+0.5*length(samples)
        width = 8+0.5*length(samples)
    }
}else if ((is.null(opt$height)) || (is.null(opt$width))){
    stop("height option must with width option, heatmap will be skipped\n", call. = FALSE)
}

    save_pheatmap <- function(x, filename) {
    stopifnot(!missing(x))
    stopifnot(!missing(filename))
    #pdf
    pdf(paste0(output_dir, "/", filename, ".heatmap.pdf"), height = height, width = width, onefile = F)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
    #png
    png(paste0(output_dir, "/", filename, ".heatmap.png"), height = height*240 , width = width*260 , res = 250)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
    }

    save_pheatmap(p, picname)


#==========cluster dataframe==========
if (length(samples) > 2) {

    if ( (opt$rowcluster == T) & (opt$colcluster == F ) ){
    
        data = pheatmap(log2(expression+0.0001),scale="row", cluster_rows = T , cluster_cols = F )
        order_row = data$tree_row$order
        datat = data.frame(expression[order_row,])
        datat = data.frame(rownames(datat),datat,check.names =F)
        colnames(datat)[1] = "gene_id"
        write.table(datat,file=paste0(output_dir ,"/", picname,".genes_reorder.cluster_result.xls"),row.names=FALSE,quote = FALSE,sep='\t')
    }else if( (opt$rowcluster == F) & (opt$colcluster == T ) ) {
    
        data = pheatmap(log2(expression+0.0001),scale="row" , cluster_rows = F , cluster_cols = T )
        order_col = data$tree_col$order
        datat = data.frame(expression[,order_col])
        datat = data.frame(rownames(datat),datat,check.names =F)
        colnames(datat)[1] = "gene_id"
        write.table(datat,file=paste0(output_dir ,"/", picname,".genes_reorder.cluster_result.xls"),row.names=FALSE,quote = FALSE,sep='\t')
    }else if( (opt$rowcluster == T) & (opt$colcluster == T ) ){
        data = pheatmap(log2(expression+0.0001),scale="row",cluster_rows = T , cluster_cols = T)
        order_row = data$tree_row$order
        order_col = data$tree_col$order
        datat = data.frame(expression[order_row,order_col])
        datat = data.frame(rownames(datat),datat,check.names =F)
        colnames(datat)[1] = "gene_id"
        write.table(datat,file=paste0(output_dir ,"/", picname,".genes_reorder.cluster_result.xls"),row.names=FALSE,quote = FALSE,sep='\t')
    }else {
        print ("Done")
    }
}

if ( file.exists("Rplots.pdf") ){ file.remove("Rplots.pdf") }
