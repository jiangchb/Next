#!/usr/bin/env Rscript
#############################################
#Author: Afan 
#Creat Time: 2019-12-17
#############################################
usage = "\
usage:
    Rscript further_volcano.R  -i [DEG file]  -l [genelist file] -o ./Volcano
input file format:
    *genelist.xls*:
    gene_id
    CD133
    OCT4
"
cat(usage)


#==========parameter import==========
suppressPackageStartupMessages(library(optparse))
option_list = list(
    make_option( c("-i", "--input" ), type = "character",
        help = "The input differential gene file(force).  e.g. *-vs-*-all.gene.xls" ),
    make_option( c("-p", "--pval" ), type = "double",default = 0.05,
        help = "pValue ratio threshold, default: 0.05 ."),
    make_option( c("-f", "--foldchange" ), type = "double",default = 2,
        help = "foldchange threshold, default: 2 ."),
    make_option( c("-l", "--genelist" ), type = "character" ,
        help = "Genelist to display the gene symbol.  e.g. genelist.xls"),
    make_option(c("-t", "--title"), type = "character", default = NULL,
        help = "Graphic title and outputfile information: Group_A-vs-Group_B . "),
    make_option( c("-o", "--outputdir" ),type="character", default = NULL,
        help="the output directory of Volcano results, default: A-vs-B." )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#==========import library==========
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(oebio))


#==========parameter check==========
if (! is.null(opt$input) ) {
    DEG <- read.delim(normalizePath(opt$input), header=T, sep="\t", check.names=F, quote="")
}else {
    print_help(opt_parser)
    stop("diff file must be supplied, Volcano plot for groups will be skipped\n", call. = FALSE)
}

if ( !is.null(opt$title) ){
    title = opt$title
}else {
    title = "Volcano"
}

#读取相应的组名

groupname <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", basename(opt$input)))#opt$input

Case_names <- unlist(strsplit(groupname,"-vs-"))[1]
Control_names <- unlist(strsplit(groupname,"-vs-"))[2]

colnames(DEG)[3] <- "Base_control"
colnames(DEG)[4] <- "Base_case"

if ( is.null(opt$outputdir) ){
    output_dir = groupname
    dir.create(output_dir)
}else{
    if ( file.exists(opt$outputdir) ){
        output_dir = opt$outputdir
    }else{
        output_dir = opt$outputdir
        dir.create(output_dir)
    }
}

#==========Limitations==========


if ( !is.null(DEG$pval) ){
    DEG = plyr::rename(DEG, c("pval"="pValue"))
}
if ( !is.null(DEG$foldChange) ){
    DEG = plyr::rename(DEG, c("foldChange"="FoldChange"))
}
columnlist = unlist( strsplit( "log2FoldChange,pValue", ",", perl = T) )
for ( i in columnlist ){
    if ( !i %in% colnames(DEG) ){stop("NO specified column found!")}
}


#==========Volcano function==========
#remove extremum

rownames(DEG)=DEG[,1]
DEG$pValue[which(DEG$pValue < 5E-300 )] = 5E-300

#replace the "-/Inf"
if( "Inf" %in% DEG$log2FoldChange | "-Inf" %in% DEG$log2FoldChange){
    tmp <- DEG[which(DEG$log2FoldChange != "Inf" & DEG$log2FoldChange != "-Inf"),]
    max = max(tmp$log2FoldChange)
    min = min(tmp$log2FoldChange)
    DEG$log2FoldChange <- as.numeric(sub("-Inf", min, DEG$log2FoldChange))
    DEG$log2FoldChange <- as.numeric(sub("Inf", max, DEG$log2FoldChange))
}

#preparation by screening
logFC_cutoff = log(opt$foldchange, 2)
DEG$change = as.factor(ifelse(DEG$pValue < opt$pval & abs(DEG$log2FoldChange) > logFC_cutoff, ifelse(DEG$log2FoldChange > logFC_cutoff ,'Up','Down'),'No change'))
#remove the NaN


#add labels
if ( !is.null(opt$genelist) ) {
    labels <- read.delim(normalizePath(opt$genelist), header=T, sep="\t", check.names=F, quote="")
    DEG$label <-''
    DEG[match(labels[, 1], DEG[, 1]),]$label = as.character(labels[,1])
    #
    if (length(rownames(labels)) <= 15) {
        #sort data
        tmp1 <- DEG[which(DEG$label == ""),]
        tmp2 <- DEG[which(DEG$label != ""),]
        DEG <- rbind(tmp1,tmp2)
        #
        g = ggplot(data=DEG, aes(x=log2(DEG$Base_control), y=log2(DEG$Base_case), color=change)) + 
        geom_point(size=1) + 
        theme_set(theme_set(theme_bw(base_size=15))) +
        theme(legend.title=element_blank()) + 
        labs(title = paste0("Transcriptomes of ",Case_names," versus ", Control_names)) + 
        xlab(paste0("Log2[FPKM] ", Control_names)) + 
        ylab(paste0("Log2[FPKM] ", Case_names)) + 
        theme(legend.position="bottom") +
        theme(plot.title = element_text(hjust = 0.5, vjust = 0.5)) + 
        guides(shape=guide_legend(override.aes=list(size=4))) +
        scale_colour_manual(values = c("Up"=c("red"), "Down"=c("blue"),"No change"=c("grey")),na.translate=FALSE) + 
        geom_label_repel(aes(label = label), size = 4, vjust=-0.5, force = 1,color = "black") + 
        theme(panel.grid = element_blank()) 
    }else {
        #sort data
        tmp1 <- DEG[which(DEG$label == ""),]
        tmp2 <- DEG[which(DEG$label != ""),]
        tmp2$change <- "Labels"
        DEG <- rbind(tmp1,tmp2)
        g = ggplot(data=DEG, aes(x=log2(DEG$Base_control), y=log2(DEG$Base_case), color=change)) + 
        geom_point(size=1) + 
        theme_set(theme_set(theme_bw(base_size=15))) + 
        theme(legend.title=element_blank()) + 
        labs(title = paste0("Transcriptomes of ",Case_names," versus ", Control_names)) + 
        xlab(paste0("Log2[FPKM] ", Control_names)) + 
        ylab(paste0("Log2[FPKM] ", Case_names)) + 
        theme(legend.position="bottom") +
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(shape=guide_legend(override.aes=list(size=4))) +
        scale_colour_manual(values = c("Up"=c("red"), "Down"=c("blue"),"No change"=c("grey"),"Labels"=c("black"))) + 
        theme(panel.grid = element_blank()) 
    }
}else {
        g = ggplot(data=DEG, aes(x=log2(DEG$Base_control), y=log2(DEG$Base_case), color=change)) + 
        geom_point(size=1) + 
        theme_set(theme_set(theme_bw(base_size=15))) + 
        theme(legend.title=element_blank()) + 
        labs(title = paste0("Transcriptomes of ",Case_names," versus ", Control_names)) + 
        xlab(paste0("Log2[FPKM] ", Control_names)) + 
        ylab(paste0("Log2[FPKM] ", Case_names)) + 
        theme(legend.position="bottom") +
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(shape=guide_legend(override.aes=list(size=4))) +
        scale_colour_manual(values = c("Up"=c("red"), "Down"=c("blue"),"No change"=c("grey"))) + 
        theme(panel.grid = element_blank()) 
}

#lines
g = g + geom_abline(intercept = 0,slope = 1, linetype = "dashed", color = c("black"), size = 2)
#
ggsave(file.path(output_dir, paste(groupname,"-comet-pval-",opt$pval,"-FC-",opt$foldchange,".gene.pdf",sep="")),height=8,width=8,plot=g)
ggsave(file.path(output_dir, paste(groupname,"-comet-pval-",opt$pval,"-FC-",opt$foldchange,".gene.png",sep="")),height=8,width=8,plot=g,dpi=1000)

