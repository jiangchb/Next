#!/usr/bin/env Rscript
#############################################
#Author: Afan 
#Creat Time: 2019-12-17
#############################################
"2021/8/25 Updated"
#############################################
#Author: Chen congjia 
#Creat Time: 2021/8/25
#更新内容：
#1.支持2019,2020,2021的所有*-vs-*-all.gene.xls输入
#2.支持p值以及q值筛选的选择
#3.修复了火山图左右两边有时候不对称，现输出图都是对称的。
#4.修复了画图过程中可能出现空值的情况
#部分修改指导：
#修改字体以及大小 修改theme(text=element_text())
#修改颜色修改 scale_colour_manual()

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
    make_option( c("-q", "--qval" ), type = "character",default = "p-value",
        help = "choose q-value or p-value as standard, default: p-value ."),
    make_option( c("-l", "--genelist" ), type = "character" ,
        help = "Genelist to display the gene symbol.  e.g. genelist.xls"),
    make_option(c("-t", "--title"), type = "character", default = NULL,
        help = "Graphic title and outputfile information: Group_A-vs-Group_B . "),
    make_option( c("-o", "--outputdir" ),type="character", default = "./Volcano",
        help="the output directory of Volcano results, default: ./Volcano ." ),
    make_option( c("-o", "--outputdir" ),type="character", default = "./Volcano",
        help="the output directory of Volcano results, default: ./Volcano ." )
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
    groupname <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", basename(opt$input)))#opt$input
    Case_names <- unlist(strsplit(groupname,"-vs-"))[1]
    Control_names <- unlist(strsplit(groupname,"-vs-"))[2]
    print (groupname)
}else {
    print_help(opt_parser)
    stop("diff file must be supplied, Volcano plot for groups will be skipped\n", call. = FALSE)
}

if ( is.null(opt$outputdir) ){
    output_dir = "Volcano"
}else{
    if ( file.exists(opt$outputdir) ){
        output_dir = opt$outputdir
    }else{
        output_dir = opt$outputdir
        dir.create(output_dir)
    }
}
print (output_dir)
if ( !is.null(opt$title) ){
    title = opt$title
}else {
    title = groupname
}

#==========Limitations==========
if ( !is.null(DEG$pValue)){
    DEG = plyr::rename(DEG, c("pValue"="pVAL"))
}
if ( !is.null(DEG$qValue) && opt$qval!="p-value" ){
    DEG = plyr::rename(DEG, c("qValue"="pValue"))
    mark="q-value"
}
if (!is.null(DEG$`q-value`) && opt$qval!="p-value" ) {
    DEG = plyr::rename(DEG, c("q-value"="pValue"))
    mark="q-value"
}
if ( !is.null(DEG$pval) && opt$qval=="p-value"){
    DEG = plyr::rename(DEG, c("pval"="pValue"))
    mark="p-value"
}
if ( !is.null(DEG$`p-value`) && opt$qval=="p-value"){
    DEG = plyr::rename(DEG, c("p-value"="pValue"))
    mark="p-value"
}
if ( !is.null(DEG$pVAL) && opt$qval=="p-value"){
    mark="p-value"
}

if ( !is.null(DEG$foldChange) ){
    DEG = plyr::rename(DEG, c("foldChange"="FoldChange"))
}


columnlist = unlist( strsplit( "log2FoldChange,pValue", ",", perl = T) )
for ( i in columnlist ){
    if ( !i %in% colnames(DEG) ){stop("NO specified column found!")}
}

l <- list(a = mark)
print (l)

eq <- substitute(-log[10]~a, l)

#==========Volcano function==========
#remove extremum
rownames(DEG)=DEG[,1]
DEG$pValue[which(DEG$pValue < 5E-300 )] = 5E-300

#replace the "-/Inf"
if( "Inf" %in% DEG$log2FoldChange | "-Inf" %in% DEG$log2FoldChange){
    tmp <- DEG[which(DEG$log2FoldChange != "Inf" & DEG$log2FoldChange != "-Inf"),]
    max = max(tmp$log2FoldChange)
    print (max)
    min = min(tmp$log2FoldChange)
    DEG$log2FoldChange <- as.numeric(sub("-Inf", min, DEG$log2FoldChange))
    DEG$log2FoldChange <- as.numeric(sub("Inf", max, DEG$log2FoldChange))
}

#去除NA
DEG[is.na(DEG)] <- 0
print (head(DEG))
#找到范围
x=ifelse(max(c(floor(max(DEG$log2FoldChange)),abs(floor(min(DEG[,"log2FoldChange"])))))>15,15,max(c(floor(max(DEG$log2FoldChange)),abs(floor(min(DEG[,"log2FoldChange"]))))))
print (xlim)
 if ( is.na(x) ){stop("请检查表格")}

#preparation by screening
logFC_cutoff = log(opt$foldchange, 2)
DEG$change = as.factor(ifelse(DEG$pValue < opt$pval & abs(DEG$log2FoldChange) > logFC_cutoff, ifelse(DEG$log2FoldChange > logFC_cutoff ,'Up','Down'),'Filtered'))
#add labels
if ( !is.null(opt$genelist) ) {
    labels <- read.delim(normalizePath(opt$genelist), header=T, sep="\t", check.names=F, quote="")
    print (labels[,1])
    print (head(DEG[,1]))
    DEG$label <-''
    DEG[match(labels[, 1], DEG[, 1]),]$label = as.character(labels[,1])
    #
    if (length(rownames(labels)) <= 15) {
        #sort data
        tmp1 <- DEG[which(DEG$label == ""),]
        tmp2 <- DEG[which(DEG$label != ""),]
        DEG <- rbind(tmp1,tmp2)
        #
        g = ggplot(data=DEG, aes(x=DEG$log2FoldChange, y=-log10(DEG$pValue), color=change)) + 
        geom_point( size=2) + 
        theme_set(theme_set(theme_bw(base_size=15))) + 
        theme(legend.title=element_blank()) + 
        labs(title = paste0(title, " :", mark,"< ", opt$pval," && |log2FC| > ",round(log(opt$foldchange, 2), 2) )) + 
        xlab(expression(paste(log[2], " Fold change"))) + 
        ylab(as.expression(eq)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(shape=guide_legend(override.aes=list(size=4))) +
        scale_colour_manual(values = c("Up"=c("#8B0000"), "Down"=c("darkblue"),"Filtered"=c("grey")),na.translate=FALSE) + 
        geom_text_repel(aes(label = label), size = 4, force = 1,color = "black") + 
        theme(text=element_text(size=14,family="ArialMT"))+
        scale_x_continuous(limits = c(-x,x))+
        theme(panel.grid = element_blank()) 
    }else {
        #sort data
        tmp1 <- DEG[which(DEG$label == ""),]
        tmp2 <- DEG[which(DEG$label != ""),]
        tmp2$change <- "Labels"
        DEG <- rbind(tmp1,tmp2)
        g = ggplot(data=DEG, aes(x=DEG$log2FoldChange, y=-log10(DEG$pValue), color=change)) + 
        geom_point( size=2) + 
        theme_set(theme_set(theme_bw(base_size=15))) + 
        theme(legend.title=element_blank()) + 
        labs(title = paste0(title, " :",mark, "< ",opt$pval," && |log2FC| > ", round(log(opt$foldchange, 2), 2) )) + 
        xlab(expression(paste(log[2], " Fold change"))) + 
        ylab(as.expression(eq)) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(shape=guide_legend(override.aes=list(size=4))) +
        theme(text=element_text(size=14,family="ArialMT"))+
        scale_x_continuous(limits = c(-x,x))+
        scale_colour_manual(values = c("Up"=c("#8B0000"), "Down"=c("darkblue"),"Filtered"=c("grey"),"Labels"=c("black")),na.translate=FALSE) + 
        theme(panel.grid = element_blank()) 
    }
}else {
    g = ggplot(data=DEG, aes(x=DEG$log2FoldChange, y=-log10(DEG$pValue), color=change)) + 
    geom_point(size=2) + 
    theme_set(theme_set(theme_bw(base_size=15))) + 
    theme(legend.title=element_blank()) + 
    labs(title = paste0(title, " :",mark, " < ",opt$pval," && |log2FC| > ", round(log(opt$foldchange, 2), 2) )) + 
    xlab(expression(paste(log[2], " Fold change"))) + 
    ylab(as.expression(eq)) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(shape=guide_legend(override.aes=list(size=4))) +
    theme(text=element_text(size=14,family="ArialMT"))+
    scale_x_continuous(limits = c(-x,x))+
    scale_colour_manual(values = c("Up"=c("#8B0000"), "Down"=c("darkblue"),"Filtered"=c("grey")),na.translate=FALSE) + theme(panel.grid = element_blank()) 
}

#lines
g = g + geom_hline(yintercept = -log(opt$pval, 10), linetype = "solid", color = c(oe_col_qua(8)), size = 0.5) + geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "solid", color = c(oe_col_qua(8)), size = 0.5)
#
ggsave(paste0(output_dir,"/",title,"-volcano-",mark,"-",opt$pval,"-FC-",opt$foldchange,".gene.pdf",sep=""),height=10,width=10,plot=g)
ggsave(paste0(output_dir,"/",title,"-volcano-",mark,"-",opt$pval,"-FC-",opt$foldchange,".gene.png",sep=""),height=10,width=10,plot=g,dpi=1000)
print("Volcano.pdf is OK")
