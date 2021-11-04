#!/usr/bin/env Rscript
#############################################
#Author: Afan 
#Creat Time: 2019-12-17
#Modify :2020-06-09 Add reading xlsx,txt,xls,csv fomat files
# Modify :2020-08-30 Integration of three volcano tools
# Modify :2020-09-18 Add the arraylist parameter,fix the erro for "#NAME?" and NA, lable color
# Modify :2020-10-30 修改默认输出图片格式
# Modify :2020-11-04 添加标注特征是否存在相关判断,删除未知特征并输出warning文件
# Modify :2020-12-11 添加字体参数, 判断编码方式, 修改log2FC数据类型报错的bug
# Modify :2020-12-26 判断编码方式优化, 空数据矩阵检查, 添加展示基因数目参数默认为30
# Modify :2021-01-21 空文本框输入报错修复, 修改数据检查逻辑; 修改报错日志文件名'oeweb_task.log'
# Modify :2021-02-01 增加点大小参数
# Modify :2021-06-17 按照常规流程修改工具样式及表头
#############################################


#==========parameter import==========
suppressPackageStartupMessages(library(optparse))
option_list = list(
    make_option( c("-i", "--input" ), type = "character",
        help = "The input file which differences not filtered (force).  e.g. *-vs-*-all.xls" ),
    make_option( c("-p", "--pval" ), type = "double",default = 0.05,
        help = "The number of p-Value/q-value ratio threshold, [default: 0.05] ."),
    make_option( c("-f", "--foldchange" ), type = "double",default = 2,
        help = "FoldChange threshold, [default: 2] ."),
    make_option(c("-t", "--format"), type = "character", default = "pdf",
        help = "Picture format type : pdf, jpg, png, tiff [default: pdf]. "),
    make_option( c("-q", "--q_ptype" ), type = "character",default = "p-value",
        help = "q-value or p-value, default: p-value."),
    make_option(c("-c", "--colors"), type = "character", default = "GGR",
        help = "Colors choise for volcano picture : GGR, BGR, [default: GGR]. "),
    make_option( c("-g", "--genelist" ), type = "character" , default = NULL,
        help = "Genelist to display the gene symbol.  e.g. genelist.txt"),
    make_option( c("-m", "--matchlist" ), type = "character" , default = NULL,
        help = "Matchlist to display the probe gene symbol.  e.g. probelist.xls"),
    make_option(c("-v", "--vip"), type = "double", default = 0,
        help = "Whether to filter datas by VIP, 0 for no visualization. "),
    make_option( c("-a", "--arraylist" ), type = "character" , default = NULL,
        help = "Whether to filter datas by DEG for array.  e.g. array_volcano_deg_data.csv"),
    make_option(c("--fontfamily"), type = "character",  default = "Arial",
                help = "Fontfamily for picture: Arial, Times, Verdana  [default = Arial]"),
    make_option(c("--fontface"), type = "character",  default = NULL,
                help = "Fontface for picture: bold, italic  [default = NULL]"),
    make_option( c("-l", "--limit" ), type = "integer" , default = 30,
        help = "Limit the gene symbol nums to display. [default = 30]"),
    make_option( c("--pointsize"), type = "double" , default = 1,
        help = "Point size, cex 1. [default = 1]")
    # make_option( c("-o", "--outputdir" ),type="character", default = ".",
        # help="the output directory of Volcano results, [default: . ]." )
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript  merge_volcano.R  -i all.xls");
opt = parse_args(opt_parser);


#==========import library==========
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(oebio))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(plyr))


#==========functions==========
#### deal with reading xlsx,xls,csv,txt file ####
readdata <- function(name, sheet = 1, rowNames = F, skip = 0, comment="", header = T, silent = F, errorstop = T){
  data <- NULL
  # 读取文本文件
  rcsv = function(csvfile, sep="\t", encoding = "UTF-8", quote=""){
      data <- read.csv(file = csvfile, header = header, check.names = F, encoding = encoding, stringsAsFactors = F, skip = skip, row.names={if(rowNames){1}else{NULL}}, sep = sep, quote = quote, blank.lines.skip=TRUE, na.strings = "NA", strip.white=T, comment.char = comment)
  }
  # 猜测文件的字符编码格式
  file_enc = function(csvfile){
      enc <- guess_encoding(csvfile, n_max = 1000)[1,1, drop = TRUE]
      # UTF-8 格式有warning时报错
      if ( !is.na(enc) & (enc == "UTF-8" | enc == "ASCII") ) {
          rcsv(csvfile)
      }else{
          stop(paste0(csvfile," 文件编码方式非 UTF-8 ！"))
      }
  }

  data_try <- try({
    if (is.character(name)) {
      if (grepl(".txt$",name)) {
          data <- file_enc(name)
      }else if (grepl(".csv$",name)) {
          data <- rcsv(name, sep = ",", quote="\"")
      }else if (grepl(".xlsx$",name)) {
          data <- read.xlsx(xlsxFile = name, sheet = sheet, rowNames = rowNames, startRow = 1, colNames = header, na.strings= "NA", skipEmptyRows=T , skipEmptyCols=T)
      }else if (grepl(".xls$",name)) {
          rxls <- try(data <- read_xls(path = name, sheet = sheet, col_names = header, skip = skip, trim_ws=T), silent = T)
          if("try-error" %in% class(rxls) & errorstop){
              data <- file_enc(name)
          }
      }else { stop(paste0(name,"非xlsx,csv,xls,txt文件,请重新整理文件"))}
    }
  }, silent = silent)
  if ("try-error" %in% class(data_try) ) {
    stop(paste0(name,"非xlsx,csv,xls,txt文件,请重新整理文件"))
  }
  
  return(data)
}

#### dealname : null ids , special characters, limit the header ####
dealname <- function(name, data, checknames=F){
    data1 <- as.data.frame(data)
    #### 不规范输入 ####
    if( dim(data1)[1] < 1 ){
        stop("数据少于2行,未提取到需要的数据")
    }
    #### 删除行或列全为空 & 重复行 ####
    # 替换空字符
    data1[data1==''] <- NA
    #### 检查行名 ####
    if( !is.character(data1[,1]) ){
        stop(paste0(name," 第一列输入非字符串,请根据分析工具页面上的使用说明对数据文件进行修改"))
    }
    # 删除空列和空行
    data1 <- as.data.frame(data1[,colSums(is.na(data1)) != nrow(data1)])
    data1 <- as.data.frame(data1[rowSums(is.na(data1)) != ncol(data1),])
    # 检查只有表头没有数据的文件格式
    if( dim(data1)[2] == 0 ){
        stop(paste0(name," 数据为空,请根据分析工具页面上的使用说明对数据文件进行修改"))
    }
    # 去重
    data1 <- as.data.frame(data1[!duplicated(data1),])
    names(data1)[1] <- names(data)[1]
    if( "TRUE" %in% is.na(data1[,1]) ){
        stop(paste0(name," 行名有空值,请根据分析工具页面上的使用说明对数据文件进行修改"))
    }
    #### 输入文件及标注及差异 ####
    if(checknames){
        if ("TRUE" %in% is.na(data1) | "TRUE" %in% is.na(names(data1)) ) { stop(paste0(name," 数据中存在空值,请检查数据的完整性")) }
    }else{
        # 字符型转换为数值型
        data1$log2FoldChange <- as.numeric(data1$log2FoldChange)
        data1$pValue <- as.numeric(data1$pValue)
        # remove NA
        data1 <- data1[which(data1$pValue !="" ),]
        # remove extremum
        data1$pValue[which(data1$pValue < 5E-300 )] = 5E-300
        # replace the "-/Inf" and NA for excel format
        if( ("Inf" %in% data1$log2FoldChange) || ("-Inf" %in% data1$log2FoldChange) || "TRUE" %in% is.na(data1$log2FoldChange)){
            tmp <- data1[which(data1$log2FoldChange != "Inf" & data1$log2FoldChange != "-Inf" & !is.na(data1$log2FoldChange)),]
            max = max(tmp$log2FoldChange)
            min = min(tmp$log2FoldChange)
            data1$log2FoldChange <- as.numeric(sub("-Inf", min, data1$log2FoldChange))
            data1$log2FoldChange <- as.numeric(sub("Inf", max, data1$log2FoldChange))
            if("TRUE" %in% is.na(data1$log2FoldChange)){data1$log2FoldChange[which(is.na(data1$log2FoldChange))] = min}
        }
    }
    #### 特殊id ####
    ids <- gsub(pattern = "γ",replacement ="gamma" ,x = data1[,1])
    ids <- gsub(pattern = "α",replacement ="alpha" ,x = ids)
    ids <- gsub(pattern = "β",replacement ="beta" ,x = ids)
    #### 日期型基因 ####
    ids <- gsub(pattern = "43900",replacement ="MARCH10" ,x = ids)
    ids <- gsub(pattern = "43901",replacement ="MARCH11" ,x = ids)
    ids <- gsub(pattern = "43893",replacement ="MARCH3" ,x = ids)
    ids <- gsub(pattern = "43894",replacement ="MARCH4" ,x = ids)
    ids <- gsub(pattern = "43895",replacement ="MARCH5" ,x = ids)
    ids <- gsub(pattern = "43896",replacement ="MARCH6" ,x = ids)
    ids <- gsub(pattern = "43897",replacement ="MARCH7" ,x = ids)
    ids <- gsub(pattern = "43898",replacement ="MARCH8" ,x = ids)
    ids <- gsub(pattern = "43899",replacement ="MARCH9" ,x = ids)
    ids <- gsub(pattern = "44075",replacement ="SEPT1" ,x = ids)
    ids <- gsub(pattern = "44084",replacement ="SEPT10" ,x = ids)
    ids <- gsub(pattern = "44085",replacement ="SEPT11" ,x = ids)
    ids <- gsub(pattern = "44086",replacement ="SEPT12" ,x = ids)
    ids <- gsub(pattern = "44088",replacement ="SEPT14" ,x = ids)
    ids <- gsub(pattern = "44076",replacement ="SEPT2" ,x = ids)
    ids <- gsub(pattern = "44077",replacement ="SEPT3" ,x = ids)
    ids <- gsub(pattern = "44078",replacement ="SEPT4" ,x = ids)
    ids <- gsub(pattern = "44079",replacement ="SEPT5" ,x = ids)
    ids <- gsub(pattern = "44080",replacement ="SEPT6" ,x = ids)
    ids <- gsub(pattern = "44081",replacement ="SEPT7" ,x = ids)
    ids <- gsub(pattern = "44082",replacement ="SEPT8" ,x = ids)
    ids <- gsub(pattern = "44083",replacement ="SEPT9" ,x = ids)
    ids <- strsplit(x = ids, split = " /")
    i <- 1
    while(i<=length(data1[,1])){
        data1[,1][i] <- ids[[i]][1]
        i <- i+1
    }

    return(data1)
}

#### prepare data for volcano ####
pre_volcano <- function(data, pval=0.05, FC=NA, VIP=NA, arrayid=NULL, shownames=NA,limit=30,
                        ...){
    data1 <- data

    if(is.na(FC)){logFC_cutoff = 0}else{logFC_cutoff = log2(FC)}
    data1$class = as.factor(ifelse(data1["pValue"] < pval & abs(data1["log2FoldChange"]) > logFC_cutoff, 
                            ifelse(data1["log2FoldChange"] > logFC_cutoff ,"Up","Down"),"Filtered"))

    #### filter VIP ####
    if(is.na(VIP)){

    }else {
        if ( !"VIP" %in% colnames(data1) ){stop("VIP表头不存在,请检查文件表头是否正确")}
        data1$class[data1$VIP < VIP] <- "Filtered"
        data1 <- data1[order(data1$VIP),]
    }

    #### filter array non-DEG ####
    if(is.null(arrayid)){

    }else {
        extra_gene = data1[,1]
        filter_gene = c()
        for ( gene in extra_gene ){
            if(is.na(match(gene,arrayid[,1]))){
                extra_gene = extra_gene[-which( extra_gene == gene )]
                filter_gene = c(filter_gene,gene)
            }
        }
        data1[match(filter_gene,data1[,1]),]$class <- "Filtered"
    }

    tmp1 <- data1[which(data1$class == "Filtered"),]
    tmp2 <- data1[which(data1$class == "Up" | data1$class == "Down"),]
    data1 <- rbind(tmp1,tmp2)

    #### add labels ####
    if(!is.na(shownames) & shownames > limit ){
        tmp1 <- data1[which(data1$label == ""),]
        tmp2 <- data1[which(data1$label != ""),]
        tmp2$class <- "Labels"
        data1 <- rbind(tmp1,tmp2)
    }

    return(data1)
}


#### draw volcano ####
draw_volcano <- function(data, name=NA, format=NA, pval=0.05, FC=NA,
                        x=ifelse(max(c(floor(max(data[,"log2FoldChange"])),abs(floor(min(data[,"log2FoldChange"])))))+5>15,15,max(c(floor(max(data[,"log2FoldChange"])),abs(floor(min(data[,"log2FoldChange"])))))+5),
                        y=floor(-log10(min(data[,"pValue"])))+2,
                        cex=1, aspect.ratio=1, size=c(0.5,3),
                        color=c("blue", "grey", "red"),
                        plot.margin=unit(c(0.3,0.3,0.3,0.3),"in"),
                        legend.position="right",
                        legend.background=element_rect(colour="black",size=0.3),
                        VIP=NA, shownames=NA,
                        width=7, height=6, dpi=300, other=theme(), 
                        family="Arial", fontface=NULL,limit=30,
                        ...){
  data1 <- data
  
  # oeRtools function
  oe_ggplot_font <- function(p, family="Arial", fontface=NULL) {
    par <- list(fontfamily = family, fontface = fontface)
    par <- par[!sapply(par, is.null)]
    gp <- do.call(gpar, par)
    g <- ggplotGrob(p)
    ng <- grid.ls(grid.force(g), print=FALSE)$name
    txt <- ng[which(grepl("text", ng))]
    
    for (i in seq_along(txt)) {
      g <- editGrob(grid.force(g), gPath(txt[i]),
                    grep = TRUE, gp = gp)
    }
    return(g)
  }

  #### filter VIP ####
  if(is.na(VIP)){
      pp <- ggplot(data=data1, 
                  mapping=aes(x =log2FoldChange, y=-log10(pValue), color=class, text=paste0("ID:",row.names(data1))) )+
                  geom_point(size=cex) + guides(color=guide_legend(title=NULL))
  }else {
      pp <- ggplot(data=data1,
                  mapping=aes(x=log2FoldChange, y=-log10(pValue), color=class, size=VIP, text=paste0("ID:",row.names(data1)))) +
                  geom_point() + scale_size_continuous(range=size) + 
                  guides(color=guide_legend(title=NULL)) + guides(size=guide_legend(order = 1))
  }
  #### add labels ####
  if(is.na(shownames)){
      pp <- pp + scale_color_manual(values = c("Filtered"=color[2],"Up"=color[3],"Down"=color[1]), limits=c("Filtered", "Up", "Down"))

  }else if(shownames <= limit ){
      pp <- pp + scale_color_manual(values = c("Filtered"=color[2],"Up"=color[3],"Down"=color[1]), limits=c("Filtered", "Up", "Down")) + 
                 geom_text_repel(aes(x =log2FoldChange, y=-log10(pValue), label = label, fontface=fontface), size = 1.8, segment.size=0.05, color = "black", family=family)
#  vjust=0.5,
  }else if(shownames > limit){
      pp <- pp + scale_color_manual(values = c("Filtered"=color[2],"Up"=color[3],"Down"=color[1],"Labels"=c(oe_col_qua(5))), limits=c("Filtered", "Up", "Down", "Labels"))
  }
  if(opt$q_ptype == "p-value"){
             pp <- pp + scale_x_continuous(limits = c(-x,x)) + scale_y_continuous(limits = c(0,y)) + theme_bw() +
             labs(x=expression(paste(log[2], " Fold Change")),y=expression(paste("-", log[10], " p-value")), title = paste0("Volcano Plot : p-value < ", pval ," && |log2FC| > ",round(log(FC, 2), 2) )) + 
             theme(panel.grid=element_blank(),
                  plot.title=element_text(hjust = 0.5),
                  aspect.ratio=aspect.ratio,
                  plot.margin=plot.margin,
                  legend.position=legend.position,
                  legend.background=legend.background,
                  line=element_line(colour = "black"),
                  legend.margin=margin(t = 0, r = 6, b = 6, l = 2, unit = "pt"))
   }else {             
         pp <- pp + scale_x_continuous(limits = c(-x,x)) + scale_y_continuous(limits = c(0,y)) + theme_bw() +
         labs(x=expression(paste(log[2], " Fold Change")),y=expression(paste("-", log[10], " q-value")), title = paste0("Volcano Plot : q-value < ", pval ," && |log2FC| > ",round(log(FC, 2), 2)))+ 
             theme(panel.grid=element_blank(),
                  plot.title=element_text(hjust = 0.5),
                  aspect.ratio=aspect.ratio,
                  plot.margin=plot.margin,
                  legend.position=legend.position,
                  legend.background=legend.background,
                  line=element_line(colour = "black"),
                  legend.margin=margin(t = 0, r = 6, b = 6, l = 2, unit = "pt"))
   }
  if(is.na(pval)){

  }else{
      pp <- pp + geom_hline(yintercept=c(-log10(pval)), size=0.2, colour=oe_col_qua(8))
  }

  if(is.na(FC)){

  }else {
      pp <- pp + geom_vline(xintercept=c(log2(FC)), size=0.2, colour=oe_col_qua(8))+
                 geom_vline(xintercept=c(-log2(FC)), size=0.2, colour=oe_col_qua(8))
  }

  pp <- pp + other 
  pp = oe_ggplot_font(pp, family=family, fontface=fontface)

  if(is.na(name)){
      return(pp)
  }else{
    if( format=="tiff" ){
      ggsave(name,pp,dpi=dpi,width=width,height=height,compression ="zip",limitsize = F)
    }else if( format=="pdf" ){
      ggsave(name,pp,dpi=dpi,width=width,height=height,device=cairo_pdf,limitsize = F)
    }else{
      ggsave(name,pp,dpi=dpi,width=width,height=height,limitsize = F)}
  }
}


#==========parameter check==========
if (!is.null(opt$input)) {
    all <- readdata(opt$input)
    if (is.null(dim(all))){stop("输入文件数据为空！")}
    # 表头替换
    if (opt$q_ptype == "p-value") {
            mark="p-value"
            if ( !is.null(all[["pval"]]) ){
            all = plyr::rename(all, c("pval"="p-value"))
        }
            if ( !is.null(all[["P-value"]]) ){
            all = plyr::rename(all, c("P-value"="p-value"))
        }
            if ( !is.null(all[["pValue"]]) ){
            all = plyr::rename(all, c("pValue"="p-value"))
        }
    }else {
            mark="q-value"
            if ( !is.null(all$qValue)){
            all = plyr::rename(all, c("qValue"="p-value"))
            }
            if ( !is.null(all$`q-value`)) {
            all = plyr::rename(all, c("q-value"="p-value"))
            }
    }
    if ( !is.null(all[["log2(FC)"]]) ){
    all = plyr::rename(all, c("log2(FC)"="log2FoldChange"))
    }
    if ( !is.null(all[["log2FC"]]) ){
    all = plyr::rename(all, c("log2FC"="log2FoldChange"))
    }
    columnlist = unlist( strsplit( "Names,log2FoldChange,p-value", ",", perl = T) )
    for ( i in columnlist ){if ( !i %in% colnames(all) ){stop(i," 表头不存在,请根据分析工具页面上的使用说明对-输入文件-进行修改")}}
    all = plyr::rename(all, c("p-value"="pValue"))
    data <- dealname("输入文件", all, checknames=F)
}else {
    stop("输入文件必须提供", call. = FALSE)
}
# -m 和 -g 不能同时使用,但由于-g 为oeweb文本框输出文件,需使用判断区分
if (!is.null(opt$matchlist)){
    genelist <- readdata(opt$matchlist)
    if ( !"IDs" %in% colnames(genelist) ){stop("IDs表头不存在,请根据分析工具页面上的使用说明对-标记名称列表-进行修改")}
    labels <- dealname("标记名称列表", genelist, checknames=T)
}else {
    genelist <- readdata(opt$genelist)
    # 默认值 "IDs" 即文本框表头,生成txt文件,无需处理空行 空格等问题
    if( dim(genelist)[1] != 0 ){
        #### 检查行名 ####
        if( !is.character(genelist[,1]) ){
            stop("输入为非字符串类型,请根据分析工具页面上的使用说明对数据进行修改")
        }
        # 去重
        labels <- as.data.frame(genelist[!duplicated(genelist),])
        names(labels)[1] <- names(genelist)[1]
    }else{
        labels <- genelist  
    }
}
if( dim(labels)[1] != 0 ){
    if ( "TRUE" %in% (!labels$IDs %in% data[,1]) ) {
        filter <- as.data.frame(labels[!labels$IDs %in% data[,1],])[,1]
        labels <- as.data.frame(labels[labels$IDs %in% data[,1],])
        names(labels)[1] <- "IDs"
        sink("oeweb_task.log", append=TRUE, split=FALSE)
        print(paste0(filter," 标记名称列表中特征名与输入文件特征名不对应,请检查"))
        sink()
    }
    data$label <-''
    if (length(names(labels))==1){
        data[match(labels$IDs,data[,1]),]$label = as.character(labels$IDs)
    }else if(length(names(labels))==2){
        data[match(labels$IDs,data[,1]),]$label = as.character(labels[,2])
    }else {
        stop("标记名称列表文件格式不正确,请检查")
    }
    shownames <- length(rownames(data[which(data$label != ""),]))
}else {
    shownames <- NA
}
# 芯片必须
if (!is.null(opt$arraylist)){
    arraylist <- readdata(opt$arraylist)
    if ( !"Names" %in% colnames(arraylist) ){stop("Names表头不存在,请检查芯片差异基因文件表头是否正确")}
    arrayid <- dealname("芯片差异基因文件", arraylist, checknames=T)
}else {
    arrayid <- NULL
}
if ( opt$colors == "GGR" ){palette <- c("#42B540FF-#ADB6B6FF-#ED0000FF")
}else if ( opt$colors == "BGR" ){palette <- c("#00468BFF-#ADB6B6FF-#ED0000FF")
}
if( !is.na(opt$pointsize) ){
    if(opt$pointsize < 1 || opt$pointsize > 3){
        cex = 1
        sink("oeweb_task.log", append=TRUE, split=FALSE)
        print("点大小设置不正确，默认为1,取值范围在1-3之间")
        sink()
    }
    else {
        cex = opt$pointsize
    }
}


#==========run function==========
format = unlist(strsplit(x=opt$format,split=","))
for(j in 1:length(format)){
    if(opt$vip==0){ VIP <- NA }else { VIP <- opt$vip }
    pre_data = pre_volcano(data, pval=opt$pval, FC=opt$foldchange, VIP=VIP, shownames=shownames, arrayid=arrayid, limit=opt$limit)
    draw_volcano(data=pre_data, name=ifelse(is.na(format[j]),format[j],paste("volcano.",format[j],sep="")),
                pval=opt$pval, FC=opt$foldchange, cex=cex,
                color=unlist(strsplit(x = palette,split = "-")), 
                VIP=VIP, shownames=shownames, format=format[j], family=opt$fontfamily, fontface=opt$fontface, limit=opt$limit)
  }

if ( file.exists("Rplots.pdf") ){ file.remove("Rplots.pdf") }


