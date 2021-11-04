#!/usr/bin/env Rscript
#############################################
# Author: Afan 
# Creat Time: 2019-12-17
# Modify :2020-09-24 Integration of two PCA tools
# Modify :2020-10-10 分组颜色问题
# Modify :2020-10-15 新增置信区间参数, ellipselevel 水平, ellipsecolor 透明度, 对错误参数填写生成warning文件
# Modify :2020-10-30 修改默认输出图片格式
# Modify :2020-12-11 添加字体参数-PCA-3D除外, 判断编码方式, 修改遍历数据类型的bug
# Modify :2020-12-26 判断编码方式优化, 空数据矩阵检查, 导出pca中间结果
# Modify :2021-01-21 修改报错日志文件名'oeweb_task.log'
#############################################


#==========parameter import==========
suppressPackageStartupMessages(library(optparse))
option_list = list(
    make_option(c("-i", "--input"), type = "character", default = NULL,
        help = "Input counts matrix file name(force).  e.g. counts.xls" ),
    make_option(c("-d", "--group"), type = "character", default = NULL,
        help = "Group information of samplenames, or use to reject samples.  e.g. sample_group.xls" ),
    make_option(c("-m", "--method"), type = "character", default = "default",
        help = "Methods for different input data : default,array,rnaseq,micro, [default: default]. " ),
    make_option(c("-l", "--ellipselevel"), type = "double", default = "",
        help = "Whether to add ellipse ,setting the confidence interval: 0.5-0.95. " ),
    make_option(c("-a", "--ellipsealpha"), type = "double", default = 0.2,
        help = "Whether to add ellipse ,setting the transparency of the circle background color : 0-1. " ),
    make_option(c("-f", "--format"), type = "character", default = "pdf",
        help = "Picture format type : pdf, jpg, png, tiff [default: pdf]. "),
    make_option(c("--fontfamily"), type = "character",  default = "Arial",
                help = "Fontfamily for picture: Arial, Times, Verdana  [default = Arial]"),
    make_option(c("--fontface"), type = "character",  default = NULL,
                help = "Fontface for picture: bold, italic  [default = NULL]")
    # make_option(c("-o", "--outputdir"), type = "character", default = ".",
        # help = "output directory for PCA results [default: . ]. ")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript  oeweb_PCA.R  -i counts.xls -m rnaseq -d  group.xls");
opt = parse_args(opt_parser);


#==========import library==========
options(rgl.useNULL= TRUE)
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(pca3d))
suppressPackageStartupMessages(library(maptools))
suppressPackageStartupMessages(library(oebio))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(grid))


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
    }else if( dim(data1)[2] < 2 ){
      stop("数据少于2列,未提取到需要的数据")
    }
    #### 删除行或列全为空 & 重复行 ####
    # 替换空字符
    data1[data1==''] <- NA
    # 删除空行和空列
    data1 <- data1[,colSums(is.na(data1)) != nrow(data1)]
    data1 <- data1[rowSums(is.na(data1)) != ncol(data1),]
    #### 数据矩阵和分组文件 ####
    if(checknames){
        if (any(duplicated(data1[,1])) & !is.null(data1[,1])) { stop(paste0(name," 行名输入有重复,请检查")) }
        if ("TRUE" %in% is.na(data1) | "TRUE" %in% is.na(names(data1)) ) { stop(paste0(name,"数据中存在空值,请检查数据的完整性")) }
        rownames(data1) <- data1[,1]
    }else{
        #### 检查行名 ####
        if( !is.numeric(data1[,2]) ){
            stop(paste0(name," 第二列非数据型,默认以第二列开始为数据输入,请根据分析工具页面上的使用说明对数据文件进行修改"))
        }
        if( typeof(data1[,1]) != "character" ){
            sink("oeweb_task.log", append=TRUE, split=FALSE)
            print("输入数据行名输入非字符串,默认以第二列开始为数据输入,请根据分析工具页面上的使用说明对数据文件进行修改 ")
            sink()
        }
        #### 过滤数据 #### 
        j = 1
        for(i in 2:length(data1))
        {   
          if( is.numeric(data1[,i]) ){j = j + 1}else break
        }
        data1 <- data1[,1:j]
        # 检查空值
        if ("TRUE" %in% is.na(data1) | "TRUE" %in% is.na(names(data1)) ) { stop(paste0(name,"数据中存在空值,请检查数据的完整性")) }
        # 去重
        # data1 <- data1[!duplicated(data1),]
        data1 <- data1[-1]
        if (dim(data1)[2] < 3){ stop("样本个数少于3,不进行主成分分析", call. = FALSE) }
    }

    return(data1)
}


#### draw PCA ####
draw_pca <- function(data, name=NA, method="default", group=NULL,
                    format=NA, level=NA, alpha=0.2, 
                    width=10, height=10, dpi=300,
                    family="Arial", fontface=NULL,
                    ...){
    
    data1 <- as.data.frame(data)
    
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
    
    if(is.null(group)){
        # select samples
        samples <- colnames(data1)
        groups <- as.character(colnames(data1))
        colData <- data.frame(row.names = samples, Group = groups)
    }else {
        # select samples
        data2 <- subset(data1, select = c(rownames(group)))
        samples <- colnames(data2)
        groups <- as.character(group$Group) 
        # show the group name  as same as  plot-3d  by reindex coldata and counts
        tmp_coldata <- data.frame(row.names = samples, Group = groups)
        colData <- tmp_coldata[order(tmp_coldata$Group), , drop = FALSE]
        data1 <- data2[ ,c(rownames(colData))]
        groups <- c(as.character(colData$Group))
        samples <- colnames(data1)
        # limit the biological duplication for ellipselevel
        for( i in unique(groups)){
            if (length(groups[which(groups == i)]) < 4 ){
                level = NA
                sink("oeweb_task.log", append=TRUE, split=FALSE)
                print("生物学重复少于等于3时，不添加置信区间")
                sink()
                break
            }
        }
    }
    
    # 不同平台过滤方法 default,array,rnaseq,micro
    if( method == "default" ){
        pca <- prcomp(t(data1), center = T, scale. = F)

    }else if( method=="array" ){
        # 对标准差过滤,阈值为 0.01
        ind <- apply(data1, 1, sd) > 0.01
        data1 <- data1[ind, ]
        # 主成分分析前,是否进行zscore
        # mat_scale <- apply(data1, 1, scale)
        pca <- prcomp(t(data1), center = T, scale. = T)

    }else if( method=="rnaseq" ){
        data1 <- round(data1)
        dds <- DESeqDataSetFromMatrix(countData = data1, colData = colData, design = ~ Group)
        # warnings()
        data <- rlog(dds, blind = T)
        # calculate the variance for each gene
        rv <- rowVars(assay(data))
        # select the ntop genes by variance, ntop=500
        select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
        pca <- prcomp(t(assay(data)[select,]), center = T, scale. = F)
        
    }else if( method=="micro" ){
        data1 <- data1[apply(data1,1,var) != 0,]
        pca <- prcomp(t(data1), center = T, scale. = T)

    }

    # plots
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
            scale_color_manual(values = oe_col_qua(1 : length(unique(pcaData$Group))),
            labels = c(as.character(unique(factor(pcaData$Group)))), name = "Group") +
            scale_fill_manual(values = oe_col_qua(1 : length(unique(pcaData$Group))),
            labels = c(as.character(unique(factor(pcaData$Group)))), name = "Group") +
            theme_bw() +
            theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
            theme(text = element_text(size = 15)) +
            scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
            scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
            geom_vline(xintercept = 0, linetype = 4, color = "grey") +
            geom_hline(yintercept = 0, linetype = 4, color = "grey") +
            coord_equal()
    
    if( !is.null(group) & !is.na(level) ){

        if( (level < 0.5 || level > 0.997)){
            level = 0.95
            sink("oeweb_task.log", append=TRUE, split=FALSE)
            print("置信区间可设置范围为0.5-0.997,默认0.95")
            sink()
        }
        if( (alpha < 0 || alpha > 1)){
            alpha = 0.2
            sink("oeweb_task.log", append=TRUE, split=FALSE)
            print("圈图背景色透明度可设置范围为0-1,默认0.2")
            sink()
        }
        pca1 = pca1 + stat_ellipse(aes(fill = Group, group = Group),geom = "polygon", level = level, alpha = alpha) + guides(fill=F)
    }

    # PCA 2D-2 
    pca2 = pca1 + geom_text_repel(aes(x = PC1, y = PC2, label = rownames(pcaData), fontface=fontface), size = 3, family=family)
    
    # force legend with 36 groups
    if (length(unique(groups)) <= 36){
        pca1 = pca1
        pca2 = pca2 
    }else{
        pca1 = pca1 + theme(legend.position='none')
        pca2 = pca2 + theme(legend.position='none')
    }
    
    pca1 = oe_ggplot_font(pca1, family=family, fontface=fontface)
    pca2 = oe_ggplot_font(pca2, family=family, fontface=fontface)

    if(is.na(name)){
        return(pca1)
        return(pca2)
    }else{
        if( format=="tiff" ){       
            ggsave(paste0("PCA_2D_1",name), pca1, dpi=dpi, width=width, height=height, compression ="zip" ,limitsize = F)
            ggsave(paste0("PCA_2D_2",name), pca2, dpi=dpi, width=width, height=height, compression ="zip" ,limitsize = F)
        }else if( format=="pdf" ){
            ggsave(paste0("PCA_2D_1",name), pca1, dpi=dpi, width=width, height=height, device=cairo_pdf, limitsize = F)
            ggsave(paste0("PCA_2D_2",name), pca2, dpi=dpi, width=width, height=height, device=cairo_pdf, limitsize = F)
        }else{
            ggsave(paste0("PCA_2D_1",name), pca1, dpi=dpi, width=width, height=height ,limitsize = F)
            ggsave(paste0("PCA_2D_2",name), pca2, dpi=dpi, width=width, height=height ,limitsize = F)
        }
    }
    
    # PCA 3D
    if (x > 1e-2 & y > 1e-2 & z > 1e-2){

        labels <- oe_col_qua(1 : length(seq_along(levels(pcaData$group))))

        all_col <- factor(pcaData$group, level = levels(pcaData$group), labels = labels)

    #### format ####
        if(is.na(format)){
        }else if(format=="tiff"){tiff(filename=paste0("PCA_3D",name), width=width, height=height, units="in", pointsize=12, res=dpi, compression="zip")
        }else if(format=="pdf"){pdf(file=paste0("PCA_3D",name), width=width, height=height, onefile = F)
        }else if(format=="jpg"){jpeg(filename=paste0("PCA_3D",name), width=width, height=height, units = "in", pointsize = 12, res=dpi)
        }else if(format=="png"){png(filename=paste0("PCA_3D",name), width=width, height=height, units = "in", pointsize = 12, res=dpi)}
        layout(matrix(1 : 2, 1, 2), width = c(7, 2))
        par(mar = c(5, 2, 5, 2))
        plot3d <- with(pcaData, scatterplot3d(PC1, PC2, PC3, xlab = paste("PC1 (", x, "%)"), ylab = paste("PC2 (", y, "%)") , zlab = paste("PC3 (", z, "%)"), main = "PCA 3D figure", color = as.character(all_col), pch = 16, label.tick.marks = TRUE, type = 'h', angle = 45))
        Position <- plot3d$xyz.convert(pcaData$PC1,pcaData$PC2,pcaData$PC3)
        pointLabel(Position$x, Position$y, labels = gsub("Sample_", "", rownames(pca$x)), cex = 0.8, col = as.character(all_col))
        par(mar = c(0.5, 0.5, 0.5, 0.3))
        plot.new()

        # force legend with 36 groups
        if (length(unique(groups)) <= 36){
            legend("center", plot3d$xyz.convert(2, 0.5, 2), pch = 16, xjust = - 1, yjust = - 1, legend = levels(pcaData$group), col = labels)
            dev.off()
        }else{
            dev.off()
        }
    }else{
        print(paste0(x, ",", y, ",", z))
        print("x,y,z轴坐标过小,无法生成PCA3D结果")
    }

    }


#==========parameter check==========
if (! is.null(opt$input)) {
    indata <- readdata(opt$input)
    if (is.null(dim(indata))){stop("输入矩阵数据为空！")}
    data <- dealname("输入文件",indata)
}else {
    stop("输入文件必须提供", call. = FALSE)
}
if (! is.null(opt$group)) {
    group <- readdata(opt$group)
    if( !"Group" %in% names(group) ){stop("Group表头不存在,请根据分析工具页面上的使用说明对数据文件进行修改")}
    ## 删除固定表头
    phenodata <- dealname("分组信息文件", group, checknames=T)
    if( dim(phenodata)[1] > dim(data)[2] ){stop("分组信息文件行数大于输入文件列数,请重新整理该文件")}
    if( "TRUE" %in% (!row.names(phenodata) %in% colnames(data)) ){stop("分组信息文件中样本名与数据矩阵表头不对应,请检查")}
}else {
    phenodata = NULL
}


#==========run function==========
format = unlist(strsplit(x=opt$format,split=","))
for(j in 1:length(format)){
    draw_pca(data=data, name=ifelse(is.na(format[j]),format[j],paste(".",format[j],sep="")), method=opt$method, level=opt$ellipselevel, alpha=opt$ellipsealpha, group=phenodata, format=format[j], family=opt$fontfamily, fontface=opt$fontface)
}

if ( file.exists("Rplots.pdf") ){ file.remove("Rplots.pdf") }

