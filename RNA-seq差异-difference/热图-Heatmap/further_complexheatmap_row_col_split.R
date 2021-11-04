#==========import library==========
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))


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
  make_option(c("-f", "--fontsize"), type = "double",default = 5, 
              help = "Base fontsize for the plot, [default 5]", metavar = "double"),
  make_option(c("-p", "--phenotype"), type = "character", default = NULL,
              help = "Phenotype information of genes.  e.g. phenotype.xls. " ),
  make_option(c("-c", "--colors"), type = "character", default = "redwhiteblue",
              help = "colors choise for Heatmap picture : redwhiteblue, redblackgreen ,yellowblackblue , [default: redwhiteblue]. "),
  make_option(c("-g", "--height"), type = "double",  default = 8,
              help = "Height limit [default 8]", metavar = "double"),
  make_option(c("-k", "--width"), type = "double",  default = 6,
              help = "Width limit [default 6]", metavar = "double"),
  make_option(c("-o", "--outputdir"), type = "character", default = "./Heatmap",
              help = "output directory for Heatmap results ,[ default: ./Heatmap] . "),
  make_option(c("-u", "--cutrow"), type = "double", default = "1",
              help = "How many sections do you want? ,[ default:1] . ")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#==========parameter==========

#############读取表达文件###########
if (! is.null(opt$expression)) {
  expression <- read.table(normalizePath(opt$expression),row.names = 1, header=T, sep="\t", check.names=F, quote="")
}else {
  print_help(opt_parser)
  stop("expression matrix file is not supplied, heatmap for replicate groups  will be skipped\n", call. = FALSE)
}

##############标题#############
if ( !is.null(opt$title) ){
  title = opt$title
}else {
  title = "Heatmap"
}

##############角度#################
angle_list = unlist( strsplit( "0,45,90,270,315", ",", perl = T) )
if ( !opt$angle %in% angle_list ){stop("You can choose only from few predefined : 0, 45, 90, 270 and 315")}
if ( !is.null(opt$angle) ){
  angle = opt$angle
}else {
  angle = 90
}

##########颜色################
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


#######路径##########
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
#提取组名

groupname <- gsub("\\.(txt|xls)$", "", gsub("(-annotation.xls)$", "", basename(opt$expression)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Total$", "(Total)", groupname)
}

#==========Annotation================读取annotation文件
##############基因注释################
if ( !is.null(opt$phenotype) ){
  annotation_row <- read.table(normalizePath(opt$phenotype), header=T, sep="\t", row.names=1,check.names=F, quote="")
  print (annotation_row)
  t = row.names(expression)
  anno2=t(annotation_row)
  annotation_row= as.data.frame(t(subset(anno2,select = t)))
  print (annotation_row)
}else {
  annotation_row = NA
}
print (annotation_row)

#############组名注释##############
if ( !is.null(opt$group) ){
  annotation_col <- read.table(normalizePath(opt$group), header=T, sep="\t",row.names=1,check.names=F, quote="")
  print (annotation_col)
  t = colnames(expression)
  anno2=t(annotation_col)
  annotation_col= as.data.frame(t(subset(anno2,select = t)))
  print (annotation_col)
}else {
  annotation_col = NA
}
print (annotation_col)

palette <- colorRamp2(c(-1.5, 0, 1.5), c("Green", "Black", "Red"))

#palette <- colorRamp2(c(-1.5, 0, 1.5), c("steelblue", "white", "#D6675A"))




##############画图函数################

#注释
ann_colors = list(Regulation=c(Up="#D6675A",Down="steelblue"))
annotation_row_final = rowAnnotation(df = annotation_row, col = ann_colors)
annotation_col_final = HeatmapAnnotation(df = annotation_col)

#写出函数
getComplexHeatmap <- function (fpkm_matrix) {
  p=Heatmap(fpkm_matrix,show_column_names = TRUE,show_row_names = TRUE,name = "Color key", 
            row_names_gp =  gpar(fontsize = opt$fontsize, fontface="italic", col="black"),
            #clustering_distance_columns = "pearson",
            clustering_distance_rows = "euclidean",
            row_split = annotation_row$Regulation,
            column_split = annotation_col$Group,
            cluster_rows = opt$rowcluster,
            cluster_columns = T,
            top_annotation = annotation_col_final,
            left_annotation = annotation_row_final,
            #row_labels = "right",
            row_names_side = "right",
            #row
            column_title_rot = 0,
            row_title_rot = 0,
            column_dend_side = "top",
            row_dend_side = "left",
            #column_title = limit30(term_name), 
#            column_title = "gene_id",
            col = palette,
            #column_title_side = "Top",
            row_title_side = "left",
            column_title_gp = gpar(fontsize=8, fontface="bold"), 
            row_title_gp = gpar(fontsize=1, fontface="bold"),
            heatmap_legend_param=list(title= "Color Key", legend_direction="horizontal",title_rot = 90)#labels_rot = 0
  )
  return (p)
}

#################数据标准化处理###############################

#Force variance not to be zero！
ind <- apply(expression, 1, mean) > 0
expression <- expression[ind, ]
sd <- apply(expression, 1, sd)
if (length(sd[sd == 0])>0) {
  stop("sd = 0, program exit!", call. = FALSE)
}

#对数处理缩小数值范围，同时不会影响趋势
data<-log2(expression+0.0001)

#data_after_scale<-scale(data)

#################按行标准化，相当于pheatmap里的scale=“row”##############
# 批量按行中心标准化，减均值除方差，Z-score
mat_scaled = apply(data, 1, scale)
# 继续原数据表列名
rownames(mat_scaled) = colnames(data)
# 转置才与原方向一致
mat_scaled = t(mat_scaled)


##################出图#######################################


pdf(paste0(output_dir,"/",groupname,".",opt$rowcluster,"-output.pdf"),width=opt$width,height=opt$height, onefile = F)
getComplexHeatmap(mat_scaled)
dev.off()

