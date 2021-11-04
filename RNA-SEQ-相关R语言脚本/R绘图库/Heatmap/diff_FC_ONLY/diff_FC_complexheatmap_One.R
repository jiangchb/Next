suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))

#suppressPackageStartupMessages(library(optparse))
#option_list = list(
#  make_option( c("-i", "--input" ), type = "character",
#               help = "The input differential gene file(force).  e.g. *-vs-*-all.gene.xls" ),
#  make_option( c("-f", "--fontsize" ), type = "double",default = 8,
 #              help = "pValue ratio threshold, default: 8 .")
#);

#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);

path1="heatmap.reorder_cluster_result.xls"
Three <- read.delim(normalizePath(path1), header=T, sep="\t",row.names = 1, check.names=F, quote="")

Three <- Three[order(Three[,"LavendustinA-vs-DMSO"],decreasing = T),]

getComplexHeatmap <- function (fpkm_matrix,palette,title_name,legend_name,width,row_fontsize) {
  p=Heatmap(fpkm_matrix,show_column_names = TRUE,show_row_names = TRUE, 
            row_names_gp =  gpar(fontsize = row_fontsize, fontface="italic", col="black"),
            #clustering_distance_columns = "pearson",
            clustering_distance_rows = "euclidean",
            #            row_split = annotation_row$Regulation,
            #            column_split = annotation_col$Group,
            #            cluster_rows = opt$rowcluster,
            cluster_columns = F,
            cluster_rows = F,
            #            top_annotation = annotation_col_final,
            #            left_annotation = annotation_row_final,
            #row_labels = "right",
            row_names_side = "left",
            width = unit(width, "cm"),
            column_names_rot = 45,
            #row
            column_title_rot = 0,
            row_title_rot = 0,
            column_dend_side = "top",
            row_dend_side = "left",
            column_title = title_name, 
            #            column_title = "gene_id",
            col = palette,
            #column_title_side = "Top",
            row_title_side = "left",
            rect_gp = gpar(col = "white"),
            column_title_gp = gpar(fontsize=8, fontface="bold"), 
            row_title_gp = gpar(fontsize=1, fontface="bold"),
            heatmap_legend_param=list(title= legend_name,legend_direction="vertical",title_rot = 90)#labels_rot = 0
  )
  return (p)
}

palette <- colorRamp2(c(-2, 0, 2), c("RoyalBlue2", "white", "Red2"))

p1=getComplexHeatmap(as.matrix(Three),palette, "Heatmap" ,"ColorKey",3,9)
p1
pdf(paste0("./","diff_FC_SIG",".","output.pdf"),width=8,height=8, onefile = F)
draw(p1)
dev.off()
