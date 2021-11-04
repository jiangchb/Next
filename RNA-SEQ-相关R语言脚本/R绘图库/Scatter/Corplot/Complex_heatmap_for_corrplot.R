library("corrplot")
library("ComplexHeatmap")
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))

path1 <- "AllSamples.correlation.xls"
path2 <- "sample_group24.txt"
d <- read.delim(path1, sep="\t",row.names = 1, header=T, quote="")
d2 <-  read.delim(path2, sep="\t",row.names = 1, header=T, quote="")

d_for_heatmap=as.matrix(d)
ann_colors = list(Group=c(BA="#006400",CC="#8B0000",CS="#FF8C00"))
annotation_row_final = rowAnnotation(df = d2, col = ann_colors)#
annotation_col_final = HeatmapAnnotation(df = d2,col = ann_colors)

palette <- colorRamp2(c(0.5, 0.75, 1), c("white", "blue", "darkblue"))

getComplexHeatmap <- function (fpkm_matrix) {
  p=Heatmap(fpkm_matrix,show_column_names = TRUE,show_row_names = TRUE,name = "Color key", 
            row_names_gp =  gpar(fontsize = 10, fontface="bold", col="black"),
            column_names_gp = gpar(fontsize = 10, fontface="bold", col="black"),
            #clustering_distance_columns = "pearson",
            clustering_distance_rows = "euclidean",
#            row_split = annotation_row$Regulation,
#            column_split = annotation_col$Group,
            cluster_rows = T,
            cluster_columns = T,
            top_annotation = annotation_col_final,
            left_annotation = annotation_row_final,
            #row_labels = "right",
            row_names_side = "right",
            column_names_side = "bottom",
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
            column_title_gp = gpar(fontsize=6, fontface="bold"), 
            row_title_gp = gpar(fontsize=6, fontface="bold"),
            heatmap_legend_param=list(title= "Color Key", legend_direction="horizontal",title_rot = 90)#labels_rot = 0
  )
  return (p)
}

getComplexHeatmap(d_for_heatmap)

p