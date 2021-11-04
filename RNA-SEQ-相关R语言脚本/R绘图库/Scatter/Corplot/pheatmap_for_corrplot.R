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

annotation_col=read.delim(path2, header=T, sep="\t", row.names=1,check.names=F, quote="")
annotation_row=read.delim(path2, header=T, sep="\t", row.names=1,check.names=F, quote="")

palette <- colorRampPalette(c("white","darkBlue"))(n=256)

getperfectHeatmap <- function (fpkm_matrix,color1) {
  p=pheatmap(fpkm_matrix,
          annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors=ann_colors,
           treeheight_row=30, treeheight_col=30, lwd=1,
           color = palette, show_rownames = T, fontsize=7, border= T,
           cluster_rows=T , cluster_cols =T , angle_col = 90,border_color = "white",cutree_rows=1,display_numbers = T,number_color = color1
  )
  return (p)
}


pdf(paste0("./","corr_plot-output.pdf"),width=8,height=6, onefile = F)
getperfectHeatmap(d,"white")
dev.off()
