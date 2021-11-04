#!/usr/bin/env Rscript
usage = "\
usage:
	*replicate* : 	Rscript pca_normal.r  -i  [count matrix file] -s [phenodata matrix file]
	*noreplicate*:	Rscript pca_normal.r  -i [count matrix file]
			Rscript pca_normal.r  -i [count matrix file]  -l [select sample list ]			
input file format:
	*phe.xls*:
	Sample	Group
	Sample_A	group2
	Sample_B	group2
	*select_sample_list.xls*:
	Sample
	Sample_A
	Sample_B
"
Sys.time()
cat(usage)
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggrepel))
set.seed(42)
option_list = list(
make_option(c("-i", "--readcounts"), type = "character", default = NULL,
help = "read counts matrix file name(force)", metavar = "character"),
make_option(c("-s", "--group"), type = "character", default = NULL,
help = "sample grouping  information", metavar = "character"),
make_option(c("-l", "--list"), type = "character", default = NULL,
help = "select samples list", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$readcounts)) {
    print_help(opt_parser)
    stop("read counts matrix file must be supplied", call. = FALSE)
}

###PCA.pdf  for replicate samples
rep_pca_plot <- function(counts, group){
    counts <- round(read.delim(counts, row.names = 1, check.names = F,quote=""))
    phenodata <- read.table(group, row.names = 1, header = T, sep = "\t", check.names = F)
###################################################################################################
    counts_tmp<-subset(counts, select = c(rownames(phenodata)))
	counts<-counts_tmp
###################################################################################################
    samples <- colnames(counts)
	index_phenodata <- phenodata
	for (i in 1 : length(samples))
	{
		index_phenodata[i,] <- phenodata[samples[i],]
	}
	print(index_phenodata)
	##modified 20170911 ,fixed the group name showing error by define the  level of index_phenodata  
	rownames(index_phenodata)<-c(as.character(samples))
	index_phenodata$Group<-factor(index_phenodata$Group,level=unique(index_phenodata$Group))
	groups <- as.character(index_phenodata$Group) 

	##modified 20170911 ,show the group name  as same as  plot-3d  by reindex coldata and counts 
	tmp_coldata <- data.frame(row.names = samples, Group = groups)
	colData<-tmp_coldata[order(tmp_coldata$Group), , drop = FALSE]
	counts<-counts[ ,c(rownames(colData))]
	groups <-c(as.character(colData$Group))
    samples <- colnames(counts)
    print(head(counts))
    if (length(samples) < 3)
    {cat("Samples number is less than 3,stop!!", call. = FALSE)}
    else
    {

        ##########################################################################
        dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Group)
        warnings()
        ##PCA 2D
        rld <- rlog(dds, blind = T)
        pcaData <- plotPCA(rld, intgroup = c("Group"), returnData = TRUE)
        min <- round(min(pcaData$PC1, pcaData$PC2))
        max <- ceiling(max(pcaData$PC1, pcaData$PC2))
		print(attr(pcaData, "percentVar"))
        percentVar <-  round(100 *attr(pcaData, "percentVar"),2)	
        ggplot(pcaData, aes(PC1, PC2, color = Group)) +
            geom_point() +
            xlab(paste0("PC1: ", round(percentVar[1],2), "% variance")) +
            ylab(paste0("PC2: ", round(percentVar[2],2), "% variance")) +
            guides(col = guide_legend(nrow = 15)) +
            scale_color_manual(values = oe_col_qua(1 : length(seq_along(levels(pcaData$Group)))),
            labels = c(as.character(unique(pcaData$Group))),
            name = "Group") +
            scale_fill_manual(values = oe_col_qua(1 : length(seq_along(levels(pcaData$Group)))),
            labels = c(as.character(unique(pcaData$Group))),
            name = "Group") +
            theme_bw() +
            theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
            theme(text = element_text(size = 15, family = "ArialMT")) +
            scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
            scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
            geom_point(size = 4) +
            geom_vline(xintercept = 0, linetype = 4, color = "grey") +
            geom_hline(yintercept = 0, linetype = 4, color = "grey") +
            coord_equal()          
        #	ggsave("PCA_deseq2_replicate.pdf",height=4*0.05*percentVar[2],width=4*0.05*percentVar[1])
        ggsave("PCA_1.pdf", height = 8, width = 8)
        ggsave("PCA_1.png",type="cairo-png",height=8,width=8)

        ggplot(pcaData, aes(PC1, PC2, color = Group)) +
            geom_point() +
            xlab(paste0("PC1: ", round(percentVar[1],2), "% variance")) +
            ylab(paste0("PC2: ", round(percentVar[2],2), "% variance")) +
            guides(col = guide_legend(nrow = 15)) +
            scale_color_manual(values = oe_col_qua(1 : length(unique(pcaData$Group))),
            labels = c(as.character(unique(factor(pcaData$Group)))),
            name = "Group") +
            scale_fill_manual(values = oe_col_qua(1 : length(unique(pcaData$Group))),
            labels = c(as.character(unique(factor(pcaData$Group)))),
            name = "Group") +
            theme_bw() +
            theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
            theme(text = element_text(size = 15, family = "ArialMT")) +
            scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
            scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
            geom_point(size = 4) +
            geom_vline(xintercept = 0, linetype = 4, color = "grey") +
            geom_hline(yintercept = 0, linetype = 4, color = "grey") +
            coord_equal() + 
            stat_ellipse(aes(fill = Group, group = Group),geom = "polygon", level = 0.95, alpha = 0.5) +
            geom_text_repel(aes(x = PC1, y = PC2, label = rownames(pcaData)), size = 3) 
        #	ggsave("PCA_deseq2_replicate.pdf",height=4*0.05*percentVar[2],width=4*0.05*percentVar[1])
        ggsave("PCA_2.pdf", height = 8, width = 8)
        ggsave("PCA_2.png",type="cairo-png",height=8,width=8)


        #PCA 3D
        #vsd <- varianceStabilizingTransformation(dds, blind = T)
        rv <- rowVars(assay(rld))
        select <- order(rv, decreasing = T)[seq_len(min(500, length(rv)))]
        pca <- prcomp(t(assay(rld)[select,]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        max_x<-round(percentVar[1]*100,2)
        max_y<-round(percentVar[2]*100,2)
        max_z<-round(percentVar[3]*100,2)
        x <- pca$x[,1]
        y <- pca$x[,2]
        z <- pca$x[,3]
        print(pca$x)
        print(percentVar)
        if (max_x > 1e-2 & max_y > 1e-2 & max_z > 1e-2)
        {
            dat3d <- data.frame(x, y, z, groups)
            print(groups)
            all_col <- factor(dat3d$groups, level = levels(dat3d$groups), labels = oe_col_qua(1 : length(seq_along(levels(dat3d$groups)))))
            pdf("PCA_3D.pdf", height = 10, width = 12)
            layout(matrix(1 : 2, 1, 2), width = c(7, 2))
            par(mar = c(5, 2, 5, 2))
            plot3d <- with(dat3d, scatterplot3d(x, y, z, xlab = paste("PC1:", max_x, "% variance"), ylab = paste("PC2:", max_y, "% variance") , zlab = paste("PC3:", max_z, "% variance"), main = "PCA 3D figure", color = as.character(all_col), pch = 16, label.tick.marks = TRUE, type = 'h', angle = 45))
            Position <- plot3d$xyz.convert(x, y, z)
            pointLabel(Position$x, Position$y, labels = gsub("Sample_", "", rownames(pca$x)), cex = 0.8, col = as.character(all_col))
            par(mar = c(0.5, 0.5, 0.5, 0.3))
            plot.new()
            legend("center", plot3d$xyz.convert(2, 0.5, 2), pch = 16, , xjust = - 1, yjust = - 1, legend = levels(dat3d$groups), col = oe_col_qua(1 : length(seq_along(levels(dat3d$groups)))))
            dev.off()
			
            png("PCA_3D.png",height=2400,width=2800,res=300)
            layout(matrix(1 : 2, 1, 2), width = c(7, 2))
            par(mar = c(5, 2, 5, 2))
            plot3d <- with(dat3d, scatterplot3d(x, y, z, xlab = paste("PC1:", max_x, "% variance"), ylab = paste("PC2:", max_y, "% variance") , zlab = paste("PC3:", max_z, "% variance"), main = "PCA 3D figure", color = as.character(all_col), pch = 16, label.tick.marks = TRUE, type = 'h', angle = 45))
            Position <- plot3d$xyz.convert(x, y, z)
            pointLabel(Position$x, Position$y, labels = gsub("Sample_", "", rownames(pca$x)), cex = 0.8, col = as.character(all_col))
            par(mar = c(0.5, 0.5, 0.5, 0.3))
            plot.new()
            legend("center", plot3d$xyz.convert(2, 0.5, 2), pch = 16, , xjust = - 1, yjust = - 1, legend = levels(dat3d$groups), col = oe_col_qua(1 : length(seq_along(levels(dat3d$groups)))))
            dev.off()

        }
        else
        {
            print(paste0(max_x, ",", max_y, ",", max_z))
            print("max_xlim or max_ylim or max_zlim is too small,PCA_3D is skipped!!")
        }
    }
}


###PCA.pdf  for no replicate samples
norep_pca_plot <- function(counts,list) {
    counts <- round(read.delim(counts, row.names = 1, check.names = F))
	if (! is.null(list)) 
    { list <- read.table(list, row.names = 1, header = T, sep = "\t", check.names = F)
	###################################################################################################
    counts_tmp<-subset(counts, select = c(rownames(list)))
	counts<-counts_tmp }
    samples <- colnames(counts)
    if (length(samples) < 3)
    {cat("Samples number is less than 3,stop!!", call. = FALSE)}
    else
    {
        colData <- data.frame(row.names = samples, Sample = samples)
        dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ Sample)
        ###PCA 2D
        rld <- rlog(dds, blind = T)
        pcaData <- plotPCA(rld, intgroup = c("Sample"), returnData = TRUE)
        percentVar <-round(100 * attr(pcaData, "percentVar"),2)	
        min <- round(min(pcaData$PC1, pcaData$PC2))
        max <- ceiling(max(pcaData$PC1, pcaData$PC2))

        ggplot(pcaData, aes(PC1, PC2, color = Sample)) +
            geom_point() +
            xlab(paste0("PC1: ", round(percentVar[1],2), "% variance")) +
            ylab(paste0("PC2: ", round(percentVar[2],2), "% variance")) +
            guides(col = guide_legend(nrow = 15)) +
            scale_color_manual(values = oe_col_qua(1 : length(unique(pcaData$Sample))),
            labels = c(as.character(unique(factor(pcaData$Sample)))),
            name = "Sample") +
            scale_fill_manual(values = oe_col_qua(1 : length(unique(pcaData$Sample))),
            labels = c(as.character(unique(factor(pcaData$Sample)))),
            name = "Sample") +
            theme_bw() +
            theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
            theme(text = element_text(size = 12, family = "ArialMT")) +
            scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
            scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
            geom_point(size = 4) +
            geom_vline(xintercept = 0, linetype = 4, color = "grey") +
            geom_hline(yintercept = 0, linetype = 4, color = "grey") +
            coord_equal()
        #	ggsave("PCA_deseq2_replicate.pdf",height=4*0.05*percentVar[2],width=4*0.05*percentVar[1])
        ggsave("PCA_1.pdf", height = 8, width = 8)
        ggsave("PCA_1.png",type="cairo-png",height=8,width=8)
        ggplot(pcaData, aes(PC1, PC2, color = Sample)) +
            geom_point() +
            xlab(paste0("PC1: ", round(percentVar[1],2), "% variance")) +
            ylab(paste0("PC2: ", round(percentVar[2],2), "% variance")) +
            guides(col = guide_legend(nrow = 15)) +
            scale_color_manual(values = oe_col_qua(1 : length(unique(pcaData$Sample))),
            labels = c(as.character(unique(factor(pcaData$Sample)))),
            name = "Sample") +
            scale_fill_manual(values = oe_col_qua(1 : length(unique(pcaData$Sample))),
            labels = c(as.character(unique(factor(pcaData$Sample)))),
            name = "Sample") +
            theme_bw() +
            theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black", size = 1, fill = NA)) +
            theme(text = element_text(size = 12, family = "ArialMT")) +
            scale_y_continuous(limits = c(min * 1.3, max * 1.3)) +
            scale_x_continuous(limits = c(min * 1.3, max * 1.3)) +
            geom_point(size = 4) +
            #geom_point(size = 6, alpha = 0.6) +
            geom_vline(xintercept = 0, linetype = 4, color = "grey") +
            geom_hline(yintercept = 0, linetype = 4, color = "grey") +
            coord_equal() +
            geom_text_repel(aes(x = PC1, y = PC2, label = rownames(pcaData)), size = 3) 
        #	ggsave("PCA_deseq2_replicate.pdf",height=4*0.05*percentVar[2],width=4*0.05*percentVar[1])
        ggsave("PCA_2.pdf", height = 8, width = 8)
        ggsave("PCA_2.png",type="cairo-png",height=8,width=8)

        ###PCA 3D
        groups <- as.character(colnames(counts))
        #vsd <- varianceStabilizingTransformation(dds, blind = T)
        rv <- rowVars(assay(rld))
        select <- order(rv, decreasing = T)[seq_len(min(500, length(rv)))]
        pca <- prcomp(t(assay(rld)[select,]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        max_x<-round(percentVar[1]*100,2)
        max_y<-round(percentVar[2]*100,2)
        max_z<-round(percentVar[3]*100,2)
        x <- pca$x[,1]
        y <- pca$x[,2]
        z <- pca$x[,3]
        print(pca$x)
        print(percentVar)
        if (max_x > 1e-2 & max_y > 1e-2 & max_z > 1e-2)
        {
            dat3d <- data.frame(x, y, z, groups)
            print(length(seq_along(levels(dat3d$groups))))
            all_col <- factor(dat3d$groups, level = levels(dat3d$groups), labels = oe_col_qua(1 : length(seq_along(levels(dat3d$groups)))))
            pdf("PCA_3D.pdf", height = 10, width = 12)
            layout(matrix(1 : 2, 1, 2), width = c(7, 2))
            par(mar = c(5, 2, 5, 2))
            plot3d <- with(dat3d, scatterplot3d(x, y, z, xlab = paste("PC1:", max_x, "% variance"), ylab = paste("PC2:", max_y, "% variance"), zlab = paste("PC3:", max_z, "% variance"), main = "PCA 3D figure", color = as.character(all_col), pch = 16, label.tick.marks = TRUE, type = 'h', angle = 45))
            Position <- plot3d$xyz.convert(x, y, z)
            pointLabel(Position$x, Position$y, labels = gsub("Sample_", "", rownames(pca$x)), cex = 0.8, col = oe_col_qua(1 : length(as.character(all_col))))
            par(mar = c(0.5, 0.5, 0.5, 0.3))
            plot.new()
            legend("center", plot3d$xyz.convert(2, 0.5, 2), pch = 16, , xjust = - 1, yjust = - 1, legend = levels(dat3d$groups), col = oe_col_qua(1 : length(seq_along(levels(dat3d$groups)))))
            dev.off()

            png("PCA_3D.png",height=2400,width=2800,res=300)
            layout(matrix(1 : 2, 1, 2), width = c(7, 2))
            par(mar = c(5, 2, 5, 2))
            plot3d <- with(dat3d, scatterplot3d(x, y, z, xlab = paste("PC1:", max_x, "% variance"), ylab = paste("PC2:", max_y, "% variance"), zlab = paste("PC3:", max_z, "% variance"), main = "PCA 3D figure", color = as.character(all_col), pch = 16, label.tick.marks = TRUE, type = 'h', angle = 45))
            Position <- plot3d$xyz.convert(x, y, z)
            pointLabel(Position$x, Position$y, labels = gsub("Sample_", "", rownames(pca$x)), cex = 0.8, col = oe_col_qua(1 : length(as.character(all_col))))
            par(mar = c(0.5, 0.5, 0.5, 0.3))
            plot.new()
            legend("center", plot3d$xyz.convert(2, 0.5, 2), pch = 16, , xjust = - 1, yjust = - 1, legend = levels(dat3d$groups), col = oe_col_qua(1 : length(seq_along(levels(dat3d$groups)))))
            dev.off()

        }
        else
        {
            print(paste0(max_x, ",", max_y, ",", max_z))
            print("max_xlim or max_ylim or max_zlim is too small,PCA_3D is skipped!!")
        }
    }
}
###sample2sample.pdf
sample_2_sample <- function(counts,group,list) {
    counts <- round(read.delim(counts, row.names = 1, check.names = F,quote=""))
	###################################################################################################
	if (! is.null(group)) 
    { phenodata <- read.table(group, row.names = 1, header = T, sep = "\t", check.names = F)

    counts_tmp<-subset(counts, select = c(rownames(phenodata)))
	counts<-counts_tmp }
	###################################################################################################
	if (! is.null(list)) 
    { list <- read.table(list, row.names = 1, header = T, sep = "\t", check.names = F)
    counts_tmp<-subset(counts, select = c(rownames(list)))
	counts<-counts_tmp }
	###################################################################################################
  samples <- colnames(counts)
  if (length(samples)<3)
  {cat("Samples number is less than 3,stop!!", call.=FALSE)}
  else
  {
    colData <- data.frame(row.names=samples, condition=samples)
    dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
    rld <- rlog(dds, blind=T)
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix)<- rld$condition
    colnames(sampleDistMatrix)<- rld$condition
    #heatmap.colors <-oe_col_seq("light-blue", length(sampleDistMatrix))
    #heatmap.colors <-colorRampPalette(c('blue','white','red'))(50)
    #heatmap.colors<-colorRampPalette(rev(c("red","gray90","green")))(102)
	heatmap.colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #20180426 modified by yxf
    pdf("sample2sample_distances_heatmap.pdf", height=1+0.9*length(rld$condition), width=1+0.9*length(rld$condition), onefile=F)
    #par(mar=c(4, 4, 4, 4))
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=heatmap.colors,
             # cellwidth = 30, cellheight = 30,
			fontsize_row=10,
			fontsize_col=10,
			treeheight_row=45,
			treeheight_col=45,
             main = 'Sample-to-Sample distances' )
             #show_rownames=T,
             #main = 'Sample-to-Sample distances',fontsize=5，fontsize_row = 6)           
			 #cellwidth = 30, cellheight = 30,
             #fontsize=5, fontsize_row = 6, fontsize_col =6,

    dev.off()
    png("sample2sample_distances_heatmap.png",res=300,height=2200,width=2200)
    pheatmap(sampleDistMatrix,
	clustering_distance_rows=sampleDists,
	clustering_distance_cols=sampleDists,
	col=heatmap.colors,
	fontsize_row=10,
	fontsize_col=10,
	treeheight_row=45,
	treeheight_col=45,
	main = 'Sample-to-Sample distances')
    dev.off()

  }
}
####main funcfion
suppressPackageStartupMessages(library(rgl))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(grDevices))
suppressPackageStartupMessages(library(genefilter))
suppressPackageStartupMessages(library(pca3d))
suppressPackageStartupMessages(library(maptools))
suppressPackageStartupMessages(library(oebio))
if (! is.null(opt$list) & ! is.null(opt$group)) {
    print_help(opt_parser)
    stop("Filtering samples can be performed using any one of (-s), (-l) at a time.", call. = FALSE)
}

if (! is.null(opt$group)) {
rep_pca_plot(opt$readcounts, opt$group)
sample_2_sample(opt$readcounts,opt$group,opt$list)
}else {
norep_pca_plot(opt$readcounts,opt$list)
sample_2_sample(opt$readcounts,opt$group,opt$list)
}


