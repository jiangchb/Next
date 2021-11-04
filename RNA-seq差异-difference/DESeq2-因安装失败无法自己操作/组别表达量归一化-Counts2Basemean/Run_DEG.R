#!/usr/bin/env Rscript
#############################################
#Author: Xiufeng Yang 
#email:  yxfhenu@163.com
#Creat Time: 2017-11-19 09:45:00
#############################################
#2021.0603
#changed by lipeng 
#p-vlalue and q-value
##################################################################################
## Method0:DEseq2 function define
run_DESeq2_nopaired <- function(counts, group, FC, Pvalue, FDR, Type,thr){
    control_s <- unlist(strsplit(as.character(group$control), ','))
    case_s <- unlist(strsplit(as.character(group$case), ','))
    control_name <- c(as.character(group$control_name))
    case_name <- c(as.character(group$case_name))
    samples <- c(control_s, case_s)
    name <- c(paste(case_name, control_name, sep = "-vs-"))
    colData <- data.frame(row.names = c(control_s, case_s),
    condition = factor(rep(c(control_name, case_name), c(length(control_s), length(case_s))),
    levels = c(control_name, case_name)))
    counts_tmp <- counts[1]
    for (n in 1 : length(samples)) {
        index <- which(names(counts[]) == samples[n])
        counts_tmp <- cbind(counts_tmp, counts[index])
    }
    rownames(counts_tmp) = counts_tmp[,1]
    #过滤全为零行信息
    counts_tmp = counts_tmp[,-1]
    # 若一个样本则不用取均，若多个样本则取均(一个样本组为vector不能直接apply 20210722)
    if (is.data.frame(counts_tmp[,case_s])){case_mean <- apply(counts_tmp[,case_s], 1, mean)}
    else{case_mean <- counts_tmp[,case_s,drop = FALSE]}
    if (is.data.frame(counts_tmp[,control_s])){con_mean <- apply(counts_tmp[,control_s], 1, mean)}
    else{con_mean <- counts_tmp[,control_s, drop = FALSE]}
    tmp_mean <- cbind(case_mean,con_mean)
    filter_mean <- subset(tmp_mean,subset=case_mean>thr|con_mean>thr)
    counts_tmp1<-counts_tmp[rownames(filter_mean),]
    counts_tmp <- counts_tmp1 %>% rownames_to_column(var="id")
    #需要给行名加表头
    write.table(counts_tmp, paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""), quote = F, row.names = F, sep = "\t")
    count.all <- read.table(paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""), header = T, sep = "\t", row.names = 1, check.names = F, quote = "", comment.char='')
    countData <- apply(count.all, 2, round)
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~  condition)
    dds <- DESeq(dds)
    res <- results(dds, independentFiltering = TRUE)
    res <- as.data.frame(res[order(rownames(res)),])
    res[, "FoldChange"] <- 2 ^ res[, "log2FoldChange"]
    res <- res[, c("baseMean", "lfcSE", "stat", "FoldChange", "log2FoldChange", "pvalue", "padj")]
    if (! is.null(Pvalue)) {
        resSig <- res[which(res$pvalue < Pvalue & abs(res$log2FoldChange) > log2(FC)),]
    }
    if (! is.null(FDR)) {
        resSig <- res[which(res$padj < FDR & abs(res$log2FoldChange) > log2(FC)),]
    }
    resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
    resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
    res <- cbind(rownames(res), res)
    colnames(res) <- c(paste0(Type, "_id"), "BaseMean", "lfcSE", "stat", "FoldChange", "log2FoldChange", "p-value", "q-value")
    write.table(res, paste0(name, "-all.", Type, ".xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
    resSig = cbind(rownames(resSig), resSig)
    colnames(resSig) <- c(paste0(Type, "_id"), "BaseMean", "lfcSE", "stat", "FoldChange", "log2FoldChange", "p-value", "q-value", "Regulation")
    if (! is.null(Pvalue)) {
        write.table(resSig, paste0(name, "-diff-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".xls"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
        stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff", paste("Total_diff(", "p-value<", Pvalue, "&|log2FC|>", log2(FC), ")", sep = ""), "Case_samples", "Control_samples"), ncol = 7)
    }
    if (! is.null(FDR)) {
        write.table(resSig, paste0(name, "-diff-", "q-val-", FDR, "-FC-", FC, ".", opt$type, ".xls"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
        stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff", paste("Total_diff(", "p-value<", Pvalue, "&|log2FC|>", log2(FC), ")", sep = ""), "Case_samples", "Control_samples"), ncol = 7)
    }
    up_num <- length(which(resSig$Regulation == "Up"))
    down_num <- length(which(resSig$Regulation == "Down"))
    total <- sum(up_num, down_num)
    stat1 <- c(case_name, control_name, up_num, down_num, total,group$case,group$control)
    if (file.exists(paste0(Type, "_diff_stat.xls"))) {
        stat <- rbind(stat1)}else {stat <- rbind(stat, stat1)
    }
    write.table(stat, paste0(Type, "_diff_stat.xls"), append = file.exists(paste0(Type, "_diff_stat.xls")), quote = F, row.names = F, col.names = F, sep = "\t", na = "")
    file.remove(paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""))
}
##################################################################################
## Method1:DEseq function define
run_DESeq <- function(counts, group, FC, Pvalue, FDR, Type,thr)
{
    case_s <- unlist(strsplit(as.character(group$case), split = ','))
    control_s <- unlist(strsplit(as.character(group$control), split = ','))
    case_name <- as.character(group$case_name)
    control_name <- as.character(group$control_name)
    c_names <- vector(length = length(c(control_s, case_s)))
    samples <- c(control_s, case_s)
    counts_tmp <- counts[1]
    for (n in 1 : length(samples)) {
        index <- which(names(counts[]) == samples[n])
        counts_tmp <- cbind(counts_tmp, counts[index])
    }
    
    rownames(counts_tmp) = counts_tmp[,1]
    #过滤全小于2信息
    counts_tmp = counts_tmp[,-1]
    # 若一个样本则不用取均，若多个样本则取均(一个样本组为vector不能直接apply 20210722)
    if (is.data.frame(counts_tmp[,case_s])){case_mean <- apply(counts_tmp[,case_s], 1, mean)}
    else{case_mean <- counts_tmp[,case_s,drop = FALSE]}
    if (is.data.frame(counts_tmp[,control_s])){con_mean <- apply(counts_tmp[,control_s], 1, mean)}
    else{con_mean <- counts_tmp[,control_s, drop = FALSE]}
    tmp_mean <- cbind(case_mean,con_mean)
    filter_mean <- subset(tmp_mean,subset=case_mean>thr|con_mean>thr)
    counts_tmp1<-counts_tmp[rownames(filter_mean),]
    counts_tmp <- counts_tmp1 %>% rownames_to_column(var="id")
    
    write.table(counts_tmp, paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""), quote = F, row.names = F, sep = "\t")
    count.all <- read.table(paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""), header = T, sep = "\t", row.names = 1, check.names = F, quote = "", comment.char='')
    data.all.int <- apply(count.all, 2, round)

    for (k in 1 : length(control_s)) {c_names[k] <- as.character(group$control_name)}
    for (k in length(control_s) + 1 : length(case_s)) {c_names[k] <- as.character(group$case_name)}

    Design <- data.frame(row.names = c(control_s, case_s), libType = c_names)
    case <- which(Design$libType == case_name)
    control <- which(Design$libType == control_name)
    countTable <- data.all.int[, c(case, control)]
    libType <- Design$libType[c(case, control)]
    cds <- newCountDataSet(countTable, libType)
    cds = estimateSizeFactors(cds)
    if (group[1, 5] == "yes")
    {
        cdsnew <- tryCatch(estimateDispersions(cds), error = function(e) NULL)
        if (is.null(cdsnew))
        {
            cdsnew <- tryCatch(estimateDispersions(cds, method = c("pooled"), fitType = c("local")), error = function(e) NULL)
            if (is.null(cdsnew))
            {
                cat(paste(case_name, "-vs-", control_name, ': Both  Default parameter and  pooled ,local parameter is not suit to the Replicate Samples,Try blind ,local ,fit-only parameter', sep = ""), call. = FALSE)
                cdsnew <- estimateDispersions(cds, method = c("blind"), sharingMode = c("fit-only"), fitType = c("local"))
            }
            else
            {
                cat(paste(case_name, "-vs-", control_name, ':Using pooled ,local parameter for Replicate Samples\n', sep = ""), file = stderr())
            }
        }
        else
        {
            cat(paste(case_name, "-vs-", control_name, ':Using Default parameter for Replicate Samples\n', sep = ""), file = stderr())
        }
    }
    if (group[1, 5] == "no") {
        cdsnew <- estimateDispersions(cds, method = c("blind"), sharingMode = c("fit-only"), fitType = c("local"))
        cat(paste(case_name, "-vs-", control_name, ':Using blind ,local ,fit-only parameter  for No Replicate Samples\n', sep = ""), file = stderr())
    }
    res <- nbinomTest(cdsnew, control_name, case_name)
    colnames(res)[1] <- paste0(opt$type, "_id")
    colnames(res)[3] <- paste("BaseMean_control_", control_name, sep = "")
    colnames(res)[4] <- paste("BaseMean_case_", case_name, sep = "")
    colnames(res)[5] <- c("FoldChange")
    colnames(res)[7] <- c("pValue")
    colnames(res)[8] <- c("qValue")
    sizefactor <- as.matrix(sizeFactors(cdsnew))
    res_new <- res
####重命名输出列表的表头
    colnames(res_new)[2] <- c("BaseMean")
    colnames(res_new)[7] <- c("p-value")
    colnames(res_new)[8] <- c("q-value")
    write.table(res_new, paste(case_name, "-vs-", control_name, "-all.", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")

    zero <- which(res[, 2] == "0" &
        res[, 3] == "0" &
        res[, 4] == "0")
    if (length(zero) >= 1) { res2 <- res[- zero,]}
    else {res2 <- res}
    int <- which(round(res2[, 2]) == res2[, 2])
    if (length(int) >= 1) {
        res2[int, 2] <- as.numeric(paste(res2[int, 2], ".000000001", sep = ""))
    }


    if (! is.null(Pvalue)) {
        #pdf
        pdf(paste(case_name, "-vs-", control_name, "-MA-p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".pdf", sep = ""))
        plotMA(res2, col = ifelse((as.numeric(res2$pValue) < Pvalue) & (as.numeric(res2$FoldChange) < (1 / FC) | as.numeric(res2$FoldChange) > FC), oe_col_qua(2), oe_col_qua(8)),
        main = c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), ":p-value", "<", Pvalue, "&& |log2FC|>", round(log2(FC), 2), sep = "")))
        dev.off()
        #png
        png(paste(case_name, "-vs-", control_name, "-MA-p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".png", sep = ""), height = 2400, width = 2400, res = 300)
        plotMA(res2, col = ifelse((as.numeric(res2$pValue) < Pvalue) & (as.numeric(res2$FoldChange) < (1 / FC) | as.numeric(res2$FoldChange) > FC), oe_col_qua(2), oe_col_qua(8)),
        main = c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), ":p-value", "<", Pvalue, "&& |log2FC|>", round(log2(FC), 2), sep = "")))
        dev.off()
        #diff get
        diff.pos <- which(as.numeric(res$pValue) < Pvalue & (as.numeric(res$FoldChange) < (1 / FC) | as.numeric(res$FoldChange) > FC))
        if (length(diff.pos) == 0) {
            cat(paste("No gene/transcript filtered by for ", as.character(group$case_name), "-vs-", as.character(group$control_name), ":p-value", "<", Pvalue, "&&|log2FC|>", log2(FC), sep = ""), call. = FALSE)
        }
    }

    if (! is.null(FDR)) {
        #pdf
        pdf(paste(case_name, "-vs-", control_name, "-MA-q-val-", FDR, "-FC-", FC, ".", opt$type, ".pdf", sep = ""))
        plotMA(res2, col = ifelse((as.numeric(res2$qValue) <= FDR) & (as.numeric(res2$FoldChange) < (1 / FC) | as.numeric(res2$FoldChange) > FC), oe_col_qua(2), oe_col_qua(8)),
        main = c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), ":q-value", "<", FDR, "&& |log2FC|>", round(log2(FC), 2), sep = "")))
        dev.off()
        #png
        png(paste(case_name, "-vs-", control_name, "-MA-q-val-", FDR, "-FC-", FC, ".", opt$type, ".png", sep = ""), height = 2400, width = 2400, res = 300)
        plotMA(res2, col = ifelse((as.numeric(res2$qValue) <= FDR) & (as.numeric(res2$FoldChange) < (1 / FC) | as.numeric(res2$FoldChange) > FC), oe_col_qua(2), oe_col_qua(8)),
        main = c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), ":q-value", "<", FDR, "&& |log2FC|>", round(log2(FC), 2), sep = "")))
        dev.off()
        diff.pos <- which(as.numeric(res$qValue) <= FDR & (as.numeric(res$FoldChange) < (1 / FC) | as.numeric(res$FoldChange) > FC))
        if (length(diff.pos) == 0) {
            cat(paste("No gene/transcript filtered by for ", as.character(group$case_name), "-vs-", as.character(group$control_name), ":q-value", "<", FDR, "&&|log2FC|>", round(log2(FC), 2), sep = ""), call. = FALSE)
        }
    }

    res.diff <- res[diff.pos,]
    up.pos <- which(res.diff[, 3] > res.diff[, 4])
    down.pos <- which(res.diff[, 3] < res.diff[, 4])
    up_down <- matrix(NA, nrow = dim(res.diff)[1])
    up_down[up.pos, 1] <- "Down"
    up_down[down.pos, 1] <- "Up"
    res.diff.up_down <- cbind(res.diff, up_down)
    colnames(res.diff.up_down)[2] <- c("BaseMean")
    colnames(res.diff.up_down)[9] <- c("Regulation")
    colnames(res.diff.up_down)[7] <- c("p-value")
    colnames(res.diff.up_down)[8] <- c("q-value")
    if (! is.null(Pvalue))
    {
        write.table(res.diff.up_down, paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")
        stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff", paste("Total_diff(", "p-value<", Pvalue, "&|log2FC|>", round(log2(FC), 2), ")", sep = ""), "Case_samples", "Control_samples"), ncol = 7)
    }
    if (! is.null(FDR))
    {
        write.table(res.diff.up_down, paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "q-val-", FDR, "-FC-", FC, ".", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")
        stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff", paste("Total_diff(", "q-value<", FDR, "&|log2FC|>", round(log2(FC), 2), ")", sep = ""), "Case_samples", "Control_samples"), ncol = 7)
    }
    up_num <- length(which(up_down[, 1] == "Up"))
    down_num <- length(which(up_down[, 1] == "Down"))
    total <- sum(up_num, down_num)
    if (length(control_s) > 1 | length(case_s) > 1) {
        stat1 <- c(case_name, control_name, up_num, down_num, total,group$case, group$control)
    }
    else {
        stat1 <- c(case_name, control_name, up_num, down_num, total,group$case, group$control)
    }

    if (file.exists(paste0(Type, "_diff_stat.xls"))) {stat <- rbind(stat1)}else {stat <- rbind(stat, stat1)}
    write.table(stat, paste0(Type, "_diff_stat.xls"), append = file.exists(paste0(Type, "_diff_stat.xls")), quote = F, row.names = F, col.names = F, sep = "\t", na = "")
    file.remove(paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""))
}

##################################################################################
##Method2 DEseq2_paired#下方均已修改20210603
run_DESeq2_paired <- function(counts, group, FC, Pvalue, FDR, Type,thr){
    control_s <- unlist(strsplit(as.character(group$control), ','))
    case_s <- unlist(strsplit(as.character(group$case), ','))
    control_name <- c(as.character(group$control_name))
    case_name <- c(as.character(group$case_name))
    samples <- c(control_s, case_s)
    name <- c(paste(case_name, control_name, sep = "-vs-"))
    colData <- data.frame(row.names = c(control_s, case_s),
    type = factor(seq(1 : length(control_s)), seq(1 : length(case_s))),
    condition = factor(rep(c(control_name, case_name), c(length(control_s), length(case_s))),
    levels = c(control_name, case_name)))
    counts_tmp <- counts[1]
    for (n in 1 : length(samples)) {
        index <- which(names(counts[]) == samples[n])
        counts_tmp <- cbind(counts_tmp, counts[index])
    }
    
    #过滤全为零行信息
    counts_tmp = counts_tmp[,-1]
    # 若一个样本则不用取均，若多个样本则取均(一个样本组为vector不能直接apply 20210722)
    if (is.data.frame(counts_tmp[,case_s])){case_mean <- apply(counts_tmp[,case_s], 1, mean)}
    else{case_mean <- counts_tmp[,case_s,drop = FALSE]}
    if (is.data.frame(counts_tmp[,control_s])){con_mean <- apply(counts_tmp[,control_s], 1, mean)}
    else{con_mean <- counts_tmp[,control_s, drop = FALSE]}
    tmp_mean <- cbind(case_mean,con_mean)
    filter_mean <- subset(tmp_mean,subset=case_mean>thr|con_mean>thr)
    counts_tmp1<-counts_tmp[rownames(filter_mean),]
    counts_tmp <- counts_tmp1 %>% rownames_to_column(var="id")
    
    write.table(counts_tmp, paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""), quote = F, row.names = F, sep = "\t")
    count.all <- read.table(paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""), header = T, sep = "\t", row.names = 1, check.names = F, quote = "", comment.char='')
    countData <- apply(count.all, 2, round)
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ type + condition)
    dds <- DESeq(dds)
    res <- results(dds, independentFiltering = TRUE)
    #    png(paste0(name,"-MA-", opt$type,".png"), height=800, width=800, res=100)
    #    plotMA(res, main="DESeq2")
    #    dev.off()
    #    pdf(paste0(name,"-MA-", opt$type,".pdf"))
    #    plotMA(res, main="DESeq2")
    #    dev.off()
    res <- as.data.frame(res[order(rownames(res)),])
    res[, "FoldChange"] <- 2 ^ res[, "log2FoldChange"]
    res <- res[, c("baseMean", "lfcSE", "stat", "FoldChange", "log2FoldChange", "pvalue", "padj")]
    if (! is.null(Pvalue)) {
        resSig <- res[which(res$pvalue < Pvalue & abs(res$log2FoldChange) > log2(FC)),]
    }
    if (! is.null(FDR)) {
        resSig <- res[which(res$padj < FDR & abs(res$log2FoldChange) > log2(FC)),]
    }
    resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
    resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
    res <- cbind(rownames(res), res)
    colnames(res) <- c(paste0(Type, "_id"), "BaseMean", "lfcSE", "stat", "FoldChange", "log2FoldChange", "p-value", "q-value")
    write.table(res, paste0(name, "-all.", Type, ".xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
    resSig = cbind(rownames(resSig), resSig)
    colnames(resSig) <- c(paste0(Type, "_id"), "BaseMean", "lfcSE", "stat", "FoldChange", "log2FoldChange", "p-value", "q-value", "Regulation")
    if (! is.null(Pvalue)) {
        write.table(resSig, paste0(name, "-diff-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".xls"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
        stat <- matrix(c("Case", "Control", "Up_diff", "Down_diff", paste("Total_diff(", "p-value<", Pvalue, "&|log2FC|>", log2(FC), ")", sep = ""), "Case_samples", "Control_samples"), ncol = 7)
    }
    if (! is.null(FDR)) {
        write.table(resSig, paste0(name, "-diff-", "q-val-", FDR, "-FC-", FC, ".", opt$type, ".xls"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
        stat <- matrix(c("Case","Control", "Up_diff", "Down_diff", paste("Total_diff(", "q-value<", FDR, "&|log2FC|>", log2(FC), ")", sep = ""),"Case_samples", "Control_samples"), ncol = 7)
    }
    up_num <- length(which(resSig$Regulation == "Up"))
    down_num <- length(which(resSig$Regulation == "Down"))
    total <- sum(up_num, down_num)
    stat1 <- c(case_name,control_name, up_num, down_num, total,group$case, group$control)
    if (file.exists(paste0(Type, "_diff_stat.xls"))) {
        stat <- rbind(stat1)}else {stat <- rbind(stat, stat1)
    }
    write.table(stat, paste0(Type, "_diff_stat.xls"), append = file.exists(paste0(Type, "_diff_stat.xls")), quote = F, row.names = F, col.names = F, sep = "\t", na = "")
    file.remove(paste(case_name, "-vs-", control_name, "-counts_tmp.xls", sep = ""))}
##################################################################################
##diff_stat_barplot
diff_stat_barplot <- function(diff_stat, type){
  library(ggplot2)
  library(stringr)
  data <- read.table(diff_stat, header = T, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = F);
#  for (i in 1 : length(data$Control))
#  {
#    if (! is.na(data$Case[i]))
#    {
#      data$Control[i] <- data$Control[i]
#    }
#    if (! is.na(data$Case[i])) {
#      data$Case[i] <- data$Case[i]
#    }
#  }
  stat <- matrix(nrow = 2 * length(data$Control), ncol = 3)
  for (i in 1 : length(data$Control)) {
    stat[(2 * i - 1),] <- cbind(paste0(data$Case[i], "-vs-", data$Control[i]), "Up", data$Up_diff[i])
    stat[(2 * i),] <- cbind(paste0(data$Case[i], "-vs-", data$Control[i]), "Down", data$Down_diff[i])
  }
  data2 <- as.data.frame(stat)
  colnames(data2) <- c("Group", "Type", "Gene_number")
  write.table(data2, paste0(type, "_diff_stat_barplot.xls"), quote = F, row.names = F, col.names = T, sep = "\t", na = "")
  data2 <- read.table(paste0(type, "_diff_stat_barplot.xls"), header = T, sep = "\t", quote = "", check.names = F);
  uniq <- data2$Group[! duplicated(data2$Group)]
  data2$Group <- factor(data2$Group, levels = uniq)
  data2$Type <- factor(data2$Type, levels = c("Up", "Down"));
  p = ggplot(data = data2, aes(x = Group, y = Gene_number, fill = Type,width=0.1*length(colnames(data2)))) +
    geom_bar(stat = "identity", position = "dodge" ) +
    scale_color_manual(values = rev(oe_col_qua(1 : length(unique(data2$Type)))),
                       labels = c(as.character(unique(factor(data2$Type)))),
                       name = "Type") +
    scale_fill_manual(values =rev(oe_col_qua(1: length(unique(data2$Type)))),
                      labels = c(as.character(unique(factor(data2$Type)))),
                      name = "Type") +
    geom_text(aes(label = data2$Gene_number), size = 3, colour = 'black', vjust = -0.6, hjust = 0.5, position = position_dodge(0.4)) +
    xlab("") +
    ylab(paste0("Differently Expressed ", type, " number")) +
    labs(title = paste0("Statistic of Differently Expressed ", type)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),plot.margin = margin(0.5, 2, 0.1, 3, "cm")) +
    theme(panel.border = element_rect(fill = NA, colour = "black"),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(axis.text.y = element_text(size = 9, color = "black")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,  vjust=1,size = 9, color = "black")) +
    theme(legend.text = element_text(size = 9)) +
    theme(panel.grid = element_blank())
  ggsave(paste0(type, "_diff_stat_barplot.pdf"), height = 4+0.5*length(colnames(data2)), width = 4+0.5*length(colnames(data2)), plot = p)
  ggsave(paste0(type, "_diff_stat_barplot.png"), type = "cairo-png", height = 6, width = 8+0.2*length(colnames(data2)), plot = p)
  file.remove(paste0(type, "_diff_stat_barplot.xls"))
}

##################################################################################
## annotation  function define
add_expression_annotation <- function(group, FC, Pvalue, FDR, expression, annotation){
    ex_name <- sub(".txt", "", sub(".xls", "", opt$expression))
    count_name <- sub(".txt", "", sub(".xls", "", opt$readcounts))
    if (! file.exists(paste(ex_name, "_anno.xls", sep = "")) | ! file.exists(paste(count_name, "_anno.xls", sep = ""))) {
        express_anno <- merge(expression, annotation, by = 1, all.x = T)
        colnames(express_anno)[1] <- paste(opt$type, "_id", sep = "")
        write.table(express_anno, paste(ex_name, "_anno.xls", sep = ""), row.names = F, quote = F, na = "", sep = "\t")
        counts_anno <- merge(counts, annotation, by = 1, all.x = T)
        colnames(counts_anno)[1] <- paste(opt$type, "_id", sep = "")
        write.table(counts_anno, paste(count_name, "_anno.xls", sep = ""), row.names = F, quote = F, na = "", sep = "\t")
    }

    all <- read.table(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-all", ".", opt$type, ".xls", sep = ""), header = T, sep = "\t", check.names = F, quote = "", comment.char='')
    if (! is.null(Pvalue)) {
        diff <- read.table(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".xls", sep = ""), header = T, sep = "\t", check.names = F, comment.char='')}
    if (! is.null(FDR)) {
        diff <- read.table(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "q-val-", FDR, "-FC-", FC, ".", opt$type, ".xls", sep = ""), header = T, sep = "\t", check.names = F, comment.char='')}

    case_s_names <- unlist(strsplit(as.character(group$case), split = ','))
    control_s_names <- unlist(strsplit(as.character(group$control), split = ','))
    samples <- c(control_s_names, case_s_names)
    expression.tmp <- expression[1]
    for (n in 1 : length(samples))
    {
        index <- which(names(expression[]) == samples[n])
        expression.tmp <- cbind(expression.tmp, expression[index])
    }
    all.anno <- annotation[match(all[, 1], annotation[, 1]),]
    all.expression <- expression.tmp[match(all[, 1], expression.tmp[, 1]),]
    colnames(all.expression) <- paste("Expression", colnames(all.expression), sep = "_")
    all <- cbind(all, all.expression[, 2 : ncol(all.expression)], all.anno[, 2 : ncol(all.anno)])
    write.table(all, paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-all", ".", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")

    diff.anno <- annotation[match(diff[, 1], annotation[, 1]),]
    diff.expression <- expression.tmp[match(diff[, 1], expression.tmp[, 1]),]
    colnames(diff.expression) <- paste("Expression", colnames(diff.expression), sep = "_")
    diff <- cbind(diff, diff.expression[, 2 : ncol(diff.expression)], diff.anno[, 2 : ncol(diff.anno)])
    if (! is.null(Pvalue))
    {
        write.table(diff, paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")
    }
    if (! is.null(FDR))
    {
        write.table(diff, paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "q-val-", FDR, "-FC-", FC, ".", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")
    }
}
##################################################################################
## heatmap_plot function define
heatmap_plot <- function(expression, group, FC, Pvalue, FDR){
    colors <- oe_col_div("rwb", 256)
    case_s_names <- unlist(strsplit(as.character(group$case), split = ','))
    control_s_names <- unlist(strsplit(as.character(group$control), split = ','))
    samples <- c(as.character(control_s_names), as.character(case_s_names))
    if (! is.null(Pvalue)) {
        f <- c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type, ".xls", sep = ""))
    }
    if (! is.null(FDR)) {
        f <- c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-diff-", "q-val-", FDR, "-FC-", FC, ".", opt$type, ".xls", sep = ""))
    }
    f_rn <- rownames(read.delim(f, header = T, sep = "\t", row.names = 1, comment.char=''))
    diff.expression <- expression[match(f_rn, expression[, 1]), samples]
    if (length(diff.expression) < 2 || length(f_rn) < 2 ) {
        cat("less than 2 diff gene/transcript,heatmap cluster is skipped!!")
    }
    else {
        rownames(diff.expression) = gsub('>.*','',f_rn)
        ind <- apply(diff.expression, 1, mean) > 0
        diff.expression <- diff.expression[ind,]
        if (! is.null(Pvalue)) {
            picname = c(paste0(as.character(group$case_name), "-vs-", as.character(group$control_name), "-heatmap-", "p-val-", Pvalue, "-FC-", FC, ".", opt$type))
            type = c("p-value")
            thresh = Pvalue
        }
        if (! is.null(FDR)) {
            picname = c(paste0(as.character(group$case_name), "-vs-", as.character(group$control_name), "-heatmap-", "q-val-", FDR, "-FC-", FC, ".", opt$type))
            type <- c("q-value")
          thresh = FDR
        }
        case_rep = rep(as.character(group$case_name), length(case_s_names))
        control_rep = rep(as.character(group$control_name), length(control_s_names))
        anno = c(control_rep, case_rep)
        annotation_col = data.frame(condition = factor(anno))
        rownames(annotation_col) = samples
        breaksList = seq(0, 5, by = 0.1)
        data <- log2(diff.expression + 0.0001)
        if (length(f_rn) < 30) {
            xx<-pheatmap(data, annotation_col = annotation_col,
            scale = "row",
            main = c(paste0(as.character(group$case_name), "-vs-", as.character(group$control_name), ":", type, "<", thresh, "&& |log2FC|>", round(log2(FC), 2))),
            color = colors,
            show_rownames = T)
        }
        else {
            xx<-pheatmap(data, annotation_col = annotation_col,
            scale = "row",
            main = c(paste0(as.character(group$case_name), "-vs-", as.character(group$control_name), ":", type, "<", thresh, "&& |log2FC|>", round(log2(FC), 2))),
            color = colors,
            show_rownames = F)
        }
        save_pheatmap <- function(x, filename) {
        stopifnot(!missing(x))
        stopifnot(!missing(filename))
        #pdf
        pdf(paste0(filename, ".pdf"))
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
        #png
	    png(paste0(filename, ".png"), height=2500, width=2500, res=300)
	    grid::grid.newpage()
	    grid::grid.draw(x$gtable)
	    dev.off()
         }
         save_pheatmap(xx, picname)
         if(file.exists("Rplots.pdf")){ file.remove("Rplots.pdf") }
    }
}
##################################################################################
## volcano_plot function define
volcano_plot <- function(group, FC, Pvalue, FDR){
    fv <- function (d, p, n1, n2, alpha, f1, name, title, cex = 0.4, pch = 19,
    cols = c(oe_col_qua(8), oe_col_qua(3), oe_col_qua(2), oe_col_qua(8)),
    ltys = c(1, 3), use_legend = TRUE, ...)
    {
        f2 <- abs(d) < log2(FC)
        col <- rep(cols[1], length(d))
        col[! f2 & p < alpha & d > 0] <- cols[3]
        col[! f2 & p < alpha & d < 0] <- cols[2]
        col[! f2 & p >= alpha] <- cols[4]
        xlab = expression(paste(log[2], " FoldChange"))
        if (name == "p-value") {ylab = expression(paste("-", log[10], "p-value"))}
        if (name == "q-value") {ylab = expression(paste("-", log[10], "q-value"))}
        plot(d, - log10(p), cex = cex, pch = pch, xlab = xlab, ylab = ylab, col = col, main = title, ...)
        abline(h = - log10(alpha), v = - f1, col = oe_col_qua(8))
        abline(h = - log10(alpha), v = f1, col = oe_col_qua(8))
        legend("topleft", c("Filtered", "Up" , "Down"), pch = pch,
        col = c(oe_col_qua(8), oe_col_qua(2), oe_col_qua(3)), inset = 0.025, cex = 0.8, , bg = "white")
    }
    data <- read.delim(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-all.", opt$type, ".xls", sep = ""), header = T, row.names = 1, sep = "\t", quote = "", comment.char='')
    #picname <-c(paste("Volcano-","-vs-",".pdf",sep=""))
    if (! is.null(Pvalue)) {
        data <- data[which(data['p.value'] != ""),]
        P <- data['p.value']
        picname <- c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-volcano", "-p-val-", Pvalue, "-FC-", FC, ".", opt$type, sep = ""))
        type = c("p-value")
        thresh = Pvalue
    }
    if (! is.null(FDR)) {
        data <- data[which(data['q.value'] != ""),]
        P <- data['q.value']
        picname <- c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), "-volcano", "-q-val-", FDR, "-FC-", FC, ".", opt$type, sep = ""))
        type <- c("q-value")
        thresh = FDR
    }
    log2FC <- data$log2FoldChange
    df <- data.frame(P, log2FC)
    colnames(df)[1] = "P"
    tmp <- df[which(df$log2FC != "Inf" & df$log2FC != "-Inf"),]
    max = max(tmp$log2FC)
    min = min(tmp$log2FC)
    df$log2FC <- as.numeric(sub("-Inf", min, df$log2FC))
    df$log2FC <- as.numeric(sub("Inf", max, df$log2FC))
    ##pdf
    pdf(paste0(picname, ".pdf"), height = 8, width = 7)
    cat(picname)
    fv(unlist(df$log2FC), unlist(df$P), n1, n2, alpha = c(thresh), f1 = round(log2(FC), 2), name = c(type),
    title = c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), ":", type, "<", thresh, "&& |log2FC|>", round(log2(FC), 2))))
    dev.off()
    ##png
    png(paste0(picname, ".png"), height = 2500, width = 2200, res = 300)
    fv(unlist(df$log2FC), unlist(df$P), n1, n2, alpha = c(thresh), f1 = round(log2(FC), 2), name = c(type),
    title = c(paste(as.character(group$case_name), "-vs-", as.character(group$control_name), ":", type, "<", thresh, "&& |log2FC|>", round(log2(FC), 2))))
    dev.off()
}



########################################################################################################################
##                                             main function moduel                                               ######
########################################################################################################################
########################################################################################################################
#parameter  import
suppressPackageStartupMessages(library("optparse"))
option_list = list(
make_option(c("-e", "--expression"), type = "character", default = NULL,
help = "expression matrix file name", metavar = "character"),
make_option(c("-r", "--readcounts"), type = "character", default = NULL,
help = "read counts matrix file name", metavar = "character"),
make_option(c("-d", "--group"), type = "character", default = NULL,
help = "comparison information file name", metavar = "character"),
make_option(c("-a", "--annotation"), type = "character", default = NULL,
help = "annotation file name", metavar = "character"),
make_option(c("-p", "--pvalue"), type = "double",
help = "pvalue ratio threshold.Filtering can be performed using any one of (-p), (-f) at a time", metavar = "double"),
make_option(c("-f", "--fdr"), type = "double",
help = "fdr ratio threshold.Filtering can be performed using any one of (-p), (-f) at a time", metavar = "double"),
make_option(c("-c", "--foldchange"), type = "double", default = 2.0,
help = "foldchange threshold [default %default]", metavar = "double"),
make_option(c("-t", "--type"), type = "character", default = NULL,
help = "transcript or gene type:mRNA,circRNA,lncRNA,miRNA,gene", metavar = "character"),
make_option(c("-o", "--outputdir"), type = "character", default = "./",
help = "output directory for results", metavar = "character"));
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

########################################################################################################################
## get this R script file's own directory
args <- commandArgs(trailingOnly = F)
script.dir <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))


#######import library###################################################################################################
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(oebio))
########parameter check
if (is.null(opt$group) |
    is.null(opt$readcounts) |
    (is.null(opt$pvalue) & is.null(opt$fdr))) {
    print_help(opt_parser)
    stop("read counts matrix file, comparison information file must be supplied", call. = FALSE)
}
if (! is.null(opt$pvalue) & ! is.null(opt$fdr)) {
    print_help(opt_parser)
    stop("Filtering can be performed using any one of (-p), (-f) at a time.", call. = FALSE)
}
if (! is.null(opt$expression)) {
    express_file_abpath <- normalizePath(opt$expression)
    expression <- read.table(normalizePath(opt$expression), header = T, check.names = FALSE, quote = "", comment.char='')
    colnames(expression)[1] <- paste(opt$type, "_id")
}else {
    cat("expression matrix file is not supplied, heatmap for replicate groups  will be skipped\n")
}
if (! is.null(opt$annotation)) {
    annotation <- read.table(normalizePath(opt$annotation), sep = "\t", header = T, quote = "", na.strings = "", fill = TRUE, comment.char='')
}
if (is.null(opt$expression) | is.null(opt$annotation)) {
    cat("expression matrix file or annotation files is not supplied, annotaion for DEG  will be skipped\n")
}

########data import#####################################################################################################
counts <- read.table(normalizePath(opt$readcounts), header = T, sep = "\t", check.names = FALSE, quote = "", comment.char='')
group <- read.table(normalizePath(opt$group), header = T, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "", colClasses = "character")
DESeq_group <- group[which(group$method == "DESeq"),]
DESeq2_group <- group[which(group$method == "DESeq2" & group$paired == "no"),]
smallRNA_no_rep_group <- group[which(group$method == "Audic clavery"),]
paired_group <- group[which(group$paired == "yes" &
    group$replicate == "yes" &
    group$method == "DESeq2"),]
rep_group <- group[which(group$replicate == "yes"),]
express_file_abpath<-normalizePath(opt$expression)
counts_file_abpath <-normalizePath(opt$readcounts)
group_file_abspath <-normalizePath(opt$group)

#####################outputdir check####################################################################################

if (file.exists(opt$outputdir) == "FALSE") {
    dir.create(opt$outputdir)
}
setwd(opt$outputdir)
if (file.exists(paste0(opt$type, "_diff_stat.xls")) == TRUE) {
    suppressPackageStartupMessages(library(ggplot2))
    file.remove(paste0(opt$type, "_diff_stat.xls"))
}

####################Method1: smallRNA no replicate######################################################################
if (length(rownames(smallRNA_no_rep_group)) > 0) {
    print("DEG  for smallRNA noreplicate groups is beginning:")
    if(!is.null(opt$pvalue)){
        cmd<-paste("perl ",script.dir,"/smallRNA_norep_DEG.pl"," -e ",express_file_abpath," -r ",counts_file_abpath," -g ", group_file_abspath," -p ",opt$pvalue," -c " ,opt$foldchange, " -t ",opt$type,sep="")
    }
    if(!is.null(opt$fdr)){
        cmd<-paste("perl ",script.dir,"/smallRNA_norep_DEG.pl"," -e ",express_file_abpath," -r ",counts_file_abpath," -g ", group_file_abspath," -f ",opt$fdr, " -c " ,opt$foldchange," -t ",opt$type,sep="")
    }
    print(cmd)
    system(cmd)
    print("DEG  for smallRNA noreplicate groups is OK")
}#######################Method2: DEseq##################################################################################
if (length(rownames(DESeq_group)) > 0) {
    print("DEG  using DESeq is beginning:")
    registerDoParallel(cores = length(rownames(DESeq_group)))
    suppressPackageStartupMessages(library(DESeq))
    foreach(i = 1 : nrow(DESeq_group)) %dopar% run_DESeq(counts, DESeq_group[i,], opt$foldchange, opt$pvalue, opt$fdr, opt$type,2)
    print("DEG using DESeq is OK")
    doParallel::stopImplicitCluster()
}
###################Method3:DESeq2 nopaired Method################################################################################
if (length(rownames(DESeq2_group)) > 0) {
        print ("DESeq2_group")
    print("DEG  using DESeq2 no_paired is beginning:")
    registerDoParallel(cores = length(rownames(DESeq2_group)))
    suppressPackageStartupMessages(library(DESeq2))
        print ("library DESeq2")
    foreach(i = 1 : nrow(DESeq2_group)) %dopar% run_DESeq2_nopaired(counts, DESeq2_group[i,], opt$foldchange,opt$pvalue, opt$fdr, opt$type,2)
    print("DEG using DESeq2 no_paired is OK")
    doParallel::stopImplicitCluster()
}

###################Method3:paired Method################################################################################
if (length(rownames(paired_group)) > 0) {
    print("DEG  for paired groups is beginning:")
    print(paired_group)
    registerDoParallel(cores = length(rownames(paired_group)))
    suppressPackageStartupMessages(library(DESeq2))
    foreach(i = 1 : nrow(paired_group)) %dopar% run_DESeq2_paired(counts, paired_group[i,], opt$foldchange, opt$pvalue, opt$fdr, opt$type,2)
    doParallel::stopImplicitCluster()
    print("DEG  for paired groups is OK")
}
###################diff_stat_barplot####################################################################################
print("barplot  for DEG statistic results is beginning:")
diff_stat_barplot(paste0(opt$type, "_diff_stat.xls"), opt$type)
print("barplot for DEG statistic results is OK:")
###################heatmap#############################
if (! is.null(opt$expression)) {
    registerDoParallel(cores = length(rownames(rep_group)))
    if (length(rownames(rep_group)) > 0)
    {   print("heatmap is beginning:")
        foreach(i = 1 : nrow(rep_group)) %dopar% heatmap_plot(expression, rep_group[i,], opt$foldchange, opt$pvalue, opt$fdr)
        print("heatmap is OK")
        doParallel::stopImplicitCluster()
    }
}
###################volcano##############################################################################################
registerDoParallel(cores = length(rownames(group)))
print("volcano is beginning:")
foreach(i = 1 : nrow(group)) %dopar% volcano_plot(group[i,], opt$foldchange, opt$pvalue, opt$fdr)
print("volcano is OK")
doParallel::stopImplicitCluster()

###################annotation###########################################################################################
if (! is.null(opt$expression) & ! is.null(opt$annotation)) {
    registerDoParallel(cores = length(rownames(group)))
    print("annotation is  beginning:")
    foreach(i = 1 : nrow(group)) %dopar% add_expression_annotation(group[i,], opt$foldchange, opt$pvalue, opt$fdr, expression, annotation)
    print("annotation is OK")
    doParallel::stopImplicitCluster()
}
