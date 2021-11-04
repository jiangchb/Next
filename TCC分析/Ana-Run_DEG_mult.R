#!/usr/bin/env Rscript
#############################################
#Author: Xuan Tang 
#email:  tangxuan@oebiotech.com
#Creat Time: 2019-8-28 16:30:00
#############################################

##################################################################################
## Method:TCC function define
run_TCC <- function(counts, sample_group, FDR, Type){
    samples <- unlist(as.character(sample_group$Sample))
    counts_tmp <- counts[1]
    for (n in 1 : length(samples)) {
        index <- which(names(counts[]) == samples[n])
        counts_tmp <- cbind(counts_tmp, counts[index])
    }
    write.table(counts_tmp, "allsample-counts_tmp.xls", quote = F, row.names = F, sep = "\t")
    data.all <- read.table("allsample-counts_tmp.xls", header = T, sep = "\t", row.names = 1, check.names = F, quote = "")

    group_tmp <- unique(sample_group[,2])
    group <- c()
    for (i in 1:nrow(sample_group)){
    group <- c(group, which(sample_group[i,2] == group_tmp))
    }

    tcc <- new("TCC", data.all, group)
    tcc <- calcNormFactors(tcc, norm.method = "tmm", test.method = "edger", iteration = 1)
    tcc <- estimateDE(tcc, test.method = "bayseq", FDR = FDR)
    result <- getResult(tcc, sort = TRUE)
    result <- as.data.frame(result)
    colnames(result) <- c("gene_id", "a.value", "m.value", "pValue", "qValue","rank","estimatedDEG")
    result_o <- select(result, "gene_id", "pValue", "qValue","rank","estimatedDEG")
    write.table(result_o, paste0("multdiff-all.", Type, ".xls"),sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
    diff_result <- result_o[which(result_o$estimatedDEG == "1") ,]
    write.table(diff_result, paste0("multdiff-diff.", Type, ".xls"),sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
}

##################################################################################
## annotation  function define
add_expression_annotation <- function(sample_group, expression, annotation){

    all <- read.table(paste("multdiff-all.", opt$type, ".xls", sep = ""), header = T, sep = "\t", check.names = F, quote = "")
    diff <- read.table(paste("multdiff-diff.", opt$type, ".xls", sep = ""), header = T, sep = "\t", check.names = F)
    

    samples <- unlist(as.character(sample_group$Sample))

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
    write.table(all, paste("multdiff-all.", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")

    diff.anno <- annotation[match(diff[, 1], annotation[, 1]),]
    diff.expression <- expression.tmp[match(diff[, 1], expression.tmp[, 1]),]
    colnames(diff.expression) <- paste("Expression", colnames(diff.expression), sep = "_")
    diff <- cbind(diff, diff.expression[, 2 : ncol(diff.expression)], diff.anno[, 2 : ncol(diff.anno)])
    write.table(diff, paste("multdiff-diff", ".", opt$type, ".xls", sep = ""), quote = F, row.names = F, sep = "\t", na = "")
    
}

#parameter  import

suppressPackageStartupMessages(library("optparse"))
option_list = list(
make_option(c("-e", "--expression"), type = "character", default = "./gene_fpkm.xls",help = "expression matrix file name", metavar = "character"),
make_option(c("-r", "--readcounts"), type = "character", default ="./gene_counts.xls",help = "read counts matrix file name", metavar = "character"),
make_option(c("-b", "--sample_group"), type = "character", default = "./sample_group.txt",help = "sample group file name", metavar = "character"),
make_option(c("-a", "--annotation"), type = "character", default = "./gene_annotation.xls",help = "annotation file name", metavar = "character"),
#make_option(c("-p", "--pvalue"), type = "double",
#help = "pvalue ratio threshold.Filtering can be performed using any one of (-p), (-f) at a time", metavar = "double"),
make_option(c("-o", "--outputdir"), type = "character", default = "./",help = "Output dir", metavar = "character"),
make_option(c("-f", "--fdr"), type = "double", default = 0.01,help = "fdr ratio threshold.", metavar = "double"),
make_option(c("-t", "--type"), type = "character", default = "gene",help = "output directory for results", metavar = "character"));
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);
#####Rscript Run_DEG.R -e fpkm -r counts -d diff_counts -a annotation -p pval -c FC -t gene -o output####
########################################################################################################################
## get this R script file's own directory
args <- commandArgs(trailingOnly = F)
script.dir <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))


#######import library###################################################################################################
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(oebio))
suppressPackageStartupMessages(library(TCC))
suppressPackageStartupMessages(library(dplyr))
########parameter check
if (is.null(opt$sample_group) |
    is.null(opt$readcounts) |
    is.null(opt$fdr)) {
    print_help(opt_parser)
    stop("read counts matrix file, sample group file, fdr must be supplied", call. = FALSE)
}

if (is.null(opt$expression) | is.null(opt$annotation)) {
    cat("expression matrix file or annotation files is not supplied, annotaion for DEG  will be skipped\n")
}

if (! is.null(opt$expression)) {
    express_file_abpath <- normalizePath(opt$expression)
    expression <- read.table(normalizePath(opt$expression), header = T, check.names = FALSE, quote = "")
    colnames(expression)[1] <- paste(opt$type, "_id")
}

if (! is.null(opt$annotation)) {
    annotation <- read.table(normalizePath(opt$annotation), sep = "\t", header = T, quote = "", na.strings = "", fill = TRUE)
}

########data import#####################################################################################################
counts <- read.table(normalizePath(opt$readcounts), header = T, sep = "\t", check.names = FALSE, quote = "")
sample_group <- read.table(normalizePath(opt$sample_group), header = T, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, quote = "")
express_file_abpath<-normalizePath(opt$expression)
counts_file_abpath <-normalizePath(opt$readcounts)
group_file_abspath <-normalizePath(opt$sample_group)

#####################outputdir check####################################################################################
if (length(unique(sample_group$Group)) < 3) {
    print("the number of group is less than 3, top.")
    q()
}

if (file.exists(opt$outputdir) == "FALSE") {
    dir.create(opt$outputdir)
}
setwd(opt$outputdir)

#######################Method: TCC##################################################################################

if (length(unique(sample_group$Group)) > 2) {
    print("DEG  using TCC is beginning:")
    
    run_TCC(counts, sample_group, opt$fdr, opt$type)
    print("TCC using TCC is OK")
    doParallel::stopImplicitCluster()
}

###################annotation###########################################################################################
if (! is.null(opt$expression) & ! is.null(opt$annotation)) {
    
    print("annotation is  beginning:")
    add_expression_annotation(sample_group, expression, annotation)
    print("annotation is OK")
    doParallel::stopImplicitCluster()
}
