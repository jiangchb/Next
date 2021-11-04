suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("oebio"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("corrplot"))
option_list = list(
make_option(c("-i", "--input"), type = "character", default = NULL,
help = "input expression  file name", metavar = "character"),
make_option(c("-s", "--group"), type = "character", default = NULL,
help = "input group  file name (optional)", metavar = "character"),
make_option(c("-t", "--type"), type = "character", default = NULL,
help = "input expression type:FPKM,TPM,RPM....", metavar = "character"),
make_option(c("-o", "--outdir"), type = "character", default = NULL,
help = "output directory", metavar = "character")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

data <- read.table(opt$input, header = T, sep = "\t", check.names = F, quote = "")

convert1 <- function(d) {
    lapply(2 : ncol(d), function(i) {
        d2 <- d[, c(1, i)]
        d2$type = colnames(d2)[2]
        colnames(d2) = c("Gene", "Expression", "Sample")
        return(d2)
    }) %>% do.call('rbind', .)
}

convert2 <- function(d, g) {
    lapply(2 : ncol(d), function(i) {
        d2 <- d[, c(1, i)]
        d2$type = colnames(d2)[2]
        d2$group = as.character(group[which(group$Sample == colnames(data[, c(1, i)])[2]),]$Group)
        colnames(d2) = c("Gene", "Expression", "Sample", "Group")
        return(d2)
    }) %>% do.call('rbind', .)
}

f <- function(x) {
    r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
}

q <- function(x) {
    subset(x, quantile(x, 0.95) < x)
}

if (! is.null(opt$group)) {show_type = "Group"}else {show_type = "Sample"}



if (is.null(opt$group)) {
    data<-data[ ,c(as.character(colnames(data)[1]), sort(as.character(colnames(data)[-1])))]
    print(head(data))
    MyData <- convert1(data)
    p<-ggplot(MyData, aes(factor(Sample), log10(Expression + 1),
    fill = factor(Sample))) +
        stat_summary(fun.data = f, geom = "boxplot",
        position = position_dodge(1)) +
        stat_summary(aes(color = factor(Sample)), fun.y = q, geom = "point",
        position = position_dodge(1)) +
        scale_color_manual(values = oe_col_qua(1 : length(unique(MyData$Sample))),
        labels = c(unique(MyData$Sample)),
        name = "Sample") +
        scale_fill_manual(values = oe_col_qua(1 : length(unique(MyData$Sample))),
        labels = c(unique(MyData$Sample)),
        name = "Sample") +
        theme_bw() +
        labs(x = "Sample") +
        ylab(paste0("log10(", opt$type , "+1)")) +
        labs(title = paste0("Boxplot for " , opt$type , " Values")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.margin = unit(c(0.5, 0, 0, 2), "cm"))
    ggsave(paste0(opt$type,"_boxplot.pdf"), height = 6, width = 8+0.2 * length(colnames(data)),plot=p)
    ggsave(paste0(opt$type,"_boxplot.png"), type="cairo-png",height = 6, width = 6+0.2 * length(colnames(data)),plot=p)
}else
{   group <- read.table(opt$group, header = T, sep = "\t", check.names = F, quote = "")
    data<-data[ ,c(as.character(colnames(data)[1]),as.character(group$Sample))]
    MyData <- convert2(data, group)
    group_labels<-c(as.character(sort(unique(MyData$Group))))
    
    write.table(MyData, paste0("MyData",".xls"), row.names = F, quote = F, sep = "\t")

    p<-ggplot(MyData, aes(Sample, log10(Expression + 1),
    fill = Group)) +
        stat_summary(fun.data = f, geom = "boxplot",
        position = position_dodge(1)) +
        stat_summary(aes(color = factor(Group)), fun.y = q, geom = "point",
        position = position_dodge(1)) +
        facet_grid(. ~ Group, scales = "free") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        scale_color_manual(values = oe_col_qua(1 : length(group_labels)), labels = group_labels,name = "Group") +
        scale_fill_manual(values = oe_col_qua(1 : length(group_labels)), labels = group_labels,name = "Group") +
        labs(x = "") +
        theme_bw() +
        ylab(paste0("log10(", opt$type , "+1)")) +
        labs(title = paste0("Boxplot for " , opt$type , " Values")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, color = "black")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(plot.margin = unit(c(0.2, 2, 0.2, 2), "cm"))
    ggsave(paste0(opt$type,"_boxplot.pdf"), height = 6, width = 8+0.2 * length(colnames(data)),plot=p)
    ggsave(paste0(opt$type,"_boxplot.png"), type="cairo-png",height = 6, width = 6+0.2 * length(colnames(data)),plot=p)
}


d <- read.table(opt$input, header = T, sep = "\t", row.names = 1, check.names = F, quote = "")
m <- dim(d)[2]
col <- colnames(d)
expression.stat <- matrix(NA, ncol = 9, nrow = m)
colnames(expression.stat) <- c("Sample", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max", "Sd", "Sum")
for (i in 1 : m) {
    data <- d[[i]]
    stats <- boxplot.stats(data)$stats
    expression.stat[i, 1] <- col[i]
    expression.stat[i, 3] <- stats[2]
    expression.stat[i, 6] <- stats[4]
    expression.stat[i, 5] <- mean(data)
    expression.stat[i, 4] <- median(data)
    expression.stat[i, 2] <- min(data)
    expression.stat[i, 7] <- max(data)
    expression.stat[i, 9] <- sum(data)
    expression.stat[i, 8] <- sd(data)
}
write.table(expression.stat, paste0("report.", opt$type, ".xls"), row.names = F, quote = F, sep = "\t")
