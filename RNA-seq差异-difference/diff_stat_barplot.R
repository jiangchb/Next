suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(oebio))
suppressPackageStartupMessages(library(ggrepel))

##diff_stat_barplot
diff_stat_barplot <- function(diff_stat, type){
  library(ggplot2)
  library(stringr)
  data <- read.table(diff_stat, header = T, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = F);
  for (i in 1 : length(data$Control))
  {
    if (! is.na(str_match(data$Case[i], "\\(")[1] == "("))
    {
      data$Control[i] <- substr(data$Control[i], str_locate_all(data$Control[i], "\\(")[[1]] + 1, str_locate_all(data$Control[i], "\\)")[[1]] - 1)
    }
    if (! is.na(str_match(data$Case[i], "\\(")[1] == "(")) {
      data$Case[i] <- substr(data$Case[i], str_locate_all(data$Case[i], "\\(")[[1]] + 1, str_locate_all(data$Case[i], "\\)")[[1]] - 1)
    }
  }
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
    geom_bar(stat = "identity", position = "dodge",width= 0.9 ) +
    scale_color_manual(values = rev(oe_col_qua(1 : length(unique(data2$Type)))),
                       labels = c(as.character(unique(factor(data2$Type)))),
                       name = "Type") +
    scale_fill_manual(values =rev(oe_col_qua(1: length(unique(data2$Type)))),
                      labels = c(as.character(unique(factor(data2$Type)))),
                      name = "Type") +
    geom_text_repel(aes(label = data2$Gene_number), force =1 ,size = 3, colour = 'black', vjust = -0.6, hjust = 0.5, position = position_dodge(0.4)) +
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
  ggsave(paste0(type, "_diff_stat_barplot.pdf"), height = 6+0.5*length(colnames(data2)), width = 8+0.5*length(colnames(data2)), plot = p)
  ggsave(paste0(type, "_diff_stat_barplot.png"), type = "cairo-png", height = 6, width = 8+0.2*length(colnames(data2)), plot = p)
  file.remove(paste0(type, "_diff_stat_barplot.xls"))
}

suppressPackageStartupMessages(library("optparse"))
option_list = list(
make_option(c("-t", "--type"), type = "character", default = NULL,
help = "transcript or gene type:mRNA,circRNA,lncRNA,miRNA,gene", metavar = "character"));
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

## get this R script file's own directory
args <- commandArgs(trailingOnly = F)
script.dir <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))


diff_stat_barplot(paste0(opt$type, "_diff_stat.xls"), opt$type)
