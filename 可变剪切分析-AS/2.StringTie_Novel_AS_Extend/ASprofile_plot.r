
#!/usr/bin/env Rscript
library("optparse")
option_list = list(
make_option(c("-i", "--input"), type = "character", default = NULL, help = "input file name", metavar = "character"),
make_option(c("-o", "--outpath"), type = "character", default = NULL, help = "outfile directory", metavar = "character"),
make_option(c("-s", "--sample"), type = "character", default = NULL, help = "sample name", metavar = "character")
);
opt_parser = OptionParser(option_list = option_list, epilogue = "Rscript ASprofile_plot.r -i T1_AS_result.xls -o result/ -s T1");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$outpath)) {
    print_help(opt_parser)
    stop("--input --outpath --sample must be supplied", call. = FALSE)
}
if (! file.exists(opt$outpath)) {dir.create(opt$outpath, recursive = T)}
opt$outpath <- gsub("/$", "", opt$outpath)

library(ggplot2)
library(oebio)
outPDFfile <- paste(opt$outpath, "/", opt$sample, '.AS_stat.pdf', sep = '')
outPNGfile <- paste(opt$outpath, "/", opt$sample, '.AS_stat.png', sep = '')

data <- read.table(opt$input, header = T, sep = "\t", quote = "", check.names = F)
data[, 2] <- gsub("_.+", "", data[, 2])
tab <- as.data.frame(table(data[, 2]))
colnames(tab) <- c("Event_type", "Freq")
max <- max(tab$Freq)
p = ggplot(data = tab, aes(x = Event_type, y = Freq, fill = Event_type)) +
    coord_flip() +
    geom_bar(stat = "identity", position = position_dodge(0.6), width = 0.6) +
    geom_text(aes(label = Freq), hjust = - 0.2, vjust = 0.4, size = 3.5) +
    scale_color_manual(values = alpha(oe_col_qua(1:length(unique(tab$Event_type))),0.8),
    labels = c(as.character(unique(tab$Event_type))))+
    scale_fill_manual(values = oe_col_qua(1 : length(unique(tab$Event_type))),
    labels = c(as.character(unique(tab$Event_type))))+
    guides(fill = FALSE) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(limit = c(0, max + max * 0.1)) +
    xlab("") +
    ylab("") +
    labs(title = paste0("Alternative splicing frequency statistics (", opt$sample, ")")) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 13, family = "ArialMT")) +
    theme(axis.text.y = element_text(size = 10, color = "black")) +
    theme(axis.text.x = element_text(hjust = 1, size = 10, color = "black"))

ggsave(outPDFfile, width = 7, height = 7, plot = p)
ggsave(outPNGfile, type = "cairo-png", width = 7, height = 7, plot = p)

