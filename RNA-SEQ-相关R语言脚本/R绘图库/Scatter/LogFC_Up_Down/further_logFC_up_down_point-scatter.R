##带双 logFC 信息的二维散点图示例
#读取另一个作图数据
##四象限双比较组绘制

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
option_list = list(
    make_option( c("-i", "--path1" ), type = "character",
        help = "The input differential gene file(force).  e.g. *-vs-*-all.gene.xls" ),
    make_option( c("-f", "--path2" ), type = "character",
        help = "The input differential gene file(force).  e.g. *-vs-*-all.gene.xls" ),
    make_option( c("-p", "--pval" ), type = "double",default = 0.05,
        help = "pValue ratio threshold, default: 0.05 ."),
#    make_option( c("-l", "--genelist" ), type = "character" ,
#        help = "Genelist to display the gene symbol.  e.g. genelist.xls"),
#    make_option(c("-t", "--title"), type = "character", default = NULL,
#        help = "Graphic title and outputfile information: Group_A-vs-Group_B . "),
    make_option( c("-o", "--outpath" ),type="character", default = "./Scatter",
        help="the output directory of results, default: ./Scatter ." )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}


#提取组名1
groupname <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", basename(opt$path1)))#opt$input
Case_names <- unlist(strsplit(groupname,"-vs-"))[1]
Control_names <- unlist(strsplit(groupname,"-vs-"))[2]

#提取组名2
groupname2 <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", basename(opt$path2)))#opt$input
Case_names2 <- unlist(strsplit(groupname2,"-vs-"))[1]
Control_names2 <- unlist(strsplit(groupname2,"-vs-"))[2]

one <- read.delim(normalizePath(opt$path1), header=T, sep="\t", check.names=F, quote="")
zero <- read.delim(normalizePath(opt$path2), header=T, sep="\t", check.names=F, quote="")

one = plyr::rename(one, c("log2FoldChange"="logFC_A"))
one = plyr::rename(one, c("q-value"="FDR_A"))
zero = plyr::rename(zero, c("log2FoldChange"="logFC_B"))
zero = plyr::rename(zero, c("q-value"="FDR_B"))

two <- inner_join(one, zero, by=c("gene_id"))

#标记显著性（默认 p < 0.05）
two[which(two$FDR_A < opt$pval & two$FDR_B < opt$pval),'type1'] <- 'sign'
two[which(two$FDR_A >= opt$pval | two$FDR_B >= opt$pval),'type1'] <- 'no'

#标记差异倍数（默认 |log2FC| >= 1）
two[which(two$logFC_A <= -1 & two$logFC_B <= -1),'type2'] <- 'a_down.b_down'
two[which(two$logFC_A >= 1 & two$logFC_B <= -1),'type2'] <- 'a_up.b_down'
two[which(two$logFC_A <= -1 & two$logFC_B >= 1),'type2'] <- 'a_down.b_up'
two[which(two$logFC_A >= 1 & two$logFC_B >= 1),'type2'] <- 'a_up_b_up'
two[is.na(two$type2),'type2'] <- 'no'

#合并显著性和差异倍数，用于标记差异基因
two$type3 <- paste(two$type1, two$type2, sep = '.')

#排序，为了使作图时显著的点绘制在前方（减少被遮盖）
two$type3 <- factor(two$type3, levels = c('sign.a_down.b_down', 'sign.a_up.b_down', 'sign.a_down.b_up', 'sign.a_up_b_up', 'no.a_down.b_down', 'no.a_up.b_down', 'no.a_down.b_up', 'no.a_up_b_up', 'sign.no', 'no.no'))
two <- two[order(two$type3, decreasing = TRUE), ]

#ggplot2 作图，点颜色定义为基因差异类型，点大小表示基因表达量 CPM 值
p <- ggplot(two, aes(logFC_A, logFC_B)) +
  geom_point(aes(color = type3, size = 5), alpha = 0.6, show.legend = FALSE) +
  scale_size(range = c(0, 4)) +
  scale_color_manual(limits = c('sign.a_down.b_down', 'sign.a_up.b_down', 'sign.a_down.b_up', 'sign.a_up_b_up', 'no.a_down.b_down', 'no.a_up.b_down', 'no.a_down.b_up', 'no.a_up_b_up', 'sign.no', 'no.no'), values = c('red', 'orange', 'purple', 'blue', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  geom_vline(xintercept = c(-1, 1), lty = 2) +
  geom_hline(yintercept = c(-1, 1), lty = 2) +
  labs(x = paste0('log2FC',"(",groupname,")"), y = paste0('log2FC',"(",groupname2,')'))

#在合适的位置添加文字标记（当然，选择 AI、PS 后续添加也很方便）
p <- p +
  annotate('text', label = paste0(groupname,' Down\n',groupname2, ' Down'), -4, -8) +
  annotate('text', label = paste0(groupname,' Down\n',groupname2, ' Up'), -4, 8) +
  annotate('text', label = paste0(groupname,' Up\n',groupname2, ' Down'), 4, -8) +
  annotate('text', label = paste0(groupname,' Up\n',groupname2, ' Up'), 4, 8)

ggsave(paste0(opt$outpath,"/",'Scatter-FC2-q-value-0.05.pdf'), p, width = 7, height = 7)
ggsave(paste0(opt$outpath,"/",'Scatter-FC2-q-value-0.05.png'), p, width = 7, height = 7)
print(paste0(opt$outpath, "Scatter", ".png(pdf) is OK"))