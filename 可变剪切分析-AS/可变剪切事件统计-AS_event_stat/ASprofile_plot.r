
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

path1 <- "LINC00278_AS_result.xls"
path2 <- "PCDNA3_1_AS_result.xls"

data <- read.table(path1, header = T, sep = "\t", quote = "", check.names = F)
data[, 2] <- gsub("_.+", "", data[, 2]) #去除每个事件后面的后缀_on,_off
tab <- as.data.frame(table(data[, 2]))#对于第二列也就是type那一列进行计数统计
colnames(tab) <- c("Event_type", "Freq")
max <- max(tab$Freq)

data2 <- read.table(path2, header = T, sep = "\t", quote = "", check.names = F)
data2[, 2] <- gsub("_.+", "", data2[, 2]) #去除每个事件后面的后缀_on,_off
tab2 <- as.data.frame(table(data2[, 2]))#对于第二列也就是type那一列进行计数统计
colnames(tab2) <- c("Event_type", "Freq")
max <- max(tab2$Freq)

#===================================计算总数=======================#
colsum2=sum(tab2$Freq)
colsum1=sum(tab$Freq)

tab["Sum"]=colsum1
tab2["Sum"]=colsum2

#===============================标记组名========================#
tab["Sample"]="LINC00278"
tab2["Sample"]="PCDNA3_1"

bindtab=rbind(tab,tab2)

write.table(bindtab, file=("./LINC00278-vs-PCDNA3_1.AS_stat.xls"), sep="\t", quote=FALSE,
            col.names=TRUE, row.names=FALSE)

#=============================绘制多个样本堆叠分布=================================#

color_tep<-c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray',"darkblue")
library(ggthemes)
library(RColorBrewer)
library(ggsci)
color1 <- brewer.pal(12, "Set3")
p1=ggplot(bindtab,aes(x=Sample,weight=Freq,fill=Event_type))+
  geom_bar(position="stack",width=0.3)+labs(title="Distribution of AS",y="Number",fill="AS Event type")+
  theme(plot.title = element_text(hjust = -2))+
#  scale_color_manual(values = color1,
#                    labels = c(as.character(unique(bindtab$Event_type))))+
#  scale_fill_manual(values = color1,
#                     labels = c(as.character(unique(bindtab$Event_type))))+
#  geom_text(label = paste(bindtab$Freq)) +
  theme_classic()+  
  scale_fill_simpsons()+
  theme(legend.position="right")+
#  theme(legend.text=element_text(size=20))+
  theme(plot.title = element_text(size=20,face = "bold"),axis.text.y=element_text(size=15),axis.text.x=element_text(size=10),axis.title.y=element_text(size=15,face="bold"),axis.title.x=element_text(size=15,face="bold")) 
p1

outPDFfile <- paste( "./","LINC00278-vs-PCDNA3_1", '.AS_stat3.pdf', sep = '')
outPNGfile <- paste( "./","LINC00278-vs-PCDNA3_1", '.AS_stat3.png', sep = '')

ggsave(outPDFfile, width = 7, height = 7, plot = p1)
ggsave(outPNGfile, type = "cairo-png", width = 7, height = 7, plot = p1)



#==========================绘制单个样本事件分布======================================#
p = ggplot(data = tab, aes(x = Event_type, y = Freq, fill = Event_type)) +
    coord_flip() +
    geom_bar(stat = "identity", position = position_dodge(0.6), width = 0.6) +
    geom_text(aes(label = Freq), hjust = - 0.2, vjust = 0.4, size = 3.5) +
#    scale_color_manual(values = alpha(oe_col_qua(1:length(unique(tab$Event_type))),0.8),
#    labels = c(as.character(unique(tab$Event_type))))+
#    scale_fill_manual(values = oe_col_qua(1 : length(unique(tab$Event_type))),
#    labels = c(as.character(unique(tab$Event_type))))+
    guides(fill = FALSE) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(limit = c(0, max + max * 0.1)) +
    xlab("") +
    ylab("") +
    labs(title = paste0("Alternative splicing frequency statistics (", "LINC00278", ")")) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 1, size = 13, family = "serif")) +
    theme(axis.text.y = element_text(size = 10, color = "black")) +
    theme(axis.text.x = element_text(hjust = 1, size = 10, color = "black"))
p

ggsave(outPDFfile, width = 7, height = 7, plot = p)
ggsave(outPNGfile, type = "cairo-png", width = 7, height = 7, plot = p)

