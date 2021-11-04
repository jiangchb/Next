suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

#path1="BA-VS-CC.All.xls"

#path1="BA-VS-CS.All.xls"

path1="CS-VS-CC.All.xls"
DEG <- read.delim(normalizePath(path1), header=T, sep="\t", check.names=F, quote="")

groupname <- gsub("\\.(txt|xls)$", "", gsub(".All.xls$", "", path1))

#=======================改名============================
DEG = plyr::rename(DEG, c("adj_pvalue"="pValue"))
DEG = plyr::rename(DEG, c("log2FC"="log2FoldChange"))

DEG[is.na(DEG)] <- 0
DEG<- subset(DEG,select=c("Protein","log2FoldChange","pValue"))


rownames(DEG)=DEG[,1]
DEG$pValue[which(DEG$pValue < 5E-300 )] = 5E-300

#replace the "-/Inf"
if( "Inf" %in% DEG$log2FoldChange | "-Inf" %in% DEG$log2FoldChange){
  tmp <- DEG[which(DEG$log2FoldChange != "Inf" & DEG$log2FoldChange != "-Inf"),]
  max = max(tmp$log2FoldChange)
  min = min(tmp$log2FoldChange)
  DEG$log2FoldChange <- as.numeric(sub("-Inf", min, DEG$log2FoldChange))
  DEG$log2FoldChange <- as.numeric(sub("Inf", max, DEG$log2FoldChange))
}

#replace the NA
#two$type3 <- factor(two$type3, levels = c('sign.a_down.b_down', 'sign.a_up.b_down', 'sign.a_down.b_up', 'sign.a_up_b_up', 'no.a_down.b_down', 'no.a_up.b_down', 'no.a_down.b_up', 'no.a_up_b_up', 'sign.no', 'no.no'))
#two <- two[order(two$type3, decreasing = FALSE), ]
#preparation by screening
logFC_cutoff = log(2, 2)

DEG$change = as.factor(ifelse(DEG$pValue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff, ifelse(DEG$log2FoldChange > logFC_cutoff ,'Up','Down'),'Filtered'))

DEG <- filter( DEG, -log10(DEG$pValue)<9)
g = ggplot(data=DEG, aes(x=DEG$log2FoldChange, y=-log10(DEG$pValue), color=change)) + 
  geom_point(size =3,alpha=0.5) + 
  theme_classic() + 
  theme(legend.title=element_blank()) + 
  labs(title = paste0( groupname, ": q-value < ","0.05"," && |log2FC| > ", "1" )) + 
  xlab(expression(paste(log[2], " Fold change"))) + 
  ylab(expression(paste("-", log[10], "q-value"))) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  guides(shape=guide_legend(override.aes=list(size=4))) +
  scale_colour_manual(values = c("Up"=c("#8B0000"), "Down"=c("darkblue"),"Filtered"=c("grey")),na.translate=FALSE) +
  theme(panel.grid = element_blank()) +
  #    scale_y_continuous(breaks=seq(0, 10, 2)) +
  scale_x_continuous(limits = c(-8,8),breaks=seq(-8, 8, 2))+ 
  theme(text=element_text(size=20,family="ArialMT")) 
g = g + geom_hline(yintercept = -log(0.05, 10), linetype = "dashed", color = c("grey"), size = 1.3) + geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = c("grey"), size = 1.3)

ggsave(paste0(groupname,".volcano.pdf"),height=10,width=10,plot=g)

