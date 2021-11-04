#!/usr/bin/env Rscript
library("optparse")
library("viridis")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript pie.r -i anno_tfs -o outdir/");
opt = parse_args(opt_parser);


###############check##############################
if(is.null(opt$input) | is.null(opt$outpath)){
  print_help(opt_parser)
  stop("--input --outpath must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

groupname <- gsub("\\.(txt|xls)$", "", gsub("^Anno_", "", gsub("-diff-pval-0.05-FC-2.gene.xls$", "", basename(opt$input))))
print (groupname)

library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)
library("ggplot2")
top10 <- function(i) { return(i[order(head(i, 10)["Number"],decreasing=T),]) }
expression<-read.table(opt$input, header = TRUE, quote="", sep="\t", check.names=F)
family=as.data.frame(table(expression$Family))
names(family)=c("Species","Number")
family <- family[order(family[,"Number"],decreasing=T),]
data <- top10(family)
print (data)
print (data$Number)
data$Species=factor(data$Species,levels=data$Species);

all<-sum(as.numeric(data$Number))
ratio<-paste(round(as.matrix(as.numeric(data$Number)/all)*100,2),"%",sep="")
myLabel = paste(data$Species,"(",as.numeric(data$Number),",",ratio,")",sep="")
print (myLabel)


###pie()画饼图
pdf(paste0(groupname,"_DEGs_classification_Top10_TFs_family",".pdf"),width=10,height=10)
pie(x=family$Number,labels=paste0(data$Species,"(",as.numeric(data$Number),",",ratio,")",sep=""),col=cividis(length(family$Species)),cex.lab=20,border="white",main=paste0(groupname,":","DEGs classification on Top10 TFs family"))
dev.off()

###显示有问题label总是很密集


##ggplot2画饼图 ---- 会有图例，我不需要这个
#p=ggplot(data, aes(x = "", y = data$Number, fill = data$Species))+
#	geom_bar(stat="identity",color="white")+
#	coord_polar(theta = "y")+ #	
#	labs(x = "", y = "", title = paste0(groupname,":","TFs distribution")) 
#	+theme(axis.text.x = element_blank(),
#        axis.ticks = element_blank(),
#        panel.grid = element_blank())
#	+theme_bw() 
#	+theme(axis.ticks = element_blank()) 
#	+theme(axis.text.x = element_blank()) 

#	+theme(plot.title = element_text(hjust = 0.5)) 
#	+theme(legend.title = element_blank(), legend.position = "right") 
#	+scale_fill_discrete(breaks = data$Species, labels = myLabel) 
#	+theme(panel.grid=element_blank()) 
#	+theme(panel.border=element_blank())


#ggsave(paste0(opt$outpath,"/", "TFs.Piechart.", "pdf"),width=7,height=6,plot=p)
#ggsave(paste0(opt$outpath,'/', "TFs.Piechart.", "png"),type="cairo-png",width=7,height=6,plot=p)
