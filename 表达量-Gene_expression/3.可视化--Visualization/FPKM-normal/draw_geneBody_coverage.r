#!/usr/bin/env Rscript
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character"),
	make_option(c("-s", "--sample"), type="character", default=NULL, help="sample name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list, epilogue="Rscript draw_geneBody_coverage.r -i Sample_Eg_2.geneBodyCoverage.txt -o outdir/ -s Sample_Eg_2");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$outpath)){
	print_help(opt_parser)
	stop("--input --outpath --sample must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)


library(ggplot2)
outPDFfile <- paste(opt$outpath, "/", opt$sample, '.geneBodyCoverage.pdf', sep='')
outPNGfile <- paste(opt$outpath, "/", opt$sample, '.geneBodyCoverage.png', sep='')
  
data <- read.table(opt$input, header=F, sep="\t", quote="", check.names=F,row.names=1)
dat<-data.frame(t(data))
colnames(dat) <- c("Percentile","Coverage")
dat$Coverage=(dat$Coverage-min(dat$Coverage))/(max(dat$Coverage)-min(dat$Coverage))

p <- ggplot(dat,aes(x=Percentile, y=Coverage))+
 geom_line(size=0.5,col="#7FC97F")+
 scale_x_continuous(breaks=seq(0,100,20),labels=seq(0,100,20))+
 #ylim(0,1)+
 xlab("Gene body percentile (5'->3')")+ylab("Coverage")+
 labs(title=paste0("Coverage of Curve (",opt$sample,")"))+
 guides(fill=FALSE)
 
 p <- p + theme(panel.border=element_rect(fill=NA,colour="black"))
 p <- p + theme(
         panel.background = element_rect(fill="transparent",colour=NA),
         panel.grid.minor = element_blank(),
	     panel.grid.major = element_blank(),
         plot.background = element_rect(fill="transparent",colour=NA),
         plot.title = element_text(hjust = 0.5, size=13)
 )

ggsave(filename=outPDFfile,plot=p, height=6.5, width=7)
ggsave(filename=outPNGfile,type="cairo-png",plot=p, height=6.5, width=7)
