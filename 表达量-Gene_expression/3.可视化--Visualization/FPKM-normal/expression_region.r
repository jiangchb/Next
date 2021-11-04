#!/usr/local/bin/r
library(ggplot2)
library(grid)
library("optparse") 
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input density_bar  file name", metavar="character"),
	make_option(c("-d", "--group"), type="character", default=NULL,
	      help="group file name", metavar="character"),
	make_option(c("-u", "--type"), type="character", default="FPKM", help="fpkm or rpm, default=FPKM", metavar="character"),
	make_option(c("-w", "--width"), type="integer", default=7, help="picture width, default=7", metavar="integer"),
	make_option(c("-t", "--height"), type="integer", default=6, help="picture height, default=6", metavar="integer")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
	print_help(opt_parser)
	stop("--input must be supplied", call.=FALSE)
}

ddd<-read.delim(opt$input,header=TRUE,quote="")
if(! is.null(opt$group)){
	group<-read.delim(opt$group,header=TRUE,quote="")
	gp<-as.character(group$Sample)
	print(gp)
	ddd$Sample=factor(ddd$Sample,levels=gp)
}else{
	ddd$Sample=factor(ddd$Sample)
}

### expression_region
#p=ggplot(ddd, aes(x=Sample, fill=expression_region)) +
p=ggplot(ddd, aes(x=Sample, y=number, fill=expression_region)) +
  geom_bar(stat = 'identity', position="stack",width=0.7)+
  geom_text(aes(label=ddd$number),hjust=0.5, vjust=-0.5, size=5, position = position_stack())+
  theme_bw()+
  theme(legend.title=element_blank())+
  theme(axis.text.x=element_text(angle=45,color="black",size=12,hjust=1,vjust = 1)) +
  theme(axis.text.y=element_text(angle = 00, size=12, color="black"))+
  theme(panel.grid =element_blank())+
  theme(legend.key.height=unit(0.5,"cm"),legend.text=element_text(size=12))+
  theme(plot.title = element_text(hjust = 0.5, size=15))+
  ylab("Transcript number")+xlab("")+labs(title = "Transcript expression in each sample")
ggsave(paste(opt$type,"_region.pdf",sep=""), width=opt$width, height=opt$height, plot=p)
ggsave(paste(opt$type,"_region.png",sep=""), type="cairo-png", width=opt$width, height=opt$height, plot=p)
print(paste(opt$type,"_region.pdf(png) is OK",sep=""))

