#!/usr/local/bin/r
library(ggplot2)
library(grid)
library("optparse") 
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input density_bar  file name", metavar="character"),
	make_option(c("-d", "--group"), type="character", default=NULL,
	      help="group file name", metavar="character"),
	make_option(c("-u", "--type"), type="character", default="FPKM",
	      help="fpkm or rpm, default=fpkm", metavar="character")
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

### expression_density
max_exp<-max(ddd$expression)
min_exp<-min(ddd$expression)
xlimit_max<-log10(max_exp)+5
xlimit_min<--xlimit_max
if(max_exp<xlimit_max && min_exp>xlimit_min){
	xlimit_max<-max_exp
	xlimit_min<-min_exp
}
#ddd$expression[ddd$expression==0]=min(ddd$expression[ddd$expression!=0])*0.001
#p=ggplot(ddd, aes(x =log10(expression),fill=Sample))+
p=ggplot(ddd, aes(x =log10(expression),colour=Sample))+
  geom_density(alpha =0.3)+
#  coord_cartesian(xlim=c(xlimit_min, xlimit_max))+
  theme_bw()+
  theme(panel.grid =element_blank())+
  theme(axis.text=element_text(size=9))+
  theme(axis.title=element_text(size=13))+
  labs(title = "Transcript expression density in each sample")+
  theme(legend.text= element_text(size=9,color="black", vjust=0.5, hjust=0.5))+
  theme(legend.title= element_text(size=9))+theme(plot.title = element_text(hjust = 0.5,size=13))+
  xlab(paste("log10(",opt$type,")",sep=""));
ggsave(paste(opt$type,"_density.pdf",sep=""), width=12, height=7, plot=p)
ggsave(paste(opt$type,"_density.png",sep=""), type="cairo-png", width=12, height=7, plot=p)
print(paste(opt$type,"_density.pdf(png) is OK",sep=""))
