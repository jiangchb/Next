#!/usr/bin/env Rscript
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character"),
	make_option(c("-s", "--sample"), type="character", default=NULL, help="sample name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript run_Saturation.r -i Sample_Eg_2.eRPKM.xls -o outdir/ -s Sample_Eg_2");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$outpath)){
	print_help(opt_parser)
	stop("--input --outpath --sample must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)


library(ggplot2)
library(reshape2)
outPDFfile <- paste(opt$outpath, "/", opt$sample, '_saturation.pdf', sep='')
outPNGfile <- paste(opt$outpath, "/", opt$sample, '_saturation.png', sep='')
  
qq <- c(1,3,15,60)
pp <- seq(5,100,5)
r <- 0.15
data <- read.table(opt$input, header=F, sep="\t", quote="", check.names=F)
dat1 <- subset(data,data[[26]] <= qq[1] & data[[26]] > 0)
dat2 <- subset(data,data[[26]] <= qq[2] & data[[26]] > qq[1])
dat3 <- subset(data,data[[26]] <= qq[3] & data[[26]] > qq[2])
dat4 <- subset(data,data[[26]] <= qq[4] & data[[26]] > qq[3])
dat5 <- subset(data,data[[26]] > qq[4])
  
count1 <- rep(0,20)
num1 <- nrow(dat1)   
rate1 <- rep(0,20)
x <- matrix(0,ncol=20,nrow=num1)
for(i in seq(7,26)) {
  x[ ,i-6] <- dat1[ ,i]/dat1[ ,26]
  for(j in seq(1,num1)) {
    if(abs(x[j,i-6] - 1) <= r) {
      count1[i-6] <- count1[i-6] + 1
    }
  }
  rate1[i-6] <- count1[i-6]/num1
}
  
count2 <- rep(0,20)
num2 <- nrow(dat2)
rate2 <- rep(0,20)
x <- matrix(0,ncol=20,nrow=num2)
for(i in seq(7,26)) {
  x[ ,i-6] <- dat2[ ,i]/dat2[ ,26]
  for(j in seq(1,num2)) {
    if(abs(x[j,i-6]-1) <= r) {
	  count2[i-6] <- count2[i-6] + 1
	}
  }
  rate2[i-6] <- count2[i-6]/num2
}

count3 <-rep(0,20)
num3 <- nrow(dat3)
rate3 <- rep(0,20)
x <-matrix(0,ncol=20,nrow=num3)
for(i in seq(7,26)) {
  x[ ,i-6] <- dat3[ ,i]/dat3[ ,26]
  for(j in seq(1,num3)) {
    if(abs(x[j,i-6]-1) <= r) {
	  count3[i-6] <- count3[i-6] + 1
	}
  }
  rate3[i-6] <- count3[i-6]/num3
}
  
count4 <- rep(0,20)
num4 <- nrow(dat4)
rate4 <- rep(0,20)
x <- matrix(0,ncol=20,nrow=num4)
for(i in seq(7,26)) {
  x[ ,i-6] <- dat4[ ,i]/dat4[ ,26]
  for(j in seq(1,num4)) {
    if(abs(x[j,i-6]-1) <= r) {
	  count4[i-6] <- count4[i-6] + 1
	}
  }
  rate4[i-6] <- count4[i-6]/num4
}

count5 <- rep(0,20)
num5 <- nrow(dat5)
rate5 <- rep(0,20)
x <- matrix(0,ncol=20,nrow=num5)
for(i in seq(7,26)) {
  x[ ,i-6] <- dat5[ ,i]/dat5[ ,26]
  for(j in seq(1,num5)) {
    if(abs(x[j,i-6]-1) <= r) {
	  count5[i-6] <- count5[i-6] + 1
	}
  }
   rate5[i-6] <- count5[i-6]/num5
}

coln1 = paste(0,"~",qq[1])
coln2 = paste(qq[1],"~",qq[2])
coln3 = paste(qq[2],"~",qq[3])
coln4 = paste(qq[3],"~",qq[4])
coln5 = paste(">",qq[4])

df <- cbind(rate1,rate2,rate3,rate4,rate5)
dat <- melt(df)
colnames(dat) <- c('x','rpkm','y')

p <- ggplot(dat,aes(x=x*5,y=y,group=rpkm))+
 geom_line(aes(colour=rpkm),size=1)+
 geom_point(aes(colour=rpkm),size=I(2))+
 scale_x_continuous(breaks=seq(0,100,10),labels=seq(0,100,10))+
 ylim(0,1)+
 xlab("Percent of Mapped reads (%)")+ylab("Fraction of Genes within 15% of Final Values")+
 labs(title=paste0("Saturation of Curve (",opt$sample,")"))+
 scale_colour_discrete(guide=guide_legend(reverse=TRUE),name="FPKM intervals",breaks=c("rate1","rate2","rate3","rate4","rate5"),labels=c(coln1,coln2,coln3,coln4,coln5))
 
 p <- p + theme(panel.border=element_rect(fill=NA,colour="black"))
 p <- p + theme(
         panel.background = element_rect(fill="transparent",colour=NA),
	 panel.grid.minor = element_blank(),
	 panel.grid.major = element_blank(),
	 plot.background = element_rect(fill="transparent",colour=NA),
	 plot.title = element_text(hjust = 0.5, size=13)
 )

ggsave(filename=outPDFfile,plot=p, height=7, width=8)
ggsave(filename=outPNGfile,type="cairo-png",plot=p, height=7, width=8)
