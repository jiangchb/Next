library(getopt);
library(oebio);

spec = matrix(c(
	'fq1','a',0,'character',
	'fq2','b',0,'character',
	'len','c',0,'numeric',
	'name','d',0,'character',
	'key','e',0,'character',
	'od','f',0,'character',
	'help' , 'g', 0, 'logical'
	), byrow=TRUE, ncol=4);
opt = getopt(spec);
print_usage <- function(spec=NULL){
	cat(getopt(spec, usage=TRUE));
	cat("Usage example: \n")
	cat("
Usage example: 
	Rscript QC.r --fq1  --fq2  --len  --key  --od
	
Usage:
	--fq1	fastq1 stat data
	--fq2	fastq2 stat data
	--len	reads length, default=150
	--name  sample name
	--key	rawdata or cleandata, default=cleandata
	--od	output dir
	--help		usage
\n")
	q(status=1);
}

if ( !is.null(opt$help) ) { print_usage(spec) }
if ( is.null(opt$fq1) ) { print_usage(spec) }
if ( is.null(opt$fq2) ) { print_usage(spec) }
if ( is.null(opt$len) ) { opt$len=150 }
if ( is.null(opt$name) ) { print_usage(spec) }
if ( is.null(opt$key) ) { opt$key="cleandata" }
if ( is.null(opt$od) ) { print_usage(spec) }

read1<-read.table(opt$fq1,header=TRUE)
read2<-read.table(opt$fq2,header=TRUE)

quan<-c(read1$mean_qual,read2$mean_qual)
quan<-10^(quan/10*-1)
max_quan=ceiling(max(quan)/0.001)*0.001
Aper<-c(read1$A_count,read2$A_count)
Tper<-c(read1$T_count,read2$T_count)
Cper<-c(read1$C_count,read2$C_count)
Gper<-c(read1$G_count,read2$G_count)
Nper<-c(read1$N_count,read2$N_count)
max_per=ceiling(mean(c(Aper,Cper,Gper,Tper,Nper)))*2
axix_at<-c(1,seq(0,opt$len*2,50)[-1])
temp<-c(1,seq(0,opt$len,50)[-1])
axix_txt<-c(temp,temp[-1])
#axix_at<-c(1,50,100,150,200,250,300)
#axix_txt<-c(1,50,100,150,50,100,150)


pdf(file=paste(opt$od,"/",opt$name,".",opt$key,".qual.pdf",sep=""))
#barplot(quan*100,col='springgreen2',space=0,ylab='Error rate(%)',border= NA,ylim=c(0,max_quan*100),xlab='Reads position(bp)',main='Quality distribution')
barplot(quan*100,col=oe_col_qua(2),space=0,ylab='Error rate(%)',border= NA,ylim=c(0,max_quan*100),xlab='Reads position(bp)',main=paste0('Quality distribution of ',opt$name))
axis(1,labels=axix_txt,at=axix_at)
abline(v=length(quan)/2+1, lty=2,col='darkgray')
box()
dev.off()
png(file=paste(opt$od,"/",opt$name,".",opt$key,".qual.png",sep=""),res=300,width=2000, height=2000)
#barplot(quan*100,col='springgreen2',space=0,ylab='Error rate(%)',border= NA,ylim=c(0,max_quan*100),xlab='Reads position(bp)',main='Quality distribution')
barplot(quan*100,col=oe_col_qua(2),space=0,ylab='Error rate(%)',border= NA,ylim=c(0,max_quan*100),xlab='Reads position(bp)',main=paste0('Quality distribution of ',opt$name))
axis(1,labels=axix_txt,at=axix_at)
abline(v=length(quan)/2+1, lty=2,col=oe_col_qua(8))
box()
dev.off()

pdf(file=paste(opt$od,"/",opt$name,".",opt$key,".base.pdf",sep=""))
plot(Aper,col=oe_col_qua(2),type = 'l',xlab='Reads position(bp)',ylab='Percent(%)',ylim=c(0,max_per),xaxt="n",lty=1,lwd=1.5,main=paste0("Base distribution of ",opt$name))
axis(1,labels=axix_txt,at=axix_at)
abline(v=length(quan)/2+1, lty=2,col='darkgray')
lines(Tper,col=oe_col_qua(3),lty=1,lwd=1.5)
lines(Gper,col=oe_col_qua(4),lty=1,lwd=1.5)
lines(Cper,col=oe_col_qua(5),lty=1,lwd=1.5)
lines(Nper,col=oe_col_qua(6),lty=1,lwd=1.5)
legend("topright",c("A","T","G","C","N"),lty=c(1,1,1,1,1),col=oe_col_qua(2:6))
dev.off()

png(file=paste(opt$od,"/",opt$name,".",opt$key,".base.png",sep=""),res=300,width=2000, height=2000)
plot(Aper,col=oe_col_qua(2),type = 'l',xlab='Reads position(bp)',ylab='Percent(%)',ylim=c(0,max_per),xaxt="n",lty=1,lwd=1.5,main=paste0("Base distribution of ",opt$name))
axis(1,labels=axix_txt,at=axix_at)
abline(v=length(quan)/2+1, lty=2,col='darkgray')
lines(Tper,col=oe_col_qua(3),lty=1,lwd=1.5)
lines(Gper,col=oe_col_qua(4),lty=1,lwd=1.5)
lines(Cper,col=oe_col_qua(5),lty=1,lwd=1.5)
lines(Nper,col=oe_col_qua(6),lty=1,lwd=1.5)
legend("topright",c("A","T","G","C","N"),lty=c(1,1,1,1,1),col=oe_col_qua(2:6))
dev.off()

print('all done!')
q()
