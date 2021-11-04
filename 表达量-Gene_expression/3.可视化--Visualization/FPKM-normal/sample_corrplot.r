library(corrplot)
library(RColorBrewer)
library("optparse")
option_list = list(
        make_option(c("-i", "--input"), type="character", default=NULL,
              help="input expression  file name", metavar="character"),
        make_option(c("-o", "--output"), type="character", default=NULL,
              help="output picture  name(*.pdf)", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$output)){
	print_help(opt_parser)
	stop("--input --output must be supplied", call.=FALSE)
}

data<-read.table(opt$input,header=T,sep="\t",row.names=1, check.names=F,quote="")
matrix<-cor(data,method="pearson")
min<-min(matrix)
wid<-7+2*log2(length(colnames(data)))
hig<-7+2*log2(length(colnames(data)))
#col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
col1 <- colorRampPalette(c("red", "white", "blue"))	
pdf(paste(opt$output,".pdf",sep=""),width = wid, height = hig)
corrplot(matrix, 
   tl.cex=(1.2+0.0001*log2(length(colnames(data)))),
   #title="Correlation Coefficient Between Samples",
   mar=c(0,0,0,0),
   method="circle",
   is.corr=FALSE,
   type="full",
   tl.col="black",
   tl.srt=45,
   cl.lim=c(min,1),
   cl.cex = (1.2+0.0001*log2(length(colnames(data)))),
   addshade="positive",
   rect.col="black",    
   col=col1(255), 
   order="AOE",
   number.cex=(1.0+0.0001*log2(length(colnames(data)))),
   addCoef.col ="black",
   number.digits =4,
   diag= FALSE) 
dev.off()

png(paste(opt$output,".png",sep=""),width = 4000, height = 4000, res = 300)
corrplot(matrix,
   tl.cex=(1.2+0.0001*log2(length(colnames(data)))),
   #title="Correlation Coefficient Between Samples",
   mar=c(0,0,0,0),
   method="circle",
   is.corr=FALSE,
   type="full",
   tl.col="black",
   tl.srt=45,
   cl.lim=c(min,1),
   cl.cex = (1.2+0.0001*log2(length(colnames(data)))),
   addshade="positive",
   rect.col="black",
   col=col1(255),
   order="AOE",
   number.cex=(1.0+0.0001*log2(length(colnames(data)))),
   addCoef.col ="black",
   number.digits =4,
   diag= FALSE)
dev.off()
