library("corrplot")

path1 <- "AllSamples.correlation.xls"
d <- read.delim(path1, sep="\t",row.names = 1, header=T, quote="")

col1 <- colorRampPalette(c("red", "white", "blue"))	
pdf(paste(opt$output,".pdf",sep=""),width = wid, height = hig)

因为不适合label所以弃用
corrplot(matrix, 
         tl.cex=(1+0.0001*log2(length(colnames(data)))),
         #title="Correlation Coefficient Between Samples",
         mar=c(0,0,0,0),
         method="circle",
         is.corr=FALSE,
         type="upper",
         tl.col="black",
         tl.srt=45,
         cl.lim=c(min,1),
         cl.cex = (1.2+0.0001*log2(length(colnames(data)))),
         addshade="positive",
         rect.col="black",    
         col=col1(255), 
         order="AOE",
         number.cex=0.8,
         number.digits =4,
         diag= FALSE, ) 

dev.off()
