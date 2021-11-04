

#这边需要加上逻辑判断输入多少个样本，并且生成VennData

data <- read.table("./VennData.xls",sep="\t",header=T,quote="",check.names=FALSE)
library("UpSetR")
datnum <- length(colnames(data))-1
pdf("./VennGraph.pdf", width=6.5, height=6.5, onefile=F)
upset(data, nsets=datnum, nintersects=NA, number.angles=0,point.size=1, line.size=0.5, mainbar.y.label="Intersection gene number", sets.x.label="Total gene number", text.scale=c(1.3,1.3,1,1,1.3,1), mb.ratio=c(0.55,0.45), order.by="freq", show.numbers="yes", sets.bar.color=c("red"),main.bar.color = "blue")
dev.off()

png("./VennGraph.png", width=4200, height=2100, res=300)
upset(data, nsets=datnum, nintersects=NA, number.angles=30, point.size=2, line.size=1, mainbar.y.label="Intersection gene number", sets.x.label="Total gene number", text.scale=c(1.3,1.3,1,1,1.3,1), mb.ratio=c(0.55, 0.45), order.by="freq", show.numbers="yes", sets.bar.color=rainbow(datnum),main.bar.color = "blue")
dev.off()

#queries = list(list(query = intersects,params = list('A','B','C'), color = "orange", active = T))) 最多两个