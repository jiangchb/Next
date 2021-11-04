### get blat_results.txt

blat<-as.matrix(read.delim("blat_results.txt",header=F,sep="\t"))
ref.gene.info<-as.matrix(read.delim("Extend.ref.gene1.bed",header=F,sep="\t"))
gene.info<-as.matrix(read.delim("Extendgene1.bed",header=F,sep="\t"))

results4<-matrix(c("gene","5'or 3'","chr","strand","original region","extenden region"),ncol=6)
for(i in 1:dim(blat)[1]){
	if(as.numeric(blat[i,3])==0){
		gene<-blat[i,1]
		info<-which(blat[i,1]==ref.gene.info[,4])
		original.region<-paste(ref.gene.info[info,2],"-",ref.gene.info[info,3],collapse="")
		chr<-ref.gene.info[info,1]
		strand<-ref.gene.info[info,6]

		info1<-	which(blat[i,5]==gene.info[,4])	
		extenden.region<-paste(gene.info[info1,2],"-",gene.info[info1,3],collapse="")
		
		#if(strand=="+"){prime<-"5"}
		#if(strand=="-"){prime<-"3"}
		prime<-"5"
		
		results4<-rbind(results4,cbind(gene,prime,chr,strand,original.region,extenden.region))
		}
	if(as.numeric(blat[i,2])==as.numeric(blat[i,4])){
		gene<-blat[i,1]
		info<-which(blat[i,1]==ref.gene.info[,4])
		original.region<-paste(ref.gene.info[info,2],"-",ref.gene.info[info,3],collapse="")
		chr<-ref.gene.info[info,1]
		strand<-ref.gene.info[info,6]

		info1<-	which(blat[i,5]==gene.info[,4])	
		extenden.region<-paste(gene.info[info1,2],"-",gene.info[info1,3],collapse="")
		
		#if(strand=="+"){prime<-"3"}
		#if(strand=="-"){prime<-"5"}
		prime<-"3"
		
		results4<-rbind(results4,cbind(gene,prime,chr,strand,original.region,extenden.region))
		}
	
	}

write.table(results4,"extendgene_results.xls",quote=F,col.names=F,row.names=F,sep="\t")
