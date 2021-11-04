library("stringr")

##############
u.tmap<-as.matrix(read.delim("u.tmap",header=F,sep="\t"))
r<-u.tmap[,c(4,5,11)]

merged1<-as.matrix(read.delim("merged.gtf",header=F,sep="\t"))
m1<-merged1[,c(1,4,5,7)]

merged2<-as.matrix(read.table("merged_col9.txt",header=F,sep=";"))
m2<-merged2[,c(1,2,3,4)]

merged<-cbind(m1,m2)

results1<-matrix(c("chr","gene_id","isoform_id","+/-","blocks","Sizes","Starts","locus_start","locus_end"),ncol=9)
for(i in 1:dim(r)[1]){
	pos_gene<-which(r[i,1]==merged[,5])
	pos_isoform<-which(r[i,2]==merged[,6])
	pos<-intersect(pos_gene,pos_isoform)
	chr<-unique(merged[pos,1])
	gene_id<-unique(merged[pos,5])
	isoform_id<-unique(merged[pos,6])
	chain<-unique(merged[pos,4])
	Blocks<-length(pos)
	Sizes<-str_c(as.character(as.numeric(merged[pos,3])-as.numeric(merged[pos,2])),collapse="///")
	Starts<-str_c(as.character(as.numeric(merged[pos,2])),collapse="///")
        
	if(chain=="."){next}
	else{
	locus_start<-min(c(as.numeric(merged[pos,3]),as.numeric(merged[pos,2])))
	locus_end<-max(c(as.numeric(merged[pos,3]),as.numeric(merged[pos,2])))
	}
	
	results1<-rbind(results1,cbind(chr,gene_id,isoform_id,chain,Blocks,Sizes,Starts,locus_start,locus_end))
	}


write.table(results1[,1:7],"results_novel_isoform.txt",row.names=F,col.names=F,quote=F,sep="\t")

write.table(results1,"results_novel_isoform1.txt",row.names=F,col.names=F,quote=F,sep="\t")

results1<-results1[-1,]
bed<-results1[,c(1,8,9,3,2,4)]
bed[,2]<-as.numeric(bed[,2])-1
write.table(bed,"novel_isoform1.bed",row.names=F,col.names=F,quote=F,sep="\t")


gene.bed<-read.delim("gene.bed",header=F,sep="\t")
gene.bed1<-gene.bed

pos<-which(gene.bed[[2]]<=200)
gene.bed1[pos,2]<-0

pos<-which(gene.bed[[2]]>200)
gene.bed1[pos,2]<-gene.bed[pos,2]-200

gene.bed1[[3]]<-gene.bed[[3]]+200

write.table(gene.bed1,"gene1_200.bed",row.names=F,col.names=F,quote=F,sep="\t")
