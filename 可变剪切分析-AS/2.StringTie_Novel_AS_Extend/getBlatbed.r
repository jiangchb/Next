library("stringr")

c.refmap<-as.matrix(read.delim("c.refmap",header=F,sep="\t"))
r<-c.refmap[,c(2,4)]

rr<-matrix(c("ref_gene_id","cufflinks_gene_id","cufflinks_isoform_id"),ncol=3)
for(i in 1:dim(r)[1]){
	split<-strsplit(r[i,2],",",fixed=T)[[1]]
	
	if(length(split)==1){
		iso<-strsplit(split,"|",fixed=T)[[1]]
		rr<-rbind(rr,cbind(r[i,1],iso[1],iso[2]))		
		}
	else{
		for(j in 1:length(split)){		
		iso<-strsplit(split[j],"|",fixed=T)[[1]]
		rr<-rbind(rr,cbind(r[i,1],iso[1],iso[2]))
		}
		}
	}

rr<-rr[-1,]

merged1<-as.matrix(read.delim("merged.gtf",header=F,sep="\t"))
m1<-merged1[,c(1,4,5,7)]

merged2<-as.matrix(read.table("merged_col9.txt",header=F,sep=";"))
m2<-merged2[,c(1,2,3,4)]

merged<-cbind(m1,m2)
####

###
results1<-matrix(ncol=7,nrow=dim(rr)[1])
colnames(results1)<-c("ref_gene_id","cufflinks_gene_id","cufflinks_isoform_id","chr","start","end","chain")
for(i in 1:dim(rr)[1]){

	ref_gene_id<-rr[i,1]
	gene_id<-rr[i,2]
	isoform_id<-rr[i,3]
	gene_pos<-which(rr[i,2] == merged[,5])
	isoform_pos<-which(rr[i,3] == merged[,6])
	pos<-intersect(gene_pos,isoform_pos)
	
	chr<-unique(merged[pos,1])
	start<-min(as.numeric(merged[pos,2]))
	end<-max(as.numeric(merged[pos,3]))
	chain<-unique(merged[pos,4])
	results1[i,]<-c(ref_gene_id,gene_id,isoform_id,chr,start,end,chain)
	}



gene_name<-as.matrix(read.delim("gene_names.gff",header=F,sep="\t"))

results2<-cbind(results1,NA,NA,NA,NA)

for(j in 1:dim(results1)[1]){
pos<-match(results1[j,1],gene_name[,9])
results2[j,8]<-gene_name[pos,1]  ##ref chr
results2[j,9]<-gene_name[pos,4] ## ref start
results2[j,10]<-gene_name[pos,5] ## ref end
results2[j,11]<-gene_name[pos,7] ## ref strand
}
results2<-results2[which(results2[,8]!="NA"),]
write.table(results2,"results_tmp_Extendgene.txt",row.names=F,col.names=F,quote=F,sep="\t")
######

results3<-read.delim("results_tmp_Extendgene.txt",header=F,sep="\t")
results3<-as.matrix(results3)
p5<-results3[which(results3[,5]<results3[,9]&results3[,6]<results3[,10]),]
p51<-p5[which(as.numeric(p5[,6])-as.numeric(p5[,9]) >= 200),]
p3<-results3[which(results3[,5]>results3[,9]&results3[,6]>results3[,10]),]
p31<-p3[which(as.numeric(p3[,10])-as.numeric(p3[,5]) >= 200),]
Extend.all<-rbind(p51,p31)
write.table(Extend.all,"Extendgene_candidate.txt",row.names=F,col.names=F,quote=F,sep="\t")
##
Extendgene<-Extend.all[,c(4,5,6,3,2,11)]
Extendgene[,2]<-as.numeric(Extendgene[,2])-1
write.table(Extendgene,"Extendgene.bed",row.names=F,col.names=F,quote=F,sep="\t")

Extend.ref<-Extend.all[,c(8,9,10,1,1,11)]
Extend.ref[,2]<-as.numeric(Extend.ref[,2])-1

write.table(Extend.ref,"Extend.ref.gene.bed",row.names=F,col.names=F,quote=F,sep="\t")
