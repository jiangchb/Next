#!/usr/bin/env Rscript
# 并列方法
library("optparse")
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
  make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character"),
  make_option(c("-n", "--Number"), type="character", default=10, help="The number of terms you want to conserve", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top****.r -i enrichment-kegg-Group1-vs-Group2-Down.txt -m Down -n 10 -o outdir/
                          ");
opt = parse_args(opt_parser);
if(is.null(opt$input) | is.null(opt$outpath) | is.null(opt$mark)){
  print_help(opt_parser)
  stop("--input --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)


library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-kegg-", "", basename(opt$input)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Down$", "(Down)", groupname)
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Up$", "(Up)", groupname)
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
  groupname <- gsub("-Total$", "(Total)", groupname)
}

Number <-as.numeric(opt$Number)

print ("We will do top ")
print (Number)

#对字数的限制
limit70 <- function(s) {
  k <- as.character(s)
  if(str_length(s)>70){k <- sub("[^ ]+$", "...", substr(k,1,67))}
  return(k)	
}

#########创建一个新列，新列是要把同一个term的p值相加的#########
d <- read.delim(opt$input, sep="\t", header=T, quote="")
attach(d)
#根据原文第11列生成相同term的和
d2<-aggregate(d[5], by=list(Term), FUN=sum)
#重新命名
names(d2)<-c('Term','ListHits_sum') 
#把这个数据重新赋值到原来的表下面
#与之前的数据合并
d3=merge(d,d2,by.x="Term",by.y="Term",all=T)
#先筛选和排序
#d3 <- d3[which(d3[,"ListHits"]>2),]
d3 <- d3[order(d3[,"ListHits_sum"],decreasing=T),]
##########筛选出top多少的term名称###############
#定义函数
top20 <- function(i) { return(i[order(head(i, 20)["ListHits_sum"],decreasing=T),]) }
#先筛选和排序

top10 <- function(i) { return(i[order(head(i, Number)["ListHits_sum"],decreasing=T),]) }
#除了p值外，其他都是升序
#stopifnot(nrow(d)>0)
if(nrow(d)==0){
  print("d items = 0, program exit!")
  q()
}
#分图##################输出参考表格################
dp <- rbind(top20(d3[d3["Category"]=="biological_process", ]),
            top20(d3[d3["Category"]=="cellular_component", ]),
            top20(d3[d3["Category"]=="molecular_function", ]))
write.table(dp[,c(1,2,3,4,5,6,7,8,9,10,11,12,13)], paste0(opt$outpath, "/GO(BPCCMF).top60.", opt$mark, ".xls"), sep="\t", quote=FALSE,
            col.names=TRUE, row.names=FALSE)

########################可能没用#####################
#创建Up,down 各top10


#########################BP画图######################
d1<-d3[which(d3[,3]=="biological_process"),]
d1_l<-d1[order(d1[,13],decreasing=T),]
d1_s <- rbind(top10(d1_l[d1_l["Regulation"]=="Up", ]),
             top10(d1_l[d1_l["Regulation"]=="Down", ]))
d1_s["Term"] <- apply(d1_s["Term"], 1, limit70)
uniq1<-d1_s$Term[!duplicated(d1_s$Term)] 
d1_s$Term <- factor(d1_s$Term, levels=uniq1)

p=ggplot(data=d1_s, aes(x=Term, y=ListHits, width=0.8, fill=Regulation,space=0)) +
  coord_flip()+
  geom_bar(stat="identity",position="dodge") +
  scale_x_discrete(breaks=d1_s$Term, labels=d1_s$Term)+ 
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_brewer(palette = "Set1")+
  labs(x="", y='Number of gene', title = paste0(groupname,": ","Top ",opt$Number, "(Biological Process) Up & Down go Term"))+  
  theme_bw() + scale_colour_gradient(low = "green", high = "red") +
  theme(axis.text.y=element_text(size=11,color="black")) + 
  theme(axis.text.x=element_text(hjust=1, size=14))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.position="right")+
  theme(legend.text=element_text(size=11))+
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,4,2,18), "lines"))+
  theme(panel.grid =element_blank()) 

ggsave(paste0(opt$outpath, "/GO(Biological_process).top.", "pdf"), height=10, width=18, plot=p)
ggsave(paste0(opt$outpath, "/GO(Biological_process).top.", "png"), type="cairo-png", height=10, width=18, plot=p)
print(paste0(opt$outpath, "/GO(Biological_process).top.", "png(pdf) is OK"));


#########################cc画图######################

d2<-dp[which(dp[,3]=="cellular_component"),]
d2_l<-d2[order(d2[,13],decreasing=T),]
d2_s <- rbind(top10(d2_l[d2_l["Regulation"]=="Up", ]),
              top10(d2_l[d2_l["Regulation"]=="Down", ]))
d2_s["Term"] <- apply(d2_s["Term"], 1, limit70)
uniq2<-d2_s$Term[!duplicated(d2_s$Term)]
d2_s$Term <- factor(d2_s$Term, levels=uniq2)

p2=ggplot(data=d2_s, aes(x=Term, y=ListHits, width=0.8, fill=Regulation,space=0)) +
  coord_flip()+
  geom_bar(stat="identity",position="dodge") +
  scale_x_discrete(breaks=d2_s$Term, labels=d2_s$Term)+ 
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_brewer(palette = "Set1")+
  labs(x="", y='Number of gene', title = paste0(groupname,": ","Top ",opt$Number, "(Cellular component) Up & Down go Term"))+  
  theme_bw() + scale_colour_gradient(low = "green", high = "red") +
  theme(axis.text.y=element_text(size=11,color="black")) + 
  theme(axis.text.x=element_text(hjust=1, size=14))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.position="right")+
  theme(legend.text=element_text(size=11))+
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,4,2,18), "lines"))+
  theme(panel.grid =element_blank()) 
ggsave(paste0(opt$outpath, "/GO(Cellular_component).top.", "pdf"), height=10, width=18, plot=p2)
ggsave(paste0(opt$outpath, "/GO(Cellular_component).top.", "png"), type="cairo-png", height=10, width=18, plot=p2)
print(paste0(opt$outpath, "/GO(Cellular_component).top.", "png(pdf) is OK"));

######################################画MF################################
d3<-dp[which(dp[,3]=="molecular_function"),]
d3_l<-d3[order(d3[,13],decreasing=T),]
d3_s <- rbind(top10(d3_l[d3_l["Regulation"]=="Up", ]),
              top10(d3_l[d3_l["Regulation"]=="Down", ]))
d3_s["Term"] <- apply(d3_s["Term"], 1, limit70)
uniq3<-d3_s$Term[!duplicated(d3_s$Term)]
d3_s$Term <- factor(d3_s$Term, levels=uniq3)


p3=ggplot(data=d3_s, aes(x=Term, y=ListHits, width=0.8, fill=Regulation,space=0)) +
  coord_flip()+
  geom_bar(stat="identity",position="dodge") +
  scale_x_discrete(breaks=d3_s$Term, labels=d3_s$Term)+ 
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_brewer(palette = "Set1")+
  labs(x="", y='Number of gene', title = paste0(groupname,": ","Top ",opt$Number, "(Molecular function) Up & Down go Term"))+  
  theme_bw() + scale_colour_gradient(low = "green", high = "red") +
  theme(axis.text.y=element_text(size=11,color="black")) + 
  theme(axis.text.x=element_text(hjust=1, size=14))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.position="right")+
  theme(legend.text=element_text(size=11))+
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,4,2,18), "lines"))+
  theme(panel.grid =element_blank()) 

ggsave(paste0(opt$outpath, "/GO(Molecular_function).top.", "pdf"), height=10, width=18, plot=p3)
ggsave(paste0(opt$outpath, "/GO(Molecular_function).top.", "png"), type="cairo-png", height=10, width=18, plot=p3)
print(paste0(opt$outpath, "/GO(Molecular_function).top.", "png(pdf) is OK"));


#d1<-dp[which(dp[,4]=="Up"),]
#d1_s<-d1[order(d1[,4],decreasing=F),]
#d2<-dp[which(dp[,4]=="Down"),]
#d2_s<-d2[order(d2[,4],decreasing=F),]

#d3_s<-d3[order(d3[,4],decreasing=F),]
#d_l<-rbind(d1,d2)

