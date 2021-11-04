#!/usr/bin/env Rscript
#用来输出堆叠的同时具有up和down 的图
library("optparse")
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
  make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character"),
  make_option(c("-n", "--Number"), type="character", default=10, help="The number of terms you want to conserve", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top10X3_GO.r -i enrichment-kegg-Group1-vs-Group2-Down.txt -m Down -o outdir/  Rscript top10X3_kegg.r -i enrichment-kegg-H-vs-C-Fake.xls -m Total -o outdir/
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
d2<-aggregate(d[11], by=list(Term), FUN=sum)
#重新命名
names(d2)<-c('Term','Enrichment_score_sum') 

print (d2[2])
#把这个数据重新赋值到原来的表下面
#与之前的数据合并
d3=merge(d,d2,by.x="Term",by.y="Term",all=T)


##########筛选出top多少的term名称###############
#定义函数
top10 <- function(i) { return(i[order(head(i, Number)["Enrichment_score_sum"],decreasing=T),]) }
#先筛选和排序
#d3 <- d3[which(d3[,"ListHits"]>2),]
d3 <- d3[order(d3[,"Enrichment_score_sum"],decreasing=T),]
#除了p值外，其他都是升序
#stopifnot(nrow(d)>0)
if(nrow(d)==0){
  print("d items = 0, program exit!")
  q()
}

#创建Up,down 各top10
d_l <- rbind(top10(d3[d3["Regulation"]=="Up", ]),
            top10(d3[d3["Regulation"]=="Down", ]))

#写入table
write.table(d_l[,c(1,2,3,4,9,10,11,13)], paste0(opt$outpath, "/kegg.top.", opt$mark, ".xls"), sep="\t", quote=FALSE,
            col.names=TRUE, row.names=FALSE)

#d1<-dp[which(dp[,4]=="Up"),]
#d1_s<-d1[order(d1[,4],decreasing=F),]
#d2<-dp[which(dp[,4]=="Down"),]
#d2_s<-d2[order(d2[,4],decreasing=F),]

#d3_s<-d3[order(d3[,4],decreasing=F),]
#d_l<-rbind(d1,d2)

print ("###############check################")
print (d_l)
d_l["Term"] <- apply(d_l["Term"], 1, limit70)
print (d_l$Term)
uniq<-d_l$Term[!duplicated(d_l$Term)]
print ("uniq")
print (uniq) 
d_l$Term <- factor(d_l$Term, levels=uniq)#factor函数会输出level代表不重复的次数
print ('Unduplicated')
print (d_l$Term)
print ("################End###################")
#为了获得Up在左，DOwn在右的效果
d_l$Regulation <- factor(d_l$Regulation, levels=c("Down", "Up"));

p=ggplot(data=d_l, aes(x=Term, y=Enrichment_score, width=0.8, fill=Regulation,space=0)) +
  coord_flip()+
  geom_bar(stat="identity",position="stack") +
  scale_x_discrete(breaks=d_l$Term, labels=d_l$Term)+ 
  guides(fill=guide_legend(reverse=TRUE))+
  scale_fill_brewer(palette = "Set1")+
  labs(x="", y='Enrichment_score', title = paste0(groupname,": ","Top ",opt$Number, " Up & Down KEGG Term"))+  
  theme_bw() + scale_colour_gradient(low = "green", high = "red") +
  theme(axis.text.y=element_text(size=11,color="black")) + 
  theme(axis.text.x=element_text(hjust=1, size=14,color=d_l$color))+
  theme(legend.position=c(-0.10,0.8), legend.key.width=unit(1, "lines"))+
  theme(legend.position="right")
theme(legend.text=element_text(size=11))+
  theme(plot.title = element_text(hjust = 0.5, vjust=4, size=18, family = "ArialMT"))+
  theme(plot.margin=unit(c(2,4,2,18), "lines"))+
  theme(panel.grid =element_blank()) 

ggsave(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".pdf"), height=10, width=18, plot=p)
ggsave(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png"), type="cairo-png", height=10, width=18, plot=p)
print(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png(pdf) is OK"));