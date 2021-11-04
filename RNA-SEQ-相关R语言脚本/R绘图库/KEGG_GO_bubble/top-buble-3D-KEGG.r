#!/usr/bin/env Rscript
#三维输出R语言图
library("optparse")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-m", "--mark"), type="character", default=NULL, help="select Up, Total or Down", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top10X3_GO.r -i enrichment-go-Group1-vs-Group2-Down.txt -m Down -o outdir/");
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

limit70 <- function(s) {
	k <- as.character(s)
	if(str_length(s)>70){k <- sub("[^ ]+$", "...", substr(k,1,67))}
	return(k)	
}

top20 <- function(i) { return(i[order(head(i, 20)["expression"],decreasing=F),]) }

d <- read.delim(opt$input, sep="\t", header=T, quote="")
#d <- d[which(d[,"ListHits"]>1),]
d <- d[order(d[,"pValue"]),]
d["expression"] <- -log(d$pValue,10)

#stopifnot(nrow(d)>0)
if(nrow(d)==0){
    print("d items = 0, program exit!")
    q()
}
write.table(d[,c(1,2,3,4,8,10,11,12)], paste0(opt$outpath, "/KEGG.top20.", opt$mark, ".xls"), sep="\t", quote=FALSE,
	col.names=TRUE, row.names=FALSE)

d_l=top20(d)

d_l["Term"] <- apply(d_l["Term"], 1, limit70)
d_l$Term <- factor(d_l$Term, levels=d_l$Term)
ylabs=expression('-log'[10]*' Pvalue')
p=ggplot(d_l,aes(x=Term,y=expression))+
  geom_point(aes(color=pValue,size=ListHits))+
  scale_colour_gradientn(colours=rainbow(6)) + #need self-define
  scale_size_area(max_size = 10) + #adjust the size of the bubble
  theme(legend.justification=c(0,0),legend.position=c(1,1)) +
  coord_flip()+
  labs(x="", y=ylabs, title = paste0("Top 20 ",opt$mark, " kegg enrichment"),font.axis=2,font.lab=2,cex.axis=2) + 
  theme(axis.text.x = element_text(face="bold", color="blue", size=8))+  
  theme(axis.title.x = element_text(size = 20))+
  theme_bw() 

ggsave(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".pdf"), height=10, width=7, plot=p)
ggsave(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png"), type="cairo-png", height=10, width=7, plot=p)
print(paste0(opt$outpath, "/KEGG.top.", opt$mark, ".png(pdf) is OK"));
