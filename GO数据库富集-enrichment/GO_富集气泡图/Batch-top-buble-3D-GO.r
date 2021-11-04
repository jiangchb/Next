#!/usr/bin/env Rscript
#三维输出R语言图
library("optparse")
library("viridis")
library("dplyr")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character"),
	make_option(c("-n", "--Number"), type="character", default=10, help="The number of terms you want to conserve", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top10X3_GO.r -i enrichment-go-Group1-vs-Group2-Down.txt -m Down -o outdir/");
opt = parse_args(opt_parser);
if(is.null(opt$input) | is.null(opt$outpath)){
	print_help(opt_parser)
	stop("--input --outpath --mark must be supplied", call.=FALSE)
}
opt$outpath<-gsub("/$", "", opt$outpath)

Number <-as.numeric(opt$Number)
print ("We will do top ")
print (Number)

library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)


groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-kegg-", "", basename(opt$input)))
if(grepl("-Down$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Down$", "(Down)", groupname)
		mark <- "Down"
}
if(grepl("-Up$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Up$", "(Up)", groupname)
		mark <- "Up"
}
if(grepl("-Total$", groupname, ignore.case = F, perl = F,fixed = F, useBytes = F)){
        groupname <- gsub("-Total$", "(Total)", groupname)
		mark <- "Total"
}

limit70 <- function(s) {
	k <- as.character(s)
	if(str_length(s)>70){k <- sub("[^ ]+$", "...", substr(k,1,67))}
	return(k)	
}
topn <- function(i) { return(i[order(head(i, Number)["qValue"],decreasing=T),]) }

d <- read.delim(opt$input, sep="\t", header=T, quote="")
d <- d[which(d[,"ListHits"]>1),]
d <- filter( d , qValue != 0 ) 
d <- d[order(d[,"qValue"]),]

#d["expression"] <- -log(d$pValue,10)
d["GeneRatio"] <- d$ListHits/d$ListTotal

#stopifnot(nrow(d)>0)
if(nrow(d)==0){
    print("d items = 0, program exit!")
    q()
}

d_l=topn(d)

#d_l<-d_l[order(d_l["GeneRatio"],decreasing=F),]

d_l["Term"] <- apply(d_l["Term"], 1, limit70)

d_l$Term <- factor(d_l$Term, levels=d_l$Term)
print (d_l)



ylabs="Enrichment score"
size.lab <- "Gene Counts"

palette <- colorRampPalette(c("firebrick3","white","navy"))(n=256)

p=ggplot(d_l,aes(x=Term,y=Enrichment_score))+
  geom_point(aes(color=qValue,size=ListHits))+
  scale_colour_gradientn(colours=palette) + #need self-define
  scale_size_area(max_size = 10) + #adjust the size of the bubble
  theme(legend.justification=c(0,0),legend.position=c(1,1)) +
  coord_flip()+
  labs(x="", y=ylabs, title = paste0("Top ",Number," ",opt$mark, " GO enrichment"),size=size.lab,font.axis=2,font.lab=2,cex.axis=2) + 
  theme(axis.text.x = element_text(face="bold", color="blue", size=8))+  
  theme(axis.title.x = element_text(size = 20))+
  theme_bw() 
mark="enrichment"
opt$outpath <- paste0(opt$outpath,"/",groupname,"/")
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
write.table(d_l[,c(1,2,3,4,5,6,7,8,9,10,11)], paste0(opt$outpath, "/GO.qValue_order.", mark, ".xls"), sep="\t", quote=FALSE,
	col.names=TRUE, row.names=FALSE)
ggsave(paste0(opt$outpath, "/GO.top.", mark, ".pdf"), height=8, width=8, plot=p)
ggsave(paste0(opt$outpath, "/GO.top.", mark, ".png"), type="cairo-png", height=8, width=8, plot=p)
print(paste0(opt$outpath, "/GO.top.", mark, ".png(pdf) is OK"));
