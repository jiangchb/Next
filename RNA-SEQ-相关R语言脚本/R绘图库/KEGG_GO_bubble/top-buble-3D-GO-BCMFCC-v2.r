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

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-go-", "", basename(opt$input)))
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

top30 <- function(i) { return(i[order(head(i, 30)["pValue"],decreasing=T),]) }

d <- read.delim(opt$input, sep="\t", header=T, quote="")
d <- d[which(d[,"ListHits"]>1),]
d <- d[order(d[,"pValue"]),]

#stopifnot(nrow(d)>0)
if(nrow(d)==0){
    print("d items = 0, program exit!")
    q()
}

dp <- rbind(top30(d[d["Category"]=="biological_process", ]),
	top30(d[d["Category"]=="cellular_component", ]),
	top30(d[d["Category"]=="molecular_function", ]))
write.table(dp[,c(1,2,3,4,8,10,11)], paste0(opt$outpath, "/GO.top30.", opt$mark, ".xls"), sep="\t", quote=FALSE,
	col.names=TRUE, row.names=FALSE)

d_l=top30(d)

d_l["Term"] <- apply(d_l["Term"], 1, limit70)
d_l$Term <- factor(d_l$Term, levels=d_l$Term)
p=ggplot(d_l,aes(x=Term,y=Enrichment_score,shape=Category))+
  geom_point(aes(color=pValue,size=ListHits))+
  scale_colour_gradientn(colours=rainbow(6)) + #need self-define
  scale_size_area(max_size = 10) +#adjust the size of the bubble
  theme(legend.justification=c(0,0),legend.position=c(1,1))+
  coord_flip()+
  labs(color="pValue",x="", y="Enrichment_score", fill="pValue", title = paste0(groupname,": ","Top 30 of GO enrichment")) + 
 
  theme_bw() 

ggsave(paste0(opt$outpath, "/GO.top.", opt$mark, ".pdf"), height=10, width=10, plot=p)
ggsave(paste0(opt$outpath, "/GO.top.", opt$mark, ".png"), type="cairo-png", height=10, width=10, plot=p)
print(paste0(opt$outpath, "/GO.top.", opt$mark, ".png(pdf) is OK"));
