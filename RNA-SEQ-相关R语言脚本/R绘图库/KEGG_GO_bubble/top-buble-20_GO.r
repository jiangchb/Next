#!/usr/bin/env Rscript
library("optparse")
library("oebio")
option_list = list(
	make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
	make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top20_GO.r -i enrichment-go-Group1-vs-Group2-Down.txt -o outdir/");
opt = parse_args(opt_parser);
if (is.null(opt$input) | is.null(opt$outpath) ){
	print_help(opt_parser)
	stop("--input --outpath must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}
opt$outpath<-gsub("/$", "", opt$outpath)

library(ggplot2);

enrich <- read.delim(opt$input, head=TRUE, sep="\t", quote="")
enrich <- enrich[which(enrich["ListHits"]>2), ]

groupname <- gsub("\\.(txt|xls)$", "", gsub("^enrichment-go-", "", basename(opt$input)))

top20 <- head(enrich[order(enrich[,"pValue"]),], 20)
#stopifnot(nrow(top20)>0)
if(nrow(top20)==0){
    print("top20 items = 0,program exit!")
    q()
}

write.table(top20[,c(1,2,3,7,9,10)], paste0(opt$outpath, "/GO.top", ".xls"), sep="\t",
	quote=FALSE, col.names=TRUE, row.names=FALSE)
top20[, "term"] <- sub("^path:", "", paste(top20[,1], top20[,2], sep=": "))

xlab <- "Enrichment Score"
title <- paste0(groupname, ": ", "GO Enrichment top 20")
size.lab <- "Number"
p=qplot(Enrichment_score, term, data=top20, size=ListHits, 
	colour=pValue, xlab=xlab, ylab="")+ggtitle(title)+
	theme_bw()+
	theme(plot.title = element_text(hjust = 0.5))+
	scale_colour_gradientn(colours=rainbow(6)) +
	labs(size=size.lab) + labs(colour="pValue")
ggsave(paste0(opt$outpath, "/GO.top", ".pdf"), height=8, width=10, plot=p)
ggsave(paste0(opt$outpath, "/GO.top", ".png"), type="cairo-png", height=8, width=10, plot=p)
print(paste0(opt$outpath, "/GO.top", ".png(pdf) is OK"));
