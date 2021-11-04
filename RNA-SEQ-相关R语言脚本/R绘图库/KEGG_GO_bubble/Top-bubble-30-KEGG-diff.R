path1="LD_intersection.xls"
path2="RD_intersection.xls"

d <- read.delim(path1, sep="\t", header=T, quote="")


d<-d[order(d["qValue"],decreasing=F),]
d<-head(d, 30)

d<-d[order(d["qValue"],decreasing=T),]
d2 <- read.delim(path2, sep="\t", header=T, quote="")

d2_new<-semi_join(d2,d,by="Term")

d2_new["Group"] <- "Reversine-vs-DMSO"
d["Group"] <- "LavendustinA-vs-DMSO"

d_l <- rbind(d,d2_new)

uniq1<-d_l$Term[!duplicated(d_l$Term)] 
d_l$Term <- factor(d_l$Term, levels=uniq1)

p=ggplot(d_l,aes(x=Term,y=Enrichment_score,shape=Group))+
  geom_point(aes(color=qValue,size=ListHits))+
  scale_colour_gradientn(colours=cividis(6)) + #need self-define
  scale_size_area(max_size = 7) +#adjust the size of the bubble
  theme(legend.justification=c(0,0),legend.position=c(1,1))+
  coord_flip()+
  theme(text=element_text(size=14,family="serif")) +
  labs(color="qValue",x="", y="Enrichment_score", fill="qValue", title = paste0("Top 30 of KEGG enrichment")) + 
  
  theme_bw() 
p
ggsave(paste0( "./KEGG.top.","pdf"), height=10, width=10, plot=p)
ggsave(paste0( "./KEGG.top.", "png"), type="cairo-png", height=10, width=10, plot=p)
print(paste0("./KEGG.top.",  "png(pdf) is OK"));
write.table(d_l[,c(1,2,3,4,5,6,7,8,9)], paste0("./KEGG.top60.xls"), sep="\t", quote=FALSE,
            col.names=TRUE, row.names=FALSE)