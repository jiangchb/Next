library(ggplot2)
args <- commandArgs(T)

df <- read.delim(args[1],sep="\t",header = T)
df2 <- df[,c(1,4,6)]

df2[which(df2$regulator=="Up"), "pval"] <- -log10(df2[which(df2$regulator=="Up"), "pval"])
df2[which(df2$regulator=="Down"), "pval"] <- log10(df2[which(df2$regulator=="Down"), "pval"])
df2$just = ifelse( df2$pval<0,0,1)

p= ggplot(df2,aes(x=reorder(term,pval),y=pval))+
  geom_bar(stat = "identity",width = 0.5,aes(fill=regulator))+
  scale_fill_manual(values = c("#255AA8", "#EE422A")) +
  #scale_x_discrete(labels = labelname)+
  theme_bw()+
  xlab("")+labs(title = "KEGG oe-vs-con")+
  theme(plot.title = element_text(hjust = 0.5, size = 20))+
  geom_text(aes(x= term, y=0, label = term), size = 2.5,hjust=df2$just*1.01) +
  theme(axis.text.y=element_blank()) +theme(axis.ticks=element_blank()) +
  coord_flip()

#p = ggplot(df2,aes(x=reorder(term,pval),y=pval))+
#    geom_bar(stat = "identity",width = 0.5,aes(fill=regulator))+
#    scale_fill_manual(values = c("#255AA8", "#EE422A")) +
#    #scale_x_discrete(labels = labelname)+
#    theme_bw()+
##    labs(title = "KEGG oe-vs-con",y="-log10 pvalue")+
#	labs(title = paste0("KEGG ",args[2])) + 
#	xlab("") +
#    theme(plot.title = element_text(hjust = 0.3, size = 20))+
#    theme(legend.position=c(.9,.08)) +
#	theme(axis.text.y=element_text(size=10))+
# # geom_text(aes(x= term, y=0, label = term), size = 2.5,hjust=df2$just*1.01) +
#    coord_flip()
ggsave(paste(args[2],".pdf",sep=""),height=12,width=10,plot=p)
