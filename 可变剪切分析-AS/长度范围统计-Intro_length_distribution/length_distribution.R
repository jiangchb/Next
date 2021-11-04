path1="intron.gtf"

library("dplyr")

d <- read.delim(path1, sep="\t", header=F, quote="")

d1<-d[which(d[,3]=="exon"),]

d2 <- subset(d1,select=c("V1","V4","V5"))

d2["length"] <- (d2["V5"]-d2["V4"])+1 #LENGTH 就是终止-起始+1

#区间计算1
d2[["interval"]] <- findInterval(d2[["length"]], c(0, 100, 300))
#findInterval这个函数是用来判断在哪个区间内的，区间又后面的向量决定，可行，但略烦


#筛选掉2000以上的
d3 <- filter( as.data.frame(d2) , length <= 2000)

#使用cut函数切割，断点为0到2000，间距50，table函数进行计数
options(scipen = 1)
data <- as.data.frame(table(cut(d2$length, breaks = seq(0,2000,50))))

write.table(data,file=paste0("tmp_result.xls"),row.names=F,quote = FALSE,sep='\t')
#画图
path2="tmp_result.xls"
df <- read.delim(path2, sep="\t", header=T, quote="")

library(scales)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
df$Var1=factor(df$Var1,levels=c(as.character(df$Var1)))
p=ggplot(df,aes(Var1,Freq))+
  geom_bar(stat="identity",fill="Royalblue",width=0.6)+
  labs(x="Length range",y="Numbers",title = "Intro length distribution")+
  theme_economist_white()+
  theme(axis.text.x = element_text(angle = 80, size=5, color="black",vjust = 0.9))+
#  geom_text(aes(label=Freq),hjust=0.5, vjust=-0.5,size=2.5)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_text(vjust = 0.5))+
  theme(title=element_text(size = 8, vjust=2 , angle = 00))+
  theme(panel.grid =element_blank())
p

ggsave("length_distribution.pdf",width=5,height=3.7,plot=p)
ggsave("length_distribution.png",type="cairo-png",width=5,height=3.7,plot=p)
