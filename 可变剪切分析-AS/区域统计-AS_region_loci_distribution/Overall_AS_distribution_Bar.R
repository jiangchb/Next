path1="result.xls"

library(ggplot2)
library(stringr)
library(grid)
library(RColorBrewer)
library(ggthemes)

d <- read.delim(path1, sep="\t", header=T, quote="")

d$Label <- factor(d$Label, levels=d$Label)

p=ggplot(data=d,aes(x=Label, y=counts)) +
  geom_bar(stat="identity",position=position_dodge(0.7),width=0.5,fill="#BB0021FF") +
  labs(x="", y='Number',title = "AS region distribution")+geom_text(aes(label=counts),hjust=0.5, vjust=-0.5,size=3.5)+
  theme_classic()+  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_text(vjust = 0.5))+
  theme(title=element_text(size = 8, vjust=2 , angle = 00))+
  theme(panel.grid =element_blank())
p
ggsave("AS_region_distribution.pdf",width=5,height=4.2,plot=p)
ggsave("AS_region_distribution.png",type="cairo-png",width=5,height=4.2,plot=p)