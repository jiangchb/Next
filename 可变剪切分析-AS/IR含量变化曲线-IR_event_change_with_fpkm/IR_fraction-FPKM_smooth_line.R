path1="LINC00278_AS_stat_forplot.xls"
path2="PCDNA3_1_AS_stat_forplot.xls"

df1<- read.table(path1, header = T, sep = "\t", quote = "", check.names = F)
df1 = plyr::rename(df1, c("LINC00278"="FPKM"))

df2 <- read.table(path2, header = T, sep = "\t", quote = "", check.names = F)
df2 = plyr::rename(df2, c("PCDNA3_1"="FPKM"))


df1["Sample"]="LINC00278"
df2["Sample"]="PCDNA3_1"

df_bind=rbind(df1,df2)

#求p值

library(ggplot2)
library(ggthemes)
library(ggsci)
#绘图
p <- ggplot(df_bind, aes(x=-log10(FPKM), y=IR_Fraction,color=Sample,group=Sample)) + 
#  scale_color_manual(values = color1,
#             labels = c(as.character(unique(df_bind$Sample))))+
  scale_color_aaas()+
  stat_smooth(method="auto", se=FALSE)  + theme_stata() + theme(legend.position=c(0.2,0.85))

p

outPDFfile <- paste( "./","LINC00278-vs-PCDNA3_1", '.AS_IR_FRACTION.pdf', sep = '')
outPNGfile <- paste( "./","LINC00278-vs-PCDNA3_1", '.AS_IR_FRACTION.png', sep = '')

ggsave(outPDFfile, width = 5, height = 5, plot = p)
ggsave(outPNGfile, type = "cairo-png", width = 5, height = 5, plot = p)


