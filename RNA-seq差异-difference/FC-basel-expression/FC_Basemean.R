根据DESeq跑出来的差异文件来进行绘制，为方便看点的离散程度，进行了log10 变化

library(ggplot2)
ggplot(my_df,aes(x=fold_change,y=WT_baseMean,color=color))+
  geom_point(aes(color=my_df$theme)) +
  scale_x_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000))+
  scale_y_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000))+
  scale_color_manual(values = c("blue", "red"))+theme(legend.position = "right")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black")) +
  coord_flip()