#Author: Chen congjia
#WGCNA 批量修改散点图+线性拟合+相关性检测

library("optparse")
library('ggplot2')
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default=NULL, help="outfile directory", metavar="character")
);

opt_parser = OptionParser(option_list=option_list,epilogue = "Rscript top***.r -i enrichment-go-Group1-vs-Group2-Down.txt -m Down -o outdir/");
opt = parse_args(opt_parser);
if(is.null(opt$input) | is.null(opt$outpath)){
  print_help(opt_parser)
  stop("--input --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}

###################识别名字长短###############################
groupname <- gsub("\\.(xls|csv)$", "", basename(opt$input))
print (groupname)
name <-strsplit(groupname,"_")
if(lengths(name)==2){
  name_x <- unlist(strsplit(groupname,"_"))[2] #颜色
  print (name_x)  
  name_y <- unlist(strsplit(groupname,"_"))[1] #性状
  print (name_y)
}
if(lengths(name)==3){
  name_z <- unlist(strsplit(groupname,"_"))[2]
  name_x <- unlist(strsplit(groupname,"_"))[3] #颜色
  print (name_x)  
  name_y <- unlist(strsplit(groupname,"_"))[1] #性状
  name_y <- paste0(name_y,"_",name_z)
  print (name_y)
}
if(lengths(name)==4){
  name_o <- unlist(strsplit(groupname,"_"))[3]
  name_z <- unlist(strsplit(groupname,"_"))[2]
  name_x <- unlist(strsplit(groupname,"_"))[4] #颜色
  print (name_x)  
  name_y <- unlist(strsplit(groupname,"_"))[1] #性状
  name_y <- paste0(name_y,"_",name_z,"_",name_o)
  print (name_y)
}
#######读取文件名字，并且用在图中和读取列中##############

############重命名列###########################
the.data <- read.csv(opt$input, header = TRUE)
colnames(the.data)[3] <- "Moddule.Membership"
colnames(the.data)[4] <- "Gene_correlation"

############拟合曲线#############################
model.lm<-lm(formula = Gene_correlation ~ Moddule.Membership, data = the.data)
summary(model.lm)

"spearman"
cor1 <- cor(the.data$Moddule.Membership,the.data$Gene_correlation,method = "pearson")

l <- list(a = as.numeric(format(coef(model.lm)[1], digits = 4)),
          b = as.numeric(cor1),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
####################统计结束#############################
eq <- substitute(italic(Cor)~"="~b~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)

#italic(y) == a + b %.% italic(x)~","~

p=ggplot(data=the.data, aes(x=Moddule.Membership, y=Gene_correlation))+
  geom_point(color=name_x)+
  stat_smooth(method="lm",formula = y~x,se=FALSE)+
  labs(x=paste0("Module Membership in ",name_x," module"), y=paste0("Gene significance for ",name_y), title = paste0("ModuleMembership-vs-GeneSignificance"))+
  theme(plot.title = element_text(hjust = 0.5, vjust=2, size=15, family = "ArialMT"))+
  geom_text(aes(x = 0.7,y=0, label = as.character(as.expression(eq))), parse = TRUE)

opt$outpath=paste0(opt$outpath,'/',name_y,"/")
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)} #grdevice 必须要找到路径才能保存

ggsave(paste0(opt$outpath,"/", groupname, ".pdf"), height=10, width=10, plot=p)
ggsave(paste0(opt$outpath,'/', groupname, ".png"), type="cairo-png", height=10, width=10, plot=p)
print(paste0(opt$outpath,'/', groupname, ".png(pdf) is OK"))


