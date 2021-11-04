library("optparse")
library('ggplot2')
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file name", metavar="character"),
  make_option(c("-c", "--corinput"), type="character", default=NULL, help="input cor file name", metavar="character"),
  make_option(c("-p", "--pvalinput"), type="character", default="./moduleTraitCor.csv", help="input pval file name", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="./moduleTraitPvalue.csv", help="outfile directory", metavar="character")

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


###d读取两张表####
"pearson"
cor_data <- read.csv(normalizePath(opt$corinput), header = TRUE)
pval_data <- read.csv(normalizePath(opt$pvalinput), header = TRUE)
"cor:"
print (cor_data[which(cor_data[1]==name_x),name_y])
cor=cor_data[which(cor_data[1]==name_x),name_y]
"p-value:"
print (pval_data[which(pval_data[1]==name_x),name_y])
pval=pval_data[which(pval_data[1]==name_x),name_y]
############拟合曲线#############################
#model.lm<-lm(formula = Gene_correlation ~ Moddule.Membership, data = the.data)
#summary(model.lm)
##################################这个值并非相关性值###################################
#cor1 <- cor(the.data$Moddule.Membership,the.data$Gene_correlation,method = "pearson")

#l <- list(a = as.numeric(format(coef(model.lm)[1], digits = 4)),
#          b = as.numeric(cor1),
#          r2 = format(summary(model.lm)$r.squared, digits = 4),
#          p = format(summary(model.lm)$coefficients[2,4], digits = 4))

#eq <- substitute(italic(Cor)~"="~b~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)

#eq2 <- substitute(Cor=b,R^2=r2,P=p, l)
#italic(y) == a + b %.% italic(x)~","~
####################统计结束#############################
p=ggplot(data=the.data, aes(x=Moddule.Membership, y=Gene_correlation))+
  geom_point(color=name_x)+
  stat_smooth(method="lm",formula = y~x,se=FALSE)+
  labs(x=paste0("Module Membership in ",name_x," module"), y=paste0("Gene significance for ",name_y), 
  title =paste0("ModuleMembership-vs-GeneSignificance","\n","cor:",as.character(format(cor,digits=4)),", p value:",as.character(format(pval, digits = 4))))+
  theme(plot.title = element_text(hjust = 0.5, vjust=2, size=15, family = "ArialMT"))
  #annotate(geom = 'text', label = as.character(as.expression(eq)), x = -Inf, y = Inf,hjust=1,vjust=0) #用来固定位置 hjust vjust

opt$outpath=paste0(opt$outpath,'/',name_y,"/")
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)} #grdevice 必须要找到路径才能保存
ggsave(paste0(opt$outpath,"/", groupname, ".pdf"), height=10, width=10, plot=p)
ggsave(paste0(opt$outpath,'/', groupname, ".png"), type="cairo-png", height=10, width=10, plot=p)
print(paste0(opt$outpath,'/', groupname, ".png(pdf) is OK"))
file.copy(opt$input, opt$outpath)

