path1="LD_intersection.xls"
path2="RD_intersection.xls"


好看的hrbrthemes主题包

d <- read.delim(path1, sep="\t", header=T, quote="")

#conflicted 包可以解决冲突
install.packages("conflicted")
# 或
# install.packages("devtools")
# devtools::install_github("r-lib/conflicted")

#=======================设置ROW.NAMES============================
rownames(DEG)=DEG[,1]
colnames(A)<-c("M","N","O","P")
rownames(B)<-c("a","b","c","d","e")
colnames(B)<-c("m","n","o","p")


#=======================基本数据属性============================
nrow(A) #输出行数
ncol(A) #输出列数
rownames(A) #输出行名
colnames(A) #输出列名

dim(A) #输出矩阵维度
length(A) #输出矩阵的元素个数

#与矩阵类似，但数组的维度可以大于2
myarray <- array(vector, dimensions, dimnames)

dimnames(A)

transform<-t(A)    #t(A)表示矩阵A的转置

as.vector(A)    #将矩阵转化为向量
#=======================改名============================
DEG = plyr::rename(DEG, c("pval"="pValue"))

#========================管道符=========================
#使用管道操作符 直接把数据传递给下一个函数调用
Result <- filter( tidy_gather , Expression>1 ) %>% arrange( Expression )



# ======================数据重塑==================
library("tidyr")

expression <- read.delim("gene_expression.txt",sep="\t",header = T)
#短数据变成长数据
#gather(data=“数据框名”，key=“key名”，value=“value名”，需要转换的列1,2,3)
tidy_gather <- gather(expression,key=Samplename,value = Expression,-id)
#长数据变短数据
tidy_spread <- spread(tidy_gather,key=Samplename,value = Expression)

#按列分割
tidy_separate<- 
  separate(tidy_gather,col=Samplename,into=c("Source","Samplename"),sep="_")

#按列合并
tidy_unite <- 
  unite(tidy_separate,col=Samplename,into=c("Source","Samplename"),sep="_")
#=======================数据清洗=========================

用0替换NA
new_df2[is.na(new_df2)] <- 0
new_df2[is.na(new_df2),] <- 0 #加个逗号就不对了，就是列了

# ========================数据清洗======================
library("dplyr")
cdata["sd"] <- apply(cdata, 1, sd)

去除方差为0的行
cdata2 <- filter(cdata , sd != 0)

去除全是0的行
counts_new=counts[which(rowSums(counts) == 0),]

# ======================删掉某一列================
library(dplyr)
df2<select(df, -Expression )
#增加expression列是原来的Expression的10倍，并替换原来的Expression列
dplyr_mutate_select <- 
  mutate(  tidy_gather, expression = Expression*10  )%>%select( -Expression )

# ======================挑选子集==================
用逻辑判断出位置索引再挑选
d1<-dp[which(dp[,3]=="Up"),]
which(e==8)    #which()函数返回的是位置

根据位置索引提取矩阵
A[1,2]    #访问矩阵的第1行，第2列
A[1,]    #访问矩阵的第1行所有元素
A[,2]    #访问矩阵的第2列所有元素
A[c(1,2),c(2,3,4)]    #访问矩阵中第1,2行、第2,3,4列所有元素

subset函数
two <- subset(two,select=c("gene_id","logFC_A","logFC_B","FDR_A","FDR_B"))

subset(patientdata,age>30 & diabetes=="Type1",c(age,status),drop=F) 

或者用dplyr包
install.packages("dplyr")
library(dplyr)
#筛选表达量大于1的行
dplyr_filter <- filter( tidy_gather, Expression>1)
#筛选表达量大于1且SampleName 为 gene1的行
dplyr_filter <- filter( tidy_gather , Expression >1 & id == "gene1") #或是|

#展示指定的GeneId SampleName  Expression 列
dplyr_select <- select( tidy_separate , id , Samplename , Expression )

#replace the "-/Inf"
向量数据
a<- c(1,2,3,4,5,6,7,8,9)
b<- c(1,3,5,7,9)
a[a %in% b] #提取两个向量的交集
a[!(a %in% b)]#提取两个向量的补集

x <- 1:10
y <- 5:15
union(x,y) #并集
intersect(x, y) #交集
setdiff(x,y)#x和y的差集，以前一个（x）为准
setequal(x,y)#判断x和y是否相等

用数据里的最大值来代替inf
if( "Inf" %in% DEG$log2FoldChange | "-Inf" %in% DEG$log2FoldChange){
  tmp <- DEG[which(DEG$log2FoldChange != "Inf" & DEG$log2FoldChange != "-Inf"),]
  max = max(tmp$log2FoldChange)
  min = min(tmp$log2FoldChange)
  #sub 语句
  DEG$log2FoldChange <- as.numeric(sub("-Inf", min, DEG$log2FoldChange))
  DEG$log2FoldChange <- as.numeric(sub("Inf", max, DEG$log2FoldChange))
}

用其他表格来标记数据
labels <- read.delim(normalizePath(opt$genelist), header=T, sep="\t", check.names=F, quote="")
DEG$label <-''
DEG[match(labels[, 1], DEG[, 1]),]$label = as.character(labels[,1])

标记显著性
two[which(two$FDR_A >= 0.05 | two$FDR_B >= 0.05),'type1'] <- 'no'

取平均值
two["average_F"]=(two["logFC_A"]+two$logFC_B)/2


# ======================排序并且取前几==================
d<-d[order(d["qValue"],decreasing=F),]
d<-head(d, 30)

或者用dplyr包
install.packages("dplyr")
library(dplyr)
##按id进行排序
dplyr_arrange <- arrange(tidy_gather , id ) 
##按id进行排序的基础上按Expression的降序排列
dplyr_arrange1 <- arrange(tidy_gather,id,desc(Expression)) 

# =======================测定某变量在一定范围内的数量==================

install.packages('dplyr')
sum(between(airquality$a,41,42.5))


#======================循环创建向量=======================#
term_list=colnames(d)
term_list_label2<-c()
for (i in 1:length(term_list)) {
  print (paste0("我们现在正在处理",term_list[i]))
  term_list_label2<-c(term_list_label2,term_list[i])
  term_list_label2<-c(term_list_label2," ")
  ####处理表达量文件
  print ("================================================")
}


#=======================
对象类型判断
mode(); class()
is.numeric() #返回值为TRUE或FALSE
is.logical()
is.charactor()
is.data.frame()

对象类型转换
as.numeric() #转换为数值型
as.logical()
as.charactor()
as.matrix()
as.data.frame()

# ======================添加新列==================
install.packages("dplyr")
library(dplyr)
#增加新列ID sub
dplyr_mutate <- mutate(tidy_gather , ID=sub( "gene", "Gene", id ) )

#实际案例 PASTE0
expression <- mutate(expression, Cluster_Name = paste0("cluster", expression$Cluster))

# ======================Groupby----------------
#分组统计平均值
dplyr_groupby_summarise<-group_by( tidy_gather ,Samplename)%>% summarise ( avg=mean( Expression ) )

# ======================计数函数===================

family=as.data.frame(table(expression$Family))


# =======================apply函数=======================
#apply()函数： 函数作用于matrix或者array的margins(可以理解为matrix或者数组的每一行或者每一列)

apply(A,1,sum) #对A矩阵的行，进行求和
apply(A,2,sum)#对A矩阵的列，进行求和

#lapply()函数：遍历列表向量内的每个元素，并且使用指定函数来对其元素进行处理。输入的数据必须是list型对每列进行操作，返回列表向量
#lapply(A, list, FUN……)

apply_test<-function(x){
     myFUN<- function(x, c1, c2) {
         c(sum(x[c1],1), mean(x[c2])) 
       }
     apply(x,1,myFUN,c1='x1',c2=c('x1','x2'))
   }

#for 循环是最慢的

# =======================小实用函数========================
seq()#等间隔函数

seq(from=4, to=10, by=2)
seq(from=0, to=2, length=5)

rep()#重复函数

rep(1,times=5)
rep(c(1,2),times=3)


k <- patientdata 
colnames(k) <- c("A1","A2","A3","A4")
attach(k)  #附变量
lm(A1~A2)
detach(k)

统计类
mean(x)  sd(x)
var(x)   median(x)
quantile(x,p)
cor(x,y) t.test()
lm(y ~ x)
cor.test() -- 有p值也有相关性

# =======================两表相关联=====================

library("dplyr")
#内连接，合并数据仅保留匹配的记录
inner_join(x,y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) 

two <- inner_join(one, zero, by=c("gene_id" = "gene_id"))
#左连接，向数据集x中加入匹配的数据集y记录
left_join(x,y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...)
#右连接，向数据集y中加入匹配的数据集x记录
right_join(x,y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) 
#全连接，合并数据保留所有记录，所有行
full_join(x,y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...)
#返回能够与y表匹配的x表所有记录 
semi_join(x,y, by = NULL, copy = FALSE, ...)
#返回无法与y表匹配的x表的所有记录
anti_join(x, y, by = NULL, copy = FALSE, ...)  


cbind()    #按列合并
cbind(A,B)

rbind()    #按行合并
rbind(A,B)

# =======================堆叠图必备-stack=================

d_l <- rbind(d,d2_new)

uniq1<-d_l$Term[!duplicated(d_l$Term)] 
d_l$Term <- factor(d_l$Term, levels=uniq1) #设置factor

因子（factor)是R语言中比较特殊的一个数据类型， 它是一个用于存储类别的类型。



# =======================批量出图必备=================
读取文件名并且拆分
groupname <- gsub("\\.(txt|xls)$", "", gsub("-all.gene.xls$", "", path1))#opt$input
Case_names <- unlist(strsplit(groupname,"-vs-"))[1]
Control_names <- unlist(strsplit(groupname,"-vs-"))[2]

批量改名
# 用file.rename函数修改
file.rename(filename1,filename2)


option_list = list(
  make_option( c("-i", "--input" ), type = "character",
               help = "The input differential gene file(force).  e.g. *-vs-*-all.gene.xls" ),
  make_option( c("-p", "--pval" ), type = "double",default = 0.05,
               help = "pValue ratio threshold, default: 0.05 ."),
  make_option( c("-f", "--foldchange" ), type = "double",default = 2,
               help = "foldchange threshold, default: 2 ."),
  make_option( c("-l", "--genelist" ), type = "character" ,
               help = "Genelist to display the gene symbol.  e.g. genelist.xls"),
  make_option(c("-t", "--title"), type = "character", default = NULL,
              help = "Graphic title and outputfile information: Group_A-vs-Group_B . "),
  make_option( c("-o", "--outputdir" ),type="character", default = "./Volcano",
               help="the output directory of Volcano results, default: ./Volcano ." )
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


if(is.null(opt$input) | is.null(opt$outpath)){
  print_help(opt_parser)
  stop("--input --outpath --mark must be supplied", call.=FALSE)
}
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)}

opt$outpath=paste0(opt$outpath,'/',groupname,"/")
if(!file.exists(opt$outpath)){dir.create(opt$outpath,recursive = T)} #grdevice 必须要找到路径才能保存

ggsave(paste0(opt$outpath,"/", groupname, ".pdf"), height=10, width=10, plot=p)
ggsave(paste0(opt$outpath,'/', groupname, ".png"), type="cairo-png", height=10, width=10, plot=p)
print(paste0(opt$outpath,'/', groupname, ".png(pdf) is OK"))

write.table(datat,file=paste0(output_dir ,"/", picname,".genes_reorder.cluster_result.xls"),row.names=FALSE,quote = FALSE,sep='\t')

# =======================常用包=================
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(optparse))

  biomaRt ID转换