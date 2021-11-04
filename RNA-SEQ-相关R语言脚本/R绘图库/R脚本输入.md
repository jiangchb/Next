

# R脚本输入

## Page 10：程序包的安装及加载

```R
#Bioconductor包
# 3.5版本以下
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/",
             repos=structure(c(CRAN="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")))
source("http://www.bioconductor.org/biocLite.R") 
biocLite("")
# 3.5版本以上
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
options(BioC_mirror=‘https://mirrors.tuna.tsinghua.edu.cn/bioconductor’) 
BiocManager::install(“ROC", version = "3.8") 
```

## Page 11：程序包的安装及加载

```R
local({
  r = getOption('repos')
  options(repos = r, BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor')
})
```

## Page 12：程序包的安装及加载

```R
# 设置安装镜像
options(repos=structure(c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) 
#程序包的安装
install.packages (“tidyr”)
#程序包的加载
library (“tidyr”)
```

## Page 13：程序包的安装及加载

```R
if(!require(devtools))
  install.packages("devtools")
if(!require(ggcor))
  devtools::install_github("houyunhuang/ggcor")
```

## Page 14：程序包的安装及加载

```R
# 如果R环境是基于conda构建，可以采用conda来安装
conda install -c bioconda r-pheatmap
conda install -c bioconda bioconductor-ComplexHeatmap
```

## Page 22：设置工作目录

```R
getwd()
setwd("D:/MyRProject")
getwd()
dir()
```

---

## Page 29：R数据结构-（1）向量

```R
a<-c(1,2,3,4,5,6,7,8,9,10)    #等价于  
a<-c(1:10)

b<-c("one", "two", "three") 
c<-c("TRUE", "FALSE", "TRUE", "FALSE")
d<-c(1,2,3, "three", "hello")

e<-c(1,7,8,2,5,6)
e[3]    #第三个值
e[-2]    #除第二个值以外的其他值
e[-2:-4]
e[c(2,5,3)]
e[e<3]
which(e==8)    #which()函数返回的是位置
```

---

## Page 30：R数据结构-（1）向量

```R
set.seed(250)
a = runif(3, min=0, max=100)
floor(a)
ceiling(a)
round(a,4)
sum(a)
mean(a)
max(a)
min(a)
sort(a)

b<-c("1","a","2","b")
sort(b)

seq()#等间隔函数

seq(from=4, to=10, by=2)
seq(from=0, to=2, length=5)

rep()#重复函数

rep(1,times=5)
rep(c(1,2),times=3)
```

---

## Page 31：R数据结构-（1）向量

```R
a<- c(1,2,3,4,5,6,7,8,9)
b<- c(1,3,5,7,9)
a %in% b
!(a %in% b)
a[a %in% b] #提取两个向量的交集
a[!(a %in% b)]#提取两个向量的补集
```

## Page 32：R数据结构-（1）向量

```R
x <- 1:10
y <- 5:15
union(x,y) #并集
intersect(x, y) #交集
setdiff(x,y)#x和y的差集，以前一个（x）为准
setequal(x,y)#判断x和y是否相等

```



##  Page 33：R数据结构-（2）矩阵

```R
A <- matrix(1:20, nrow=5, ncol=4, byrow=T)
B <- matrix(1:20, nrow=5, ncol=4, byrow=F)
C <- matrix(1:20, nrow=5, ncol=4) # 默认byrow=F
```

---

##  Page 34： R数据结构-（2）矩阵

```R
A[1,2]    #访问矩阵的第1行，第2列
A[1,]    #访问矩阵的第1行所有元素
A[,2]    #访问矩阵的第2列所有元素
A[c(1,2),c(2,3,4)]    #访问矩阵中第1,2行、第2,3,4列所有元素

rownames(A)<-c("A","B","C","D","E")
colnames(A)<-c("M","N","O","P")
rownames(B)<-c("a","b","c","d","e")
colnames(B)<-c("m","n","o","p")

nrow(A)
ncol(A)
rownames(A)
colnames(A)
dim(A)
length(A)
dimnames(A)
```

---

## Page 35：R数据结构-（2）矩阵

```R
cbind()    #按列合并
cbind(A,B)

rbind()    #按行合并
rbind(A,B)

transform<-t(A)    #t(A)表示矩阵A的转置

as.vector(A)    #将矩阵转化为向量
```

---

## Page 36：R数据结构-（3）数据框

```R
patientID <- c(1, 2, 3, 4)
age <- c(25, 34, 28, 52)
diabetes <- c("Type1", "Type2", "Type1", "Type1")
status <- c("Poor", "Improved", "Excellent", "Poor")
patientdata <- data.frame(patientID, age, diabetes, status)
```

---

## Page 37：R数据结构-（3）数据框

```R
patientdata[1,3] 
patientdata[1,2]
patientdata[1:2]

patientdata[c("diabetes", "status")]
patientdata$age
patientdata[2:3,-2]
```

---

## Page 38：R数据结构-（3）数据框

```R
subset(patientdata,age > 30, c(age,status),drop=F)
subset(patientdata,age > 30, c(-age,-status),drop=F)
subset(patientdata,age>30 & diabetes=="Type1",c(age,status),drop=F) 
A
subset(A,A[,2]>6,select=c(1,3))
```

---

## Page 39：R数据结构-（3）数据框

```R
k <- patientdata 
colnames(k) <- c("A1","A2","A3","A4") 
```

---

## Page 40：R数据结构-（3）数据框

```R
lm(A1~A2)

attach(k) 

lm(A1~A2) 

detach(k) 

lm(A1~A2) 
```

---

## Page 41：R数据结构-（4）数组

```R
dim1 <- c("A1", "A2")
dim2 <- c("B1", "B2", "B3")
dim3 <- c("C1", "C2", "C3","C4")
a <- array(1:24, c(2, 3, 4), dimnames=list(dim1, dim2, dim3))

a[2,1,2]

a[1,2:3,2:3]
```

---

## Page 42：R数据结构-（5）列表

```R
g <- "My First List"
h <- c(25, 26, 18, 39)
j <- matrix(1:10, nrow=5)
k <- c("one", "two", "three")
mylist <- list(title=g, ages=h, j, k)

mylist[[2]]
mylist$age
mylist[["ages"]]
```

---

## Page 43：R数据结构-（6）因子

```R
#因子创建
fac = factor(c('A','B','A','B','AB','O','A','A'),  # 数据
             levels=c('A','B','AB','O'))  # 类别

heights = data.frame(height_cm=c(156,182,170),
                                    gender=c('f','m','f') )  # 数据框中含有字符的列为因子型
class(heights$height_cm)  # 返回numeric
class(heights$gender)  # 返回factor

#因子引用
fac[1]  # 因子水平不改变
```

## Page 44： 对象和它的模式与属性

```R
a <- c(1,2,3,4,5,6)              
b<- matrix(1:20, nrow=5, ncol=4, byrow=T)
mode(a)
mode(b)
class(a)
class(b)
```

---

## Page 45：对象类型转换

```R
mode()
class()
is.numeric() #返回值为TRUE或FALSE
is.logical()
is.charactor()
is.data.frame()
as.numeric() #转换为数值型
as.logical()
as.charactor()
as.matrix()
as.data.frame()
```

##  Page 48：数据的输入和输出

```R
FPKM<-read.delim("FPKMs.xls",header=T,sep="\t",check.names = F)
# FPKM<-read.table("D:/MyRProject/FPKMs.xls",header=T,sep="\t",check.names = F)
# 当然也可以输入绝对路径，如果您是在D盘新建了一个MyRProject文件夹的话
FPKMnew<-log(FPKM[,-1]) 
write.table(FPKMnew, "FPKMnew.txt",sep="\t",col.names=T,row.names=F,quote=F)

FPKM<-read.table("FPKMs.xls",header=T,sep="\t",row.names=1)
FPKMnew <- log(FPKM)     #取log计算
genename <- rownames(FPKM)
FPKMnew<- cbind(genename,FPKM)
write.table(FPKMnew,"logFPKM.txt",sep="\t",col.names=T,row.names=F,quote=F)
```

## Page 49：数据保存与加载

```R
# 关闭Rstudio-保存数据（注意路径）
#.Rdata 保存代码运行的数据对象
#save(FPKM,file = “.Rdata”) #默认保存
load(“.Rdata”)
# .Rhistory 保存代码运行的历史命令
load (“.Rhistory.R”)
```

  ## Page 51：数据合并与重塑- tidyr

```r
install.packages("tidyr")
library(tidyr)
#读取数据
expression <- read.delim("gene_expression.txt",sep="\t",header = T)
#短数据变成长数据
tidy_gather <- gather(expression,key=Samplename,value = Expression,-id)
#长数据变短数据
tidy_spread <- spread(tidy_gather,key=Samplename,value = Expression)
```

## Page 52 ：数据合并与重塑- tidyr

```R
#按列分割
tidy_separate <- separate(tidy_gather,col=Samplename,into=c("Source","Samplename"),sep="_")
#按列合并
tidy_unite <- unite(tidy_separate,col=Samplename,into=c("Source","Samplename"),sep="_")
```

___

## Page 53：数据合并与重塑- dplyr

```R
install.packages("dplyr")
library(dplyr)
#数据排序
##按id进行排序
dplyr_arrange <- arrange(tidy_gather,id) 
##按id进行排序的基础上按Expression的降序排列
dplyr_arrange1 <- arrange(tidy_gather,id,desc(Expression)) 
```

## Page 54：数据合并与重塑- dplyr

```R
#筛选表达量大于1的行
dplyr_filter <- filter( tidy_gather, Expression>1)
#筛选表达量大于1且SampleName 为 gene1的行
dplyr_filter <- filter( tidy_gather , Expression >1 & id == "gene1")
```

## Page 55：数据合并与重塑- dplyr

```R
tmp <-  filter( tidy_gather , Expression>1 ) 
Result <- arrange( tmp , Expression )
#使用管道操作符 直接把数据传递给下一个函数调用
Result <- filter( tidy_gather , Expression>1 ) %>% arrange( Expression )
```

   ## Page 56：数据合并与重塑- dplyr

```R
#select指定列
#展示指定的GeneId SampleName  Expression 列
dplyr_select <- select( tidy_separate,id,Samplename,Expression )
dplyr_select <- select(tidy_separate,-Source)
```

___

## Page 57：数据合并与重塑- dplyr

```r
#mutate增加新列ID
dplyr_mutate <- mutate(tidy_gather,ID=sub("gene","Gene",id))

```

## Page 58：高级软件包dplyr

```R
#增加expression列是原来的Expression的10倍，并替换原来的Expression列
dplyr_mutate_select <- mutate(tidy_gather,expression=Expression*10)%>%select(-Expression)
```

## Page 59：数据合并与重塑- dplyr

```R
#分组统计平均值
dplyr_groupby_summarise<-group_by(tidy_gather,Samplename)%>%summarise(avg=mean(Expression))
```

## Page 61：数据合并与重塑- dplyr

```R
a <- read.table("a.txt",header=T,sep="\t")
b <- read.table("b.txt",header=T,sep="\t")
c <- read.table("c.txt",header=T,sep="\t")

bind_rows(a,c)
bind_cols(a,c)
union(a,c)
setdiff(a,c)
```

## Page 62：数据合并与重塑- dplyr

```R
inner_join(a,b,by="x1")
full_join(a,b,by="x1")
left_join(a,b,by="x1")
```

##  Page 63：数据合并与重塑- dplyr

```R
right_join(a,b,by="x1")
semi_join(a,b,by="x1")
anti_join(a,b,by="x1")
```

## Page 72 ：练习题：折线图-plot()

```R
dose <- c(20, 30, 40, 45, 60)
drugA <- c(16, 20, 27, 40, 60)
pdf("plot.pdf",width=12,height=9)
par(mar=c(2,2,2,2),lwd=2, cex=1.5,cex.axis=.75, font.axis=3)

plot(dose, drugA, type="b", pch=19, lty=2, col="red")
dev.off()


dose <- c(20, 30, 40, 45, 60)
drugB <- c(15, 18, 25, 31, 40)
par(mar=c(2,2,2,2),lwd=2, cex=1.5,cex.axis=.75, font.axis=3)
plot(dose, drugB, type="b", pch=23, lty=6, col="blue", bg="green")
```

## Page 74 ：自定义函数

```R
mean_function <- function(a){
  b = sum(a)/length(a)
  return (b)
}

a <- c(10,20,30)
mean(a)
mean_function(a) 
```

------

## Page 76 ：自定义函数

```R
plot_function <- function(x,y){
  plot(x,y)
  return(x+y)
}
x <- seq(10,100,by=10)
y <- seq(1,10,by=1)
plot_function(x,y)
```

------

## Page 77 ：自定义函数

```R
sqtest <- function(x,y){
  z1=x^2;z2=y^2;z3=z1+z2;
  return(z3);
}
sqtest(3,4)
```

## Page 78 ：自定义函数

```R
# 加载外部自定义函数
#rm() # 删除某一指定对象
#rm = (list = ls()) ##清空环境变量
source('sqtest.R')  #默认当前路径
sqtest(3,4)
```

## Page 79：判断

```R
if(x<8){
  x <- x + 4
}
if(x>=8){
  x <- x - 4
}

if(x<8){
  x <- x + 4
}else{
  x <- x - 4 
}
```

------

## Page 80：循环

```R
x <- c(3,4,5)
for( n in x){
    print (n^2)
}

total <- 0
for (i in 1:10){
  total <- total + i
}
total
```

------

## Page 81 ：for循环应用—批处理

```R
#生成命名为1-10个以csv结尾的文件
for (i in 1:10){
  write.csv(i,paste(i,".csv",sep=""),row.names = F)
}
# 列出文件名
filename1 <- list.files(pattern=".csv$")

```

## Page 82 ：批量修改文件名称

```R
# 目标文件文件名
filename2 <- paste("A",filename1,sep="")
# 用file.rename函数修改
file.rename(filename1,filename2)
```

## Page 84 ：apply

```R
A<-matrix(1:12,c(3:4))
A
apply(A, MARGIN=1, FUN=sum)  # 按行操作：求和
apply(A, MARGIN=2, FUN=sum)  # 按列操作：求和
apply(A, MARGIN=1, function(x) sum(x)+4) 
```

## Page 85：lapply

```R
A<-matrix(1:12,c(3,4))
A
A.df<-data.frame(A)
is.list(A.df)
lapply(A.df, function(x) x+4)
lapply(A.df, function(x) sum(x)+4)
```

## Page 86：sapply

```R
B<-sapply(A.df, function(x) x^2)
B
C<-sapply(A.df, function(x,y) x^2+y, y=4)
C
D<-sapply(A.df, sum)
D
```

## Page 87：几种循环效率比较

```R
#清空环境变量
rm(list=ls())
#生成数据集
x <- cbind(x1=2, x2 = c(800:1, 2:1000))
head(x)
```

##  Page 88：几种循环效率比较

```R
#apply方法构建函数
apply_test<-function(x){
     myFUN<- function(x, c1, c2) {
         c(sum(x[c1],1), mean(x[c2])) 
       }
     apply(x,1,myFUN,c1='x1',c2=c('x1','x2'))
   }

#for循环方法构建函数
for_test<-function(x){
     df<-data.frame()
     for(i in 1:nrow(x)){
         row<-x[i,]
         df<-rbind(df,rbind(c(sum(row[1],1), mean(row))))
       }
   }

#向量化编程构建函数
vector_test<-function(x){
     data.frame(x1=x[,1]+1,x2=rowMeans(x))
   }
```

## Page 89：几种循环效率比较

```R
#分别统计3种方法的CPU耗时
system.time(apply_test(x))
system.time(for_test(x))
system.time(vector_test(x))
```











