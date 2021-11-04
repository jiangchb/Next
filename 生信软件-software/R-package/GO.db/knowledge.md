对于一个数据表而言，首先我们需要知道表头信息，就可以通过columns和keytypes函数来访问得到，示例如下
对于一个数据表而言，首先我们需要知道表头信息，就可以通过columns和keytypes函数来访问得到，示例如下

> keytypes(GO.db)
[1] "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"      
> columns(GO.db)
[1] "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"

从以上代码的结果可以看出，GO.db提供的数据表共有4列。

假如想要访问其中某一列的值，可以通过keys函数，示例如下

> keys(GO.db, keytype = "GOID")[1:3]
[1] "GO:0000001" "GO:0000002" "GO:0000003"

上述代码返回GOID这一列的内容。

对于数据库而言，查询是基本操作，在SQL语言中，通过select实现，对应的在R中通过select函数来实现，示例如下

> k <- keys(GO.db, keytype = "GOID")[1:3]
> select(GO.db,
   keys = k,
   columns = c("TERM","ONTOLOGY"),
   keytype="GOID")

'select()' returned 1:1 mapping between keys and columns
        GOID                             TERM ONTOLOGY
1 GO:0000001        mitochondrion inheritance       BP
2 GO:0000002 mitochondrial genome maintenance       BP
3 GO:0000003                     reproduction       B

通过返回结果可以看到，GO.db提供了一张4列的数据表，GOID表示GO编号，DEFINITION表示GO功能的详细描述信息，TERM表示功能的简单介绍，ONTOLOGY表示GO的3大类别。

除了基本的数据表之外，在这种类型的包中还会提供很多其他信息，可以通过ls函数查看，示例如下

> ls("package:GO.db")
[1] "GO"            "GO.db"         "GO_dbconn"     "GO_dbfile"     "GO_dbInfo"     "GO_dbschema"  
[7] "GOBPANCESTOR"  "GOBPCHILDREN"  "GOBPOFFSPRING" "GOBPPARENTS"   "GOCCANCESTOR"  "GOCCCHILDREN"
[13] "GOCCOFFSPRING" "GOCCPARENTS"   "GOMAPCOUNTS"   "GOMFANCESTOR"  "GOMFCHILDREN"  "GOMFOFFSPRING"
[19] "GOMFPARENTS"   "GOOBSOLETE"    "GOSYNONYM"     "GOTERM"

其中有一部分对象的类型AnnDbBimap, 示例如下

> GOTERM
TERM map for GO (object of class "GOTermsAnnDbBimap")  

这种对象类似基本数据结构中的list, 常用的操作语句示例如下

> mappedkeys(GOTERM)[1:3]
[1] "GO:0000001" "GO:0000002" "GO:0000003"

> ls(GOTERM)[1:3]
[1] "all"        "GO:0000001" "GO:0000002"


> GOTERM[["GO:0000001"]]
GOID: GO:0000001
Term: mitochondrion inheritance
Ontology: BP
Definition: The distribution of mitochondria, including the mitochondrial genome, into daughter
    cells after mitosis or meiosis, mediated by interactions between mitochondria and the
    cytoskeleton.
Synonym: mitochondrial inheritance

> get("GO:0000001", GOTERM)
GOID: GO:0000001
Term: mitochondrion inheritance
Ontology: BP
Definition: The distribution of mitochondria, including the mitochondrial genome, into daughter
    cells after mitosis or meiosis, mediated by interactions between mitochondria and the
    cytoskeleton.
Synonym: mitochondrial inheritance


> mget("GO:0000001", GOTERM)
$`GO:0000001`
GOID: GO:0000001
Term: mitochondrion inheritance
Ontology: BP
Definition: The distribution of mitochondria, including the mitochondrial genome, into daughter
    cells after mitosis or meiosis, mediated by interactions between mitochondria and the
    cytoskeleton.
Synonym: mitochondrial inheritance

1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
ls和mappedkeys函数都是用于查看这个列表的名称，只不过ls会对所有key排序;get和mget选取其中的内容，也可以像list一样，通过[[ ]]操作符直接访问。

由于和list类似，所以经常会将这些对象通过as.list转换之后，在进行操作，示例如下

> go <- as.list(GOTERM)
> go[[1]]
GOID: GO:0000001
Term: mitochondrion inheritance
Ontology: BP
Definition: The distribution of mitochondria, including the mitochondrial genome, into daughter
    cells after mitosis or meiosis, mediated by interactions between mitochondria and the
    cytoskeleton.
Synonym: mitochondrial inheritance
1
2
3
4
5
6
7
8
9
需要注意的是这个步骤是非常耗时的，实际使用时，可以先挑选子集，然后在转换成list。

很多做GO富集分析的R包都会调用GO.db, 掌握其基本操作，有助于理解其他封装好的R包。
