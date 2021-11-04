#R对于生信的重要性

在生信分析领域，R语言由于其简单易用的特点和良好的生态环境，占用重要的一席之地。其中，Bioconductor作为生信分析专用的R语言社区，提供了许多的R包。

software类型的R包用于执行某项具体的分析内容，比如edgeR, DESeq2等，AnnotationData类型的包在R中存储了对应的数据库，比如GO.db等，ExperimentData类型的包存储了实验数据，Workflow类型的包提供了完整分析的pipeline。本文主要介绍AnnotationData类型的包。

#AnnotationData类型
为了规范化开发，方便R包的使用，Bioconductor的开发者提供了几种基础的R包，用于定义几种基础信息的存储方式。

对于数据库内容的存储和使用，在AnnotationDbi这个包中统一进行了定义。由于采用了面向对象的编程方式，所有继承了这种对象的R包其使用方式是一样的。

在Bioconductor中，有以下4种类别的注释信息包，都继承了AnnotationDbi

Organism level
比如human对应的Org.Hs.eg.db, 存储了人类的基因信息

Platform level
比如hgu133plus2.db, 这种类型的包主要存储不同平台的数据，比如不同芯片的探针信息

Homology-level
比如hom.Dm.inp.db, 存储了同源信息

System-biology level
比如GO.db, 存储生物学相关的数据库

所有这些后缀为.db的R包，其本质都为一个sqlite数据库，一种轻量级的关系型数据库，只不过是通过R来进行访问。

以GO.db为例，在下载的源代码中，可以找到对应的后缀为.sqlite的数据库文件，位于extdata目录下。

关系型数据库中的基本单位是表，对于一个.db的R包而言，可以通过以下4个函数访问其中的内容

columns
keytypes
keys
select

