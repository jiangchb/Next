GO的数据挖掘

#根据GO的分类来进行数据挖掘

对于Gene ontology 而言，目前共有2万多个Go trems。 做完富集分析后，我们可能会得到几百甚至几千个富集到的GO terms, 这样的一个数据量对于人工一个个检索而言，仍然是一个艰巨的任务。为了有效的利用GO富集分析的结果，我们势必需要对结果再次进行过滤。

这样的结构我们称之为有向无环图DAG, 虽然在图这种数据结构中，节点并没有严格的层级关系，但是由于在GO这张图中，存在了祖先节点，即最上层的3个节点，其他的节点都可以看做是其子节点，从而引用了树状结构中的level的概念，定义子节点到祖先节点的路径上包含的节点数即为该节点的level，祖先节点的level为1. 示意图如下

需要注意的是，由于子节点到祖先节点的路径不止一条，所以一个子节点可能拥有用多个level, 这意味着GO terms的level不是一个值，在使用level对GO Terms进行过滤时就需要注意。

想象一下，对于一个有多个level的GO term, 我们采用哪个值来表征其level, 是取最大值，还是最小值，或者是均值，由于不同取值算法带来的不确定性，所以采用level对GO过滤会存在一定风险，特别是level很大时。比如我们只选取level > 7的GO terms, 无论是用哪个值作为level, 其过滤的结果和我们预期的都是不符合的。

GO官网对于GO level也进行了说明，参考以下链接

http://www.geneontology.org/faq/how-can-i-calculate-level-go-term

传统的费舍尔精确检验也好，GSEA也罢，这些富集分析的算法都只是为单个GO term进行分析，不会考虑该GO term在整个网状结果中的层级关系。对于这些分析的结果，采用上述的GO level 进行过滤时，只能是采用较小的level, 在一下R包中，比如goprofiler， 推荐的最小层级是level为2。

采用level对结果过滤效果有限，为了有效筛选结果，出现了Gene Ontolgy network analysis，示意图如下


根据所有富集到的GO terms, 从整个GO Graph中取出一个子图subgraph， 图中有颜色的节点为富集到的GO, 颜色的深浅有P值决定， 节点的大小由degree决定。

根据这个network, 应用图论的算法可以挖掘出其中重要的GO term，从而实现对GO富集结果的过滤。

图来源：http://www.360doc.com/content/19/1225/20/68068867_882179240.shtml


library("GO.db")---gene2level2   分析--可视化--柱状图
library("")
library("topGO")---level—all  分析--可视化--有向无环图


#根据单个富集到的GO条目来进行挖掘
传统的费舍尔精确检验 --- 可视化就是咱们常规的柱状图和气泡图
GSEA  ---- 更加准确
