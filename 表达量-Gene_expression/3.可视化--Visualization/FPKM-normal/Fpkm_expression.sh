softwarepath="/public/cluster2/works/intern/ccj_temp1/follow_up/Gene_expression_pipeline/FPKM-normal/"
#"/home/shenyitian/pipeline/denovo_RNAseq_v3/5.Quantification_diffgene/" #常规大小
#"/public/cluster2/works/intern/ccj_temp1/follow_up/Gene_expression_pipeline/FPKM-mega" #巨图

perl $softwarepath/expression_density_bar.pl -i gene_fpkm_new.xls -T fpkm -g sample_group.txt -o ./1.3.gene_expression

Rscript $softwarepath/fpkm_stati.r -i gene_fpkm_new.xls -s sample_group.txt -t FPKM -o ./1.3.gene_expression
Rscript $softwarepath/pca_normal.r -i gene_counts_new.xls -s sample_group.txt
mv ./*boxplot* ./1.3.gene_expression
mv ./*PCA* ./1.3.gene_expression
mv ./sample2sample_distances_heatmap* ./1.3.gene_expression

Rscript $softwarepath/sample_corrplot.r -i gene_counts_new.xls -o ./1.3.gene_expression/heatmap-coefficient_matrix.pdf
perl $softwarepath/diff_stat_barplot.pl -i gene_diff_stat.xls -o ./1.4.different_expressed_gene
