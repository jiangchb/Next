path="/public/cluster2/works/intern/ccj_temp1/follow_up/Gene_expression_pipeline/TPM"
#counts_anno 转换为gene_counts_new.xls

#counts 转换 tpm
/home/rna/bin/python3 $path/counts2tpm.py -c gene_counts_new.xls -l gene_len.txt -o TPM.xls
#TPM 为id

#tpm表达量图
perl $path/expression_density_bar.pl -i TPM.xls -T tpm -g sample_group.txt -o ./1.3.gene_expression/
#sed "1i/id" TPM.xls > TPM.xls
Rscript $path/fpkm_stati.r -i TPM.xls -s sample_group.txt -t TPM -o ./1.3.gene_expression/

#注释

#提取注释信息
num=`head -1 gene_counts.xls| awk -F '\t' '{for (i=1;i<=NF;i++) {if ($i == "Dbxref") {print i}}}'`
cut -f 1,$num-1000 gene_counts.xls > anno.txt
#进行注释
/home/rna/bin/python3 /public/cluster2/works/intern/ccj_temp1/follow_up/TFs/concat.py -f TPM.xls -f2 anno.txt -r1 ./1.3.gene_expression/TPM_anno.xls -index_name id

#准备所有对照组的差异文件以及TPM.xls sample.xls文件放在一个目录下面，然后运行
echo "======================== 正在进行热图绘制============================="
mkdir -p ./diff_FC_2
mkdir -p ./diff_FC_1.5

cp ./1.1.different_expressed_gene/*-vs-*-diff*.xls ./diff_FC_1.5
cp ./1.4.different_expressed_gene/*-vs-*-diff*.xls ./diff_FC_2
cp ./TPM.xls ./diff_FC_2
cp ./TPM.xls ./diff_FC_1.5
cp ./sample_group.txt ./diff_FC_1.5
cp ./sample_group.txt ./diff_FC_2
cat ./anno.txt | sed '1s/id/gene_id/' > ./diff_FC_1.5/anno.txt
cat ./anno.txt | sed '1s/id/gene_id/' > ./diff_FC_2/anno.txt

cd ./diff_FC_2
sh /public/cluster2/works/intern/ccj_temp1/follow_up/Gene_expression_pipeline/further-heatmap-for-diff.sh
for i in ./*-vs-*/; do /home/rna/bin/python3 /public/cluster2/works/intern/ccj_temp1/follow_up/TFs/concat.py -f $i/TPM.genes_reorder.cluster_result.xls -f2 ./anno.txt -r $i/TPM_diff_anno.xls -index_name gene_id; done

cd ../diff_FC_1.5
sh /public/cluster2/works/intern/ccj_temp1/follow_up/Gene_expression_pipeline/further-heatmap-for-diff.sh
for i in ./*-vs-*/; do /home/rna/bin/python3 /public/cluster2/works/intern/ccj_temp1/follow_up/TFs/concat.py -f $i/TPM.genes_reorder.cluster_result.xls -f2 ./anno.txt -r $i/TPM_diff_anno.xls -index_name gene_id; done
##准备注释文件
echo "======================== 完成 ============================="