#无法输入sample文件
#用于差异文件批量跑热图
#统计样本个数

##############注意新旧版本交替！！！！！！！！DEseq这样
for i in *-vs-*.xls; do 
sample_number=`awk '{for(i = 1; i <= 100; i++) if ($i ~ /^Expression_/) print $i }' $i | wc -w` ; #求出sample数量
cut -f 1,10-$((9+$sample_number)) $i | sed "s/Expression_//g" > fpkm_$i;
Rscript /public/cluster2/works/intern/ccj_temp1/follow_up/Gene_expression_pipeline/TPM/further_heatmap_for_fpkm2tpm.R -e TPM.xls -i fpkm_$i -o ./ -k 6 -d sample_group.txt -g 8;done

##############DEseq2 列名可能会有变化 