#无法输入sample文件
#用于差异文件批量跑热图
#统计样本个数
for i in *.xls; do 
sample_number=`awk '{for(i = 1; i <= 100; i++) if ($i ~ /^Expression_/) print $i }' $i | wc -w` ; #求出sample数量
cut -f 1,10-$((9+$sample_number)) $i | sed "s/Expression_//g" > fpkm_$i;
Rscript /public/cluster2/works/intern/ccj_temp1/follow_up/Heatmap/further_heatmap_v2.R -e fpkm_$i -l F -o Heatmap_$i -k 5 -d sample_group.txt -g 10;done
