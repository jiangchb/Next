#!/bin/bash
###
### Bowtie_anno
### Author:Congjia.chen
### Usage: Get Gene_annotation
###   sh Bowtie_anno2.sh <ANNO_PATH> <FPKM_PATH> <OUTPUT_PATH>
### Options:
###   <ANNO_PATH>   type the path where the ANNO_PATH IS
###   <FPKM_PATH>   type the path where FPKM is 
###   <OUTPUTPATH>   type the path where output is
###   -h        Show this message.

help() {
	awk -F'### ' '/^###/ { print $2 }' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi
ano_path=$1
path=$2
outpath=$3
cp $ano_path/KEGG/anno-kegg.backgroud.xls $outpath/gene_anno-kegg.background.xls
cp $ano_path/KEGG/kegg.backgroud.xls $outpath/gene_kegg.background.xls
cp $ano_path/GO/go.backgroud.xls $outpath/gene_go.background.xls
cp $ano_path/Swissprot/gene.Swissprot.blast.best.anno.xls $outpath/gene_swiss.background.xls
cp $ano_path/NR/gene.NR.blast.best.anno.xls $outpath/gene_nr.background.xls

cut -f 1,6 $outpath/gene_swiss.background.xls | sed "s/#GeneID/id/g">$outpath/swiss_tmp
cut -f 1,6 $outpath/gene_nr.background.xls | sed "s/#GeneID/id/g">$outpath/nr_tmp

# 要从genome.gff相应的位置提取出基因名字和type（比如：coding）和product信息。

#生成annotation文件对应的gene列表
#grep ">" gene.fa|sed "s/>//g"| sed '1i\id' > target_gene_list_final
cut -f 1 $path/fpkm.xls |sed  's/gene_id/id/g'> $outpath/target_gene_list_final

#生成annotation文件需要的go，kegg backgroud
cat $outpath/gene_kegg.background.xls| sed '1i\id\tKO\tKO anno' > $outpath/kegg_tmp
cat $outpath/gene_go.background.xls| sed '1i\id\tGO\tGO anno' > $outpath/go_tmp

#依次合并这些文件顺序为biotype-product-GO-KO
/home/rna/softwares/miniconda3/bin/python /public/cluster2/works/intern/chencongjia/Prok/script/Database_integrate.py -f $outpath/target_gene_list_final -f2 $outpath/nr_tmp -r $outpath/tmp1 -index_name id
/home/rna/softwares/miniconda3/bin/python /public/cluster2/works/intern/chencongjia/Prok/script/Database_integrate.py -f $outpath/tmp1 -f2 $outpath/go_tmp -r $outpath/tmp3 -index_name id
/home/rna/softwares/miniconda3/bin/python /public/cluster2/works/intern/chencongjia/Prok/script/Database_integrate.py -f $outpath/tmp3 -f2 $outpath/kegg_tmp -index_name id -r $outpath/gene_annotation.xls

rm $outpath/*tmp*
