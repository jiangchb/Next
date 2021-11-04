#!/bin/bash
###
### Gene Ontology Map
### Author:Congjia.chen
### Usage:
###   sh ref_script.sh <work_dir> <genome_dir> <groupname>
###
### Options:
###   <work_dir>   Input dir with enrichment-go-*-Total.xls
###   <genome_dir>   genome_dir with gene_go.backgroud.xls
###   -h        Show this message.

help() {
	awk -F'### ' '/^###/ { print $2 }' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi

work_dir=$1
genome_dir=$2
mkdir -p $work_dir/result
soft_path="/public/cluster2/works/intern/ccj_temp1/follow_up/GO_classification/enrich_script"
"/home/rna/RNA_works/ref_RNAseq_pipeline_v2.1/5.Quantification_DEG_Enrichment/enrich" #有参常规
"/home/lipeng/lipeng_script/pipline/Prok_Automata/5.Differential-expression/enrich/"

ori_dir=`pwd`

cd $work_dir
filename=`ls ./enrichment-go-*-Total.xls`
groupname=`echo $filename | awk {' gsub(/\.\/enrichment-go-/, "");gsub(/-Total.xls/,"");print '}`
echo $groupname

cd $ori_dir
perl $soft_path/diff_go_level2.pl -i $work_dir/enrichment-go-*-Total.xls -o $work_dir/result3 #生成GO.level2.stat.Total.xls
perl $soft_path/format_gobackgroud.pl -i $genome_dir/gene_go.backgroud.xls -s Unigene -c /home/lipeng/lipeng_script/pipline/Prok_Automata/5.Differential-expression/enrich/category.xls -o $work_dir/result3 && echo This-Work-is-Completed! #生成Unigene.GO.classification.stat.xls
perl $soft_path/compare_go_level2.pl -i $work_dir/result3/Unigene.GO.classification.stat.xls,$work_dir/result3/GO.level2.stat.Total.xls -o $work_dir/result3 -m ALL,DEG -n $groupname #画图
