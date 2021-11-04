#! /bin/bash

indir=$1
diff_file=$2
outdir=$3
DirForScriptSelf=$(cd "$(dirname "$0")";pwd)

mkdir ./tmp
mkdir ./circle_plot
cp $indir/top_300_diff*.txt ./tmp
cp $diff_file/*-vs-*-diff-*FC-*.xls .

wc -l ./tmp/* |
    awk '{print $1"\t"$2}'|
    while read a b ;do 
    if [ $a -lt 2 ];
        then echo $b" file empty,delete";
        rm $b;
    fi;
    done
cp ./tmp/* ./

ls top_300_diff* 
for i in top_300_diff*.txt;do n=$(basename $i|sed 's|_gene2gene_network.txt||;s|top_300_diff-||');
	/root/miniconda3/bin/Rscript $DirForScriptSelf/ppi_circle.r \
	-i $i \
	-d ${n}-diff-*FC-*.xls \
	-n 30 \
	-o ./circle_plot
done

mv circle_plot/* $outdir
