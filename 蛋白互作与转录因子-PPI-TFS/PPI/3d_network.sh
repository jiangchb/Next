#! /bin/bash
export PATH=/home/mxy/software/pandoc/pandoc-2.6/bin/:$PATH
indir=$1
TF=$2


cd $indir
for i in top_300_diff*.txt;do
n=${i%.txt}
awk 'BEGIN{FS=OFS="\t"}NR>1{print $1}' $i >> tmp1
awk 'BEGIN{FS=OFS="\t"}NR>1{print $3}' $i >> tmp1
sort tmp1 |uniq -c|sed 's/ /\t/g'| awk 'BEGIN{FS=OFS="\t"}{print $(NF-1),$NF}' |awk '{print $0"\t"(FNR-1)}'|awk 'BEGIN{FS=OFS="\t"}{print $2,$3,$1}'  > network.tmp && rm tmp1
awk 'BEGIN{FS=OFS="\t"}NR>1{print $1,$2}' $i >> tmp2
awk 'BEGIN{FS=OFS="\t"}NR>1{print $3,$4}' $i >> tmp2
###inside color
sort -u tmp2 | awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3;next}{if ($1==a[$1]) print $0"_TFs";else print $0}' $TF - |awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;next}{if (a[$1]==$1) print $0,b[$1]}' - network.tmp |sed 1'i\name\tid\tsize\tgroup' > ${i}_col  && rm tmp2
###relationship
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3;next}{if (a[$1]==$1) print $3,b[$1]}' network.tmp $i |awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=$1;b[$1]=$2;c[$1]=$3;next}{if (a[$1]==$1) print $2,b[$1],"1","black"}' network.tmp - | sed 1'i\prot1\tprot2\tscore\tcol'> ${i}_line && rm network.tmp
###networkD3
Rscript /public/cluster2/works/gonggy/script/8.PPi/networkD3.r -l ${i}_line -n ${i}_col -o $n.html
rm *_line *_col
done
cd -
