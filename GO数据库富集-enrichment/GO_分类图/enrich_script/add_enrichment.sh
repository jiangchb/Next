#/bin/bash
SCRIPT=$1
input=$2
mRNA_2_gene=$3

cd $input

for i in `ls -1d enrichment*/*enrichment/*/`;do
    cd $i
    for a in `ls *xls`; do
        n=${a%.xls}
        perl $SCRIPT/add_enrichment.pl ${mRNA_2_gene} ${a} ${n}.txt 
        mv ${n}.txt ${a}
    done
    cd $input
done

