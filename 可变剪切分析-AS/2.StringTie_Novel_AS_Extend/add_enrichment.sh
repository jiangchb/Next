#/bin/bash
# by wsb 2019.08.01

SCRIPT=$1
input=$2
mRNA_2_gene=$3

cd $input

for i in `ls -1d *AS_result.xls`;do
    n=${i%.xls}
    perl $SCRIPT/add_enrichment.pl ${mRNA_2_gene} ${i} ${n}.txt 
    mv ${n}.txt  ${i}
    sed -i 's:();::g' ${i}
done

