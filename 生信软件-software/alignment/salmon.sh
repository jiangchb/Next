LD_LIBRARY_PATH=/home/mxy/software/salmon/salmon-latest_linux_x86_64/lib/:$LD_LIBRARY_PATH
/home/mxy/software/salmon/salmon-latest_linux_x86_64/bin/salmon index -p 5 -t /public/mid/rna/Ref_Genome_2018/human/Homo_sapiens_assembly_GRCh38.p12/mRNA.fa -i mRNA
ls clean_data/*.R1.fq.gz |sed 's:clean_data/::;s/.R1.fq.gz//' >sample.list
# -l ISR 第一链特异性
# -l IU   非链
#-l A 自动检测
while read a
do
/home/mxy/software/salmon/salmon-latest_linux_x86_64/bin/salmon quant -i mRNA -l A -1 clean_data/${a}.R1.fq.gz -2 clean_data/${a}.R2.fq.gz -p 8  -g mRNA_2_gene.xls  -o ${a}_quant
cut -f5 ${a}_quant/quant.genes.sf |sed "s/NumReads/$a/" >${a}_quant/tmp.count
done<sample.list
top=`head -n1 sample.list`
cut -f1 ${top}_quant/quant.genes.sf |sed 's/Name/id/' >gene.list
paste gene.list */tmp.count >gene_counts.xls
cut -f1,2 ${top}_quant/quant.genes.sf |sed 's/Name/id/;s/Length/len/' >gene_len.xls
python counts2fpkm.py -c gene_counts.xls -l gene_len.xls  -o gene_fpkm.xls


