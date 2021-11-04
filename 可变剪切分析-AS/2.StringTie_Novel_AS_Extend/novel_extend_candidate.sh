#!/bin/bash
# by The Coder, 2017-10-10

genome=$1 
cpc_run=$2
R="/home/fanyucai/software/R/R-v3.4.0/bin/R"
bedtools="/home/fanyucai/software/bedtools/bedtools2/bin/bedtools"
seqtk="/home/fanyucai/software/Seqtk/seqtk-master/seqtk"
perl="/home/fanyucai/software/perl/perl-v5.24.1/bin/perl"
NGSQCToolkit_PATH="/home/fanyucai/software/NGS_QC_Toolkit/NGSQCToolkit_v2.3.3"
export LD_LIBRARY_PATH=/home/fanyucai/software/zlib/zlib-v1.2.11/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:$LD_LIBRARY_PATH 

selfpath=$(cd "$(dirname "$0")";pwd)
cut -f 9 merged.gtf | sed 's/gene_id //g; s/ transcript_id //g; s/ exon_number //g; s/ oId //g' > merged_col9.txt

$R CMD BATCH $selfpath/getNovelisoform.r
awk -F "\t" '{printf("%s\t%.0f\t%.0f\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' OFS="\t" novel_isoform1.bed > novel_isoform.bed
$bedtools sort -i novel_isoform.bed > novel_isoform_sort.bed
awk -F "\t" '{printf("%s\t%.0f\t%.0f\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' OFS="\t" gene1_200.bed > gene_200.bed
$bedtools sort -i gene_200.bed > gene_200_sort.bed
$bedtools coverage -b gene_200_sort.bed -a novel_isoform_sort.bed > overlap.txt
awk -F "\t" '($7==0){print $4}' OFS="\t" overlap.txt > novel_isoform.list
$seqtk subseq cuffmerge.fa novel_isoform.list -l 80 > cpc0.fa
$perl $NGSQCToolkit_PATH/Trimming/AmbiguityFiltering.pl -p 10 -n 180 -i cpc0.fa -out cpc.fa
awk '/^>/{sub(">", "", $1); print $1}' cpc.fa |awk 'BEGIN{FS=OFS="\t"; print "chr","gene_id","isoform_id","+/-","blocks","Sizes","Starts"} \
NR==FNR{a[$1]=1; next} (a[$3]){ print $0}' - results_novel_isoform.txt > novel_isoforms.xls

$R CMD BATCH $selfpath/getBlatbed.r
sort Extendgene.bed | uniq > Extendgene1.bed 
sort Extend.ref.gene.bed | uniq > Extend.ref.gene1.bed
awk '{print $1,$2,$3,$4,"0",$6,$2,$3}' OFS="\t" Extendgene1.bed > Extendgene2.bed
awk '{print $1,$2,$3,$4,"0",$6,$2,$3}' OFS="\t" Extend.ref.gene1.bed > Extend.ref.gene2.bed
$bedtools getfasta -s -name -fi $genome  -bed Extendgene2.bed -fo Extendgene.fa
$bedtools getfasta -s -name -fi $genome  -bed Extend.ref.gene2.bed -fo Extend.ref.gene.fa

mkdir  ../Extend_gene
mv Extendgene1.bed Extend.ref.gene1.bed Extendgene.fa Extend.ref.gene.fa Extendgene_candidate.txt ../Extend_gene/

mkdir ../Novel_transcript
mv novel_isoforms.xls  cpc.fa  ../Novel_transcript 

