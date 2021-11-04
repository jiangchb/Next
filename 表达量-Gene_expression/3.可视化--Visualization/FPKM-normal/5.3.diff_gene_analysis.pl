#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin);
use File::Basename;
use Cwd;

####################### software or env ##################################################
my $qsub_pbs="/data/software/qsub/qsub-sge.pl";
my $env="export PATH=/data/software/gcc/gcc-v6.4.0//bin/:/data/software/Anaconda3/envs/mro-v3.4.3/bin/:\$PATH";
$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
########################################################################################

my ($indir, $outdir, $work_sh, $sample_num, $annotation, $sample_group, $diff_group, $pval, $fdr, $foldchange, $cpu, $queue, $help);
GetOptions(
	"i:s"        => \$indir,
	"o:s"        => \$outdir,
	"a:s"        => \$annotation,
	"n:i"        => \$sample_num,
	"sg:s"       => \$sample_group,
	"dg:s"       => \$diff_group,
	"p:f"        => \$pval,
	"f:f"        => \$fdr,
	"c:f"        => \$foldchange,
	"e:s"        => \$work_sh,
	"t:i"        => \$cpu,
	"q:s"        => \$queue,
	"h|help!"    => \$help,
);
my $usage=<< "USAGE";
Program: $0
Description: difference gene analysis
Options:
	-i     <indir>    The input directory of uniq_counts.fpkm.*.txt         [Required]
	-o     <outdir>   The output directory of result                        [Required]
	-a     <infile>   annotation.xls file                                   [Required]
	-n     <num>      The sample number                                     [Required]
	-dg    <infile>   Input difference group file                           [Optional]
	-sg    <infile>   Input sample group file                               [Optional]
	-p     <num>      Pvalue ratio threshold.                               [Optional]
	-f     <num>      fdr(Qvalue) ratio threshold.                          [Optional]
	                  Filtering can be performed using any one of (-p), (-f) at a time
	-c     <num>      foldchange threshold. [default: 2]                    [Optional]
	-e     <outdir>   The output directory of shell scripts.[default: ./]   [Optional]
	-t     <num>      cpu num. [default: 3]                                 [Optional]
	-q     <str>      queue:all,fat,big. [default: all]                     [Optional]
	-h|help           print help info
Example:
	perl 4.3.diff_gene_analysis.pl -i bowtie2_express/ -o Quantification -a annotation.xls -n 4 -sg sample_group.txt -dg diff_group.txt -p 0.05

USAGE

die $usage if(!$indir || !$outdir || !$annotation || !$sample_num || $help);
if(!$diff_group && $sample_num>1){
	die "Error: When sample number >1 , option '-dg' is required !!!\n";
}
if((!$pval && !$fdr) || (defined $pval && defined $fdr)){
	warn "Error: Filtering can be performed using any one of (-p), (-f) at a time!\n";
	die $usage;
}
$cpu ||= 3;
$queue ||= "cu";
$foldchange ||= 2;
$work_sh ||= getcwd;
(-d $indir) || die "Error: don't find indir: $indir !\n";
(-d $outdir) || mkdir $outdir;
(-d $work_sh) || mkdir $work_sh;
(-s $annotation) || die "Error: don't find annotation file: $annotation !\n";
$indir =File::Spec->rel2abs($indir);
$outdir=File::Spec->rel2abs($outdir);
$work_sh=File::Spec->rel2abs($work_sh);
$annotation=File::Spec->rel2abs($annotation);

if(defined $sample_group){
	(-s $sample_group) || die "Error: don't find sample group file: $sample_group !\n";
	$sample_group=File::Spec->rel2abs($sample_group);
}
if(defined $diff_group){
	(-s $diff_group) || die "Error: don't find difference group file: $diff_group !\n";
	$diff_group=File::Spec->rel2abs($diff_group);
}
my @counts_fpkm=glob "$indir/uniq_counts.fpkm.*.txt";
(@counts_fpkm!=$sample_num) && die "Error: don't find $sample_num uniq_counts.fpkm.*.txt in $indir !\n";

open SH,">$work_sh/a.diff_gene.sh" || die $!;
print SH "cd $outdir && perl $Bin/merge_counts_fpkm.pl -i $indir -o $outdir && $env";
print SH " && Rscript $Bin/fpkm_stati.r -i fpkm.xls -o fpkm_boxplot  -t fpkm -s $sample_group -o $outdir";
print SH " && perl $Bin/expression_density_bar.pl -i fpkm.xls -t fpkm -g $sample_group -o $outdir" if(defined $sample_group);
print SH " && perl $Bin/expression_density_bar.pl -i fpkm.xls -t fpkm -o $outdir" if(!defined $sample_group);
if($sample_num > 2){
	print SH " && Rscript $Bin/pca_normal.r -i counts.xls -s $sample_group" if(defined $sample_group);
	print SH " && Rscript $Bin/pca_normal.r -i counts.xls" if(!defined $sample_group);
}
if($sample_num > 1){
	print SH " && Rscript $Bin/sample_corrplot.r -i fpkm.xls -o heatmap-coefficient_matrix";
	print SH " && Rscript $Bin/Run_DEG.R -e fpkm.xls -r counts.xls -d $diff_group -a $annotation -p $pval -c $foldchange -t gene -o $outdir" if(defined $pval);
	print SH " && Rscript $Bin/Run_DEG.R -e fpkm.xls -r counts.xls -d $diff_group -a $annotation -f $fdr -c $foldchange -t gene -o $outdir" if(defined $fdr);
	print SH " && perl $Bin/diffgene_venn.pl -i *-vs-*-diff-*.xls -o $outdir -a $annotation";
}else{
	print SH " && perl $Bin/add_annotation.pl -i fpkm.xls,counts.xls -a $annotation -o $outdir";
}
print SH "\n";
close SH;
system("perl $qsub_pbs --queue $queue --lines 1 --num_proc $cpu $work_sh/a.diff_gene.sh");
