#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Cwd;
use FindBin qw($Bin);
use File::Basename;

####################### software or env ################################################
my $bowtie_bin="/home/fanyucai/software/Anaconda3/bin";
my $express="export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v4.9.4/lib64/:\$LD_LIBRARY_PATH && /home/fanyucai/software/eXpress/express-1.5.1-linux_x86_64/express";
my $qsub_pbs="/home/shenyitian/scripts/qsub/qsub-sge.pl";
########################################################################################

my ($cleandatadir, $sample, $unigene, $lib_type, $gene_type, $outdir, $work_sh, $qsub, $cpu, $queue, $thread_number, $help);
GetOptions(
	"g:s"     => \$unigene,
	"s:s"     => \$sample,
	"d:s"     => \$cleandatadir,
	"t:s"     => \$lib_type,
	"w:s"     => \$gene_type,
	"o:s"     => \$outdir,
	"e:s"     => \$work_sh,
	"q:s"     => \$qsub,
	"c:i"     => \$cpu,
	"qu:s"    => \$queue,
	"p:i"     => \$thread_number,
	"h|help!" => \$help,
);

my $usage=<< "USAGE";

Program: $0
Description: run mapping(bowtie2), fpkm and counts(express) analysis
Options:
	-g      <infile>   The reference sequence(Unigene) fasta file
	-s      <str>      The sample name of the same assembly group, eg:Sample_A1,Sample_A2
	-d      <dir>      The input directory of cleandata
	-t      <str>      The library type: unstrand,RF,FR. [default: RF]
	                   unstrand: standard illumina RNA sequencing
	                   RF: rf-stranded, Strand-specific RNA sequencing(dUTP method)
	                   FR: fr-stranded, Strand-specific RNA sequencing
	-w      <str>      The gene type: Whether the reference sequence is whole. [default: yes]
	-o      <dir>      The output directory of result. [default: ./]
	-e      <dir>      The output directory of shell scripts.[default: ./]
	-q      <Y|N>      Whether qsub task? If Y, will qsub. If N, will only generate shell scripts. [default: Y]
	-c      <num>      cpu number. [default: 4]
	-qu     <str>      queue:cu. [default: cu]
	-p      <num>      max thread number. [default: 5]
	-h|help            print help info
Example:
	perl 1.bowtie2_express.pl -g Unigene.fa -s Sample_A1,Sample_A2 -d ./clean_data -o out/

USAGE

die $usage if(!$unigene || !$sample || !$cleandatadir || $help);
$outdir ||=getcwd;
$work_sh ||= getcwd;
$lib_type ||= "RF";
$gene_type ||= "yes";
$qsub ||= "Y";
$cpu ||= 4;
$queue ||= "cu";
$thread_number ||= 10;
(-s $unigene) || die "Error: don't find unigene file: $unigene !\n";
$unigene=File::Spec->rel2abs($unigene);
(-d $cleandatadir) || die "Error: don't find cleandata directory: $cleandatadir !\n";
$cleandatadir=File::Spec->rel2abs($cleandatadir);$cleandatadir=~s/\/$//g;
(-d $outdir) || mkdir $outdir;
$outdir=File::Spec->rel2abs($outdir);$outdir=~s/\/$//g;
(-d $work_sh) || mkdir $work_sh;
$work_sh=File::Spec->rel2abs($work_sh);$work_sh=~s/\/$//g;
my $lib_type_opt;
if($lib_type eq "unstrand"){
	$lib_type_opt="";
}elsif($lib_type eq "RF"){
	$lib_type_opt="--rf-stranded";
}elsif($lib_type eq "FR"){
	$lib_type_opt="--fr-stranded";
}else{
	print "Error: Please check options '-t' !\n";
	die $usage;
}
my $gene_type_opt;
if($gene_type eq "yes"){
	$gene_type_opt="";
}elsif($gene_type eq "no"){
	$gene_type_opt="--no-bias-correct";
}else{
	print "Error: Please check options '-w' !\n";
	die $usage;
}


############# check sample cleandata fastq
my @sample=split /,/,$sample;
my %sample;
for(my $i=0; $i<=$#sample; $i++){
	my $sample_R1="$cleandatadir/$sample[$i].R1.fq.gz";
	my $sample_R2="$cleandatadir/$sample[$i].R2.fq.gz";
	(-s $sample_R1) || die "Error: don't find file: $sample_R1 !\n";
	(-s $sample_R2) || die "Error: don't find file: $sample_R2 !\n";
	$sample{$sample[$i]}{R1}=$sample_R1;
	$sample{$sample[$i]}{R2}=$sample_R2;
}

########### run bowtie2 and express
my $unigene_filename=basename($unigene);
open SH,">$work_sh/a.build_index.sh" || die $!;
print SH "$bowtie_bin/bowtie2-build $unigene $outdir/$unigene_filename\n";
close SH;

# --no-bias-correct 
open SH1,">$work_sh/b.bowtie2_express.sh" || die $!;
my $str;my $str_name;
for my $s (@sample){
	$str.="$outdir/$s.bowtie2.log,";
	$str_name.="$s,";
	print SH1 "$bowtie_bin/bowtie2 --reorder -k30 -t -p 20 -x $outdir/$unigene_filename -1 $sample{$s}{R1} -2 $sample{$s}{R2} -S $outdir/$s.hits.sam 2>$outdir/$s.bowtie2.log";
	print SH1 " && $express --no-update-check $gene_type_opt $lib_type_opt -o $outdir/$s\_express $unigene $outdir/$s.hits.sam";
	print SH1 " && mv $outdir/$s\_express/results.xprs $outdir/$s.results.xprs && rm -fr $outdir/$s\_express/ && ";
	print SH1 "cut -f 2,8,11 $outdir/$s.results.xprs |awk 'NR == 1; NR > 1 {print \$0 | \"sort -k1,1n\"}' > $outdir/uniq_counts.fpkm.$s.txt && rm $outdir/$s.hits.sam\n";
}
close SH1;

$str=~s/,$//g; $str_name=~s/,$//g;
open SH2,">$work_sh/c.mapping_stat.sh" || die $!;
print SH2 "perl $Bin/mapping_stat.pl -i $str -s $str_name -o $outdir/mapping_stat_transcriptome.xls\n";
close SH2;

open SH3,">$work_sh/d.merge_gene.sh" || die $!;
print SH3 "cd $outdir && perl $Bin/merge_counts_fpkm.pl -i ./ -o ./tmp_results ";
close SH3;

if($qsub eq "Y"){
	system("perl $qsub_pbs --queue $queue --lines 1 --num_proc $cpu $work_sh/a.build_index.sh");
	system("perl $qsub_pbs --queue $queue --lines 1 --maxjob $thread_number --num_proc $cpu $work_sh/b.bowtie2_express.sh");
	system("perl $qsub_pbs --queue $queue --lines 1 --num_proc 1 $work_sh/c.mapping_stat.sh");
	system("perl $qsub_pbs --queue $queue --lines 1 --num_proc 1 $work_sh/d.merge_gene.sh");
}
