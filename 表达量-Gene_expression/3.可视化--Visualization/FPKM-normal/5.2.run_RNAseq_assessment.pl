#!/usr/bin/env perl
use warnings;
use strict;
use File::Spec;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;

####################### software or env ################################################
my $qsub_pbs="/data/software/qsub/qsub-sge.pl";
########################################################################################

my ($unigene, $sample, $samdir, $outdir, $work_sh, $cpu, $queue, $thread_number, $help);
GetOptions(
	"r:s"      => \$unigene,
	"n:s"      => \$sample,
	"s:s"      => \$samdir,
	"o:s"      => \$outdir,
	"e:s"      => \$work_sh,
	"c:i"      => \$cpu,
	"q:s"      => \$queue,
	"p:i"      => \$thread_number,
	"h|help!"  => \$help,
);

my $usage=<< "USAGE";

Program: $0
Description: RNAseq assessment--sequencing saturation evaluation
Options:
	-r    <file>    The input unigene fasta file
	-n    <str>     The sample name of the same assembly group, eg:A1,A2
	-s    <dir>     The input sam file directory
	-o    <dir>     The output directory of result
	-e    <dir>     The output directory of shell scripts. [default: ./]
	-c    <num>     The cpu number. [default: 6]
	-q    <str>     queue:all,big,cu. [default: all]
	-p    <num>     Max thread number. [default: 5]
Example:
	perl 4.2.run_RNAseq_assessment.pl -i Unigene.fa -n A1,A2,A3 -s bowtie2_express -o bowtie2_express

USAGE

die $usage if(!$unigene || !$sample || !$samdir || !$outdir || $help);
$work_sh ||= getcwd;
$cpu ||= 6;
$queue ||= "all";
$thread_number ||= 5;
(-s $unigene) || die "Error: don't find unigene file: $unigene !\n";
$unigene=File::Spec->rel2abs($unigene);
(-d $samdir) || die "Error: don't find sam file directory: $samdir !\n";
$samdir=File::Spec->rel2abs($samdir);$samdir=~s/\/$//g;
(-d $outdir) || mkdir $outdir;
$outdir=File::Spec->rel2abs($outdir);$outdir=~s/\/$//g;
(-d $work_sh) || mkdir $work_sh;
$work_sh=File::Spec->rel2abs($work_sh);$work_sh=~s/\/$//g;

open SH1,">$work_sh/a.createBEDfromUnigene.sh" || die $!;
print SH1 "perl $Bin/createBEDfromUnigene.pl -i $unigene -o $outdir/Unigene.bed\n";
close SH1;

open SH2,">$work_sh/b.RNAseq_assessment.sh" || die $!;
my @sample=split /,/,$sample;
for my $s (@sample){
	print SH2 "perl $Bin/Saturation_assessment.pl --bam $samdir/$s.hits.sam --bed $outdir/Unigene.bed --sample $s --out $outdir\n";
	print SH2 "perl $Bin/geneBody_coverage.pl --sam $samdir/$s.hits.sam --bed $outdir/Unigene.bed --sample $s --out $outdir\n";
}
close SH2;

system("perl $qsub_pbs --queue $queue --lines 1 --maxjob $thread_number --num_proc $cpu $work_sh/a.createBEDfromUnigene.sh");
system("perl $qsub_pbs --queue $queue --lines 2 --maxjob $thread_number --num_proc $cpu $work_sh/b.RNAseq_assessment.sh");
