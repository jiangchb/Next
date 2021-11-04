#!/usr/bin/env perl -w
#Write by xiaoyue.wang 2018.02.05
use strict;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin);
use File::Basename;

############################# Initialization of software ################################
my $samtools="/data/software/samtools/samtools-v1.3.1/bin/samtools";
my $python = "/data/software/Anaconda3/bin/python";
my $geneBody_coverage = "/home/fanyucai/software/python/Python-v2.7.9/bin/geneBody_coverage.py";
my $env="export PATH=/data/software/gcc/gcc-v6.4.0//bin/:/data/software/Anaconda3/envs/mro-v3.4.3/bin/:\$PATH";
$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
#########################################################################################

### Usage
my $usage =<< "USAGE";
Program: $0
Description: Calculate the RNA-seq reads coverage over gene body.
Version: 1.0
=============================================================================
Options:
		-i  --sam     input sam file
		-r  --bed     reference bed file
		-s  --sample  sample name
		-o  --out     output file directory
		-h  --help    print this help info
=============================================================================

USAGE

my ($sam,$bed,$sample,$outdir,$help);
GetOptions(
	"sam|i=s"    => \$sam,
	"bed|r=s"    => \$bed,
	"sample|s=s" => \$sample,
	"out|o=s"    => \$outdir,
	"help|h"     => \$help,
);

if( !defined($sam) || !defined($bed) || !defined($sample) || !defined($outdir) || $help ){
	die $usage;
}
(-s $sam) || die "Error: don't open sam file: $sam !\n";
(-s $bed) || die "Error: don't open bed file: $bed !\n";
(-d $outdir) || mkdir $outdir;
$sam=File::Spec->rel2abs($sam);
$bed=File::Spec->rel2abs($bed);
$outdir=File::Spec->rel2abs($outdir);$outdir =~ s/\/$//;

my $bam_name=basename($sam);$bam_name=~s/\.sam$//;
$bam_name="$outdir/$bam_name";
system("$samtools view -S $sam -b -o $bam_name.bam");
system("$samtools sort -m 2G -o $bam_name.sort.bam $bam_name.bam && rm $bam_name.bam");
system("$samtools index $bam_name.sort.bam");
system("$python $geneBody_coverage -i $bam_name.sort.bam -r $bed -o $outdir/$sample");
system("$env && Rscript $Bin/draw_geneBody_coverage.r -i $outdir/$sample.geneBodyCoverage.txt -o $outdir -s $sample && rm $outdir/$sample.geneBodyCoverage.curves.pdf $outdir/$sample.geneBodyCoverage.r $bam_name.sort.bam $bam_name.sort.bam.bai");
