#!/usr/bin/env perl -w
#Modify by xiaoyue.wang 2018.02.05
use strict;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin);

############################# Initialization of software ################################
my $python = "/data/software/Anaconda3/bin/python";
my $RPKM_saturation = "/home/fanyucai/software/python/Python-v2.7.9/bin/RPKM_saturation.py";
my $env="export PATH=/data/software/gcc/gcc-v6.4.0//bin/:/data/software/Anaconda3/envs/mro-v3.4.3/bin/:\$PATH";
$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
#########################################################################################

### Usage
sub Usage {
  print STDERR  << "USAGE";
Description: RNAseq assessment--sequencing saturation evaluation
Version: 1.0
=============================================================================
Options:
		-b  --bam     input bam or sam file
		-r  --bed     reference bed file
		-s  --sample  sample name
		-o  --out     output file directory
		-h  --help    print this help info
=============================================================================

USAGE
  exit();
}

my ($bam,$bed,$sample,$outdir,$help);
GetOptions(
	"bam|b=s"    => \$bam,
	"bed|r=s"    => \$bed,
	"sample|s=s" => \$sample,
	"out|o=s"    => \$outdir,
	"help|h"     => \$help,
);

if( !defined($bam) || !defined($bed) || !defined($sample) || !defined($outdir) || $help ){
	&Usage();
}
(-s $bam) || die "Error: don't open bam or sam file: $bam !\n";
(-s $bed) || die "Error: don't open bed file: $bed !\n";
(-d $outdir) || mkdir $outdir;
$bam=File::Spec->rel2abs($bam);
$bed=File::Spec->rel2abs($bed);
$outdir=File::Spec->rel2abs($outdir);$outdir =~ s/\/$//;

system("$python $RPKM_saturation -i $bam -r $bed -o $outdir/$sample && $env && Rscript $Bin/run_Saturation.r -i $outdir/$sample.eRPKM.xls -o $outdir -s $sample");
system("rm $outdir/$sample.saturation.pdf $outdir/$sample.saturation.r");
