#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use File::Basename;

my ($idir, $odir, $help);
GetOptions(
	"i:s"     =>\$idir,
	"o:s"     =>\$odir,
	"help|h!" =>\$help,
);

my $usage =<<"USAGE";
Program: $0 v1.0
Contact: xiaoyue.wang\@oebiotech.com
Description: generate Filtered_data_stat.xls and report_Filtered_data_stat.xls
Options:
	-i      <indir>      The input directory of sample.summary.xls file
	-o      <outdir>     The output directory of result
	-h|help              Print help info
Example:
	perl merge_fastp_summary.pl -i cleandata/ -o outdir/

USAGE

die $usage if(!$idir || !$odir || $help);
(-d $idir) || die "Error: don't find open indir: $idir !\n";
(-d $odir) || mkdir $odir;
$idir=File::Spec->rel2abs($idir);
$odir=File::Spec->rel2abs($odir);

my @infile=glob("$idir/*\.summary.xls");
(@infile==0) && die "Error: don't find *.summary.xls file in directory $idir !\n";

my @head=("Sample","RawReads","RawBases","CleanReads","CleanBases","ValidBases","Q30","GC");
open SUM,">$odir/Filtered_data_stat.xls" || die $!;
open ROP,">$odir/report_Filtered_data_stat.xls" || die $!;
print SUM join("\t",@head)."\n";
print ROP join("\t",@head)."\n";
for my $f (@infile){
	open IN,"<$f" || die $!;
	while(<IN>){
		chomp;
		next if(/^#/ || /^\s*$/ || /^Sample\tRawReads/ ||/^Sample\traw_reads/);
		if(/\(report\)/){
			s/\(report\)//g;
			print ROP "$_\n";
		}else{
			print SUM "$_\n";
		}
	}
	close IN;
}
close SUM;
close ROP;
