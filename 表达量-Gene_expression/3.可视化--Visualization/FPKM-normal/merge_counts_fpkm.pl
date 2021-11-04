#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use Math::BigFloat;
use File::Basename;

my ($indir, $outdir, $help);
GetOptions(
	"i:s"       => \$indir,
	"o:s"       => \$outdir,
	"h|help!"   => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: merge uniq_counts.fpkm.*.txt to counts.xls,fpkm.xls
Options:
	-i     <indir>    The input directory of uniq_counts.fpkm.*.txt
	-o     <outdir>   The output directory of result(counts.xls,fpkm.xls)
	-h|help           print help info
Example:
	perl merge_counts_fpkm.pl -i bowtie -o outdir/

USAGE

die $usage if(!$indir || !$outdir || $help);
(-d $indir) || die "Error: don't find indir: $indir !\n";
(-d $outdir) || mkdir $outdir;
$indir =File::Spec->rel2abs($indir);
$outdir=File::Spec->rel2abs($outdir);

my @infile=glob "$indir/uniq_counts.fpkm.*.txt";
(@infile==0) && die "Error: don't find uniq_counts.fpkm.*.txt in $indir !\n";

my %counts_fpkm;my @sample_sort;
for my $file(@infile){
	my $sample_name=basename($file);
	$sample_name=~s/^uniq_counts\.fpkm\.//g;
	$sample_name=~s/\.txt$//g;
	push(@sample_sort,$sample_name);
	(-s $file) || die "Error: don't find file: $file !\n";
	open IN,"<$file" || die $!;
	while(<IN>){
		chomp;
		next if(/^target_id/ || /^#/ || /^\s*$/);
		my @l=split /\t/;
		$counts_fpkm{$l[0]}{$sample_name}{counts}=$l[1];
		$counts_fpkm{$l[0]}{$sample_name}{fpkm}=new Math::BigFloat $l[2];
	}
	close IN;
}

open CS,">$outdir/counts.xls" || die $!;
open FM,">$outdir/fpkm.xls" || die $!;
print CS "gene_id\t".join("\t",@sample_sort)."\n";
print FM "gene_id\t".join("\t",@sample_sort)."\n";
for my $k1 (sort keys %counts_fpkm){
	print CS "$k1";
	print FM "$k1";
	for my $k2 (@sample_sort){
		print CS "\t$counts_fpkm{$k1}{$k2}{counts}";
		print FM "\t$counts_fpkm{$k1}{$k2}{fpkm}";
	}
	print CS "\n";
	print FM "\n";
}
close CS;
close FM;
