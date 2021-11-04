#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;

my ($infile, $outfile, $sample, $help);
GetOptions(
	"i:s"     => \$infile,
	"o:s"     => \$outfile,
	"s:s"     => \$sample,
	"h|help!" => \$help,
);

my $usage=<< "USAGE";

Program: $0
Description: merge and generate mapping_stat_transcriptome.xls
Options:
	-i     <str>  The input file. eg:Sample_A1.bowtie2.log,Sample_A2.bowtie2.log,Sample_A3.bowtie2.log
	-s     <str>  The sample name,One to one correspondence with the Options "-i". eg:Sample_A1,Sample_A2,Sample_A3
	-o     <str>  The output file name
	-h|help       print help info
Example:
	perl mapping_stat.pl -i Sample_A1.bowtie2.log,Sample_A2.bowtie2.log -s Sample_A1,Sample_A2 -o mapping_stat_transcriptome.xls

USAGE

die $usage if(!$infile || !$outfile || !$sample || $help);

my @infile_arr=split /,/,$infile;
my @sample_arr=split /,/,$sample;

my %map_stat;
for(my $i=0; $i<=$#sample_arr; $i++){
	$infile_arr[$i]=File::Spec->rel2abs($infile_arr[$i]);
	(-s $infile_arr[$i]) || die "Error: don't find file: $infile_arr[$i] !\n";
	open IN,"<$infile_arr[$i]" || die $!;
	my ($tmp1, $tmp2, $tmp3, $tmp4, $tmp5, $tmp6);
	while(<IN>){
		chomp;
		if(/^(\d+)\s+reads; of these:/){
			$tmp1=$1;
		}elsif(/^\s+(\d+)\s+\(\S+\)\s+aligned concordantly exactly 1 time/){
			$tmp2=$1;
		}elsif(/^\s+(\d+)\s+\(\S+\)\s+aligned concordantly >1 times/){
			$tmp3=$1;
		}elsif(/^\s+(\d+)\s+\(\S+\)\s+aligned discordantly 1 time/){
			$tmp4=$1;
		}elsif(/^\s+(\d+)\s+\(\S+\)\s+aligned exactly 1 time/){
			$tmp5=$1;
		}elsif(/^\s+(\d+)\s+\(\S+\)\s+aligned >1 times/){
			$tmp6=$1;
		}
	}
	close IN;
	die "Error: run bowtie2 is err($sample_arr[$i]) !\n" if(!defined $tmp1 || !defined $tmp2 || !defined $tmp3 || !defined $tmp4 || !defined $tmp5 || !defined $tmp6);
	my $Total_reads=$tmp1*2;
	my $Total_mapped_reads=$tmp2*2+$tmp3*2+$tmp4*2+$tmp5+$tmp6;
	my $Multiple_mapped=$tmp3*2+$tmp6;
	my $Uniquely_mapped=$tmp2*2+$tmp4*2+$tmp5;
	my $Reads_mapped_in_proper_pairs=$tmp2*2+$tmp3*2;

	my $per_Total_reads=sprintf("%.2f",($Total_reads*100)/$Total_reads);
	my $per_Total_mapped_reads=sprintf("%.2f",($Total_mapped_reads*100)/$Total_reads);
	my $per_Multiple_mapped=sprintf("%.2f",($Multiple_mapped*100)/$Total_reads);
	my $per_Uniquely_mapped=sprintf("%.2f",($Uniquely_mapped*100)/$Total_reads);
	my $per_Reads_mapped_in_proper_pairs=sprintf("%.2f",($Reads_mapped_in_proper_pairs*100)/$Total_reads);

	$map_stat{"Total reads"}{$sample_arr[$i]}="$Total_reads($per_Total_reads%)";
	$map_stat{"Total mapped reads"}{$sample_arr[$i]}="$Total_mapped_reads($per_Total_mapped_reads%)";
	$map_stat{"Multiple mapped"}{$sample_arr[$i]}="$Multiple_mapped($per_Multiple_mapped%)";
	$map_stat{"Uniquely mapped"}{$sample_arr[$i]}="$Uniquely_mapped($per_Uniquely_mapped%)";
	$map_stat{"Reads mapped in proper pairs"}{$sample_arr[$i]}="$Reads_mapped_in_proper_pairs($per_Reads_mapped_in_proper_pairs%)";
}

my @keys=("Total reads", "Total mapped reads", "Multiple mapped", "Uniquely mapped", "Reads mapped in proper pairs");
open OUT,">$outfile" || die $!;
print OUT "Term\\Sample\t".join("\t",@sample_arr)."\n";
for my $k1 (@keys){
	print OUT "$k1";
	for my $k2 (@sample_arr){
		print OUT "\t$map_stat{$k1}{$k2}";
	}
	print OUT "\n";
}
close OUT;
