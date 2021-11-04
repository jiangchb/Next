#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use File::Basename;

my %opts;
GetOptions(\%opts, "i=s", "a:s", "o=s", "help|h!");

my $usage=<< "USAGE";
Program : $0 
Usage : perl $0 [options]
	-i   <file>      The input file
	-a   <file>      The input annotation file
	-o   <dir >      The result output directory
Example :
	perl add_annotation.pl -i Quantification/counts.xls,Quantification/fpkm.xls -a annotation.xls -o Quantification/

USAGE

die $usage if(!$opts{i} || !$opts{a} || !$opts{o} || $opts{help});
my @infile=split /,/,$opts{i};
for(my $i=0; $i<@infile; $i++){
	(-s $infile[$i]) || die "Error: don't open $infile[$i] !\n";
	$infile[$i]=File::Spec->rel2abs($infile[$i]);
}

(-s $opts{a}) || die "Error: don't open file: $opts{a} !\n";
$opts{a}=File::Spec->rel2abs($opts{a});

(-d $opts{o}) || mkdir $opts{o};
$opts{o}=File::Spec->rel2abs($opts{o});

my %des;my $head;
open AN,"<$opts{a}" || die $!;
while(<AN>){
	chomp;
	next if(/^\s*$/);
	my @t=split /\t/;
	if($.==1){
		$head=join("\t",@t[1..$#t]);
	}else{
		$des{$t[0]}=join("\t",@t[1..$#t]);
	}
}
close AN;

for(my $i=0; $i<@infile; $i++){
	open IN,"<$infile[$i]" || die $!;
	my $filename=basename($infile[$i]); $filename=~s/\.\w+$//;
	open OUT,">$opts{o}/$filename\_anno.xls" || die $!;
	while(<IN>){
		chomp;
		next if(/^\s*$/);
		my @t=split /\t/,$_,2;
		if($.==1){
			print OUT "$_\t$head\n";
		}else{
			print OUT "$_\t$des{$t[0]}\n";
		}
	}
	close IN;
	close OUT;
}
