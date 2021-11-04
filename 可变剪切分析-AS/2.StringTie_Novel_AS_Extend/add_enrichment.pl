#!/usr/bin/perl
my %Transcript_genesymbol;
###读如Transcript_genesymbol.txt到%Transcript_genesymbol
open(IN,"$ARGV[0]");#Transcript_genesymbol.txt来自annotation.txt
while(<IN>)
{
    chomp();
    my ($key,$val) = split(/\s+/,$_);
    $Transcript_genesymbol{$key}=$val;
}
#print %Transcript_genesymbol;
###读入enrichment-kegg-Sample_T.VS.Sample_C.xl
open(FILE,"$ARGV[1]");#KEGG或者GO结果文件
open (OUT,">$ARGV[2]");#结果文件
my $name = <FILE>;
print OUT $name;#输入标题
while(<FILE>){
		chomp;
		my @Genes = split(/\t/,$_);
		my @transcripts=split(/;/,$Genes[-1]);
		my $before=join("\t",@Genes[0..($#Genes-1)]);
		print OUT "$before\t";
		foreach my $tra (sort @transcripts){
			if(exists $Transcript_genesymbol{$tra}){
	        		print OUT "$tra($Transcript_genesymbol{$tra});";
			}else{
				print OUT "$tra();";
			}
	
		}
		print OUT "\n";
}

