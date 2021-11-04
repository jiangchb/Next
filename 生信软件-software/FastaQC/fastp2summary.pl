#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);

################################### software/db/env #######################
my $fastp="/home/fanyucai/software/fastp/fastp";
my $env="export PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/bin/:/home/fanyucai/software/R/R-v3.4.0/bin/:\$PATH";
$env.=" && export LD_LIBRARY_PATH=/home/fanyucai/software/gcc/gcc-v6.1.0/lib64/:\$LD_LIBRARY_PATH";
###########################################################################

my ($raw1, $raw2, $clean1, $clean2, $odir, $sample, $help);
GetOptions(
	"r1:s"     =>\$raw1,
	"r2:s"     =>\$raw2,
	"c1:s"     =>\$clean1,
	"c2:s"     =>\$clean2,
	"s:s"      =>\$sample,
	"o:s"      =>\$odir,
	"help|h!"  =>\$help,
);

my $usage =<<"USAGE";
Program: $0 v1.0
Contact: xiaoyue.wang\@oebiotech.com
Description: Statistical summary file and draw qual,base picture from fastp
Options:
	-r1     <infile>     The input rawdata read1 fastq(.gz)
	-r2     <infile>     The input rawdata read2 fastq(.gz)
	-c1     <infile>     The input cleandata read1 fastq(.gz)
	-c2     <infile>     The input cleandata read2 fastq(.gz)
	-s      <name>       The sample name
	-o      <outdir>     The output directory of result
	-h|help              Print help info
Example:
	perl fastp2summary.pl -r1 AA.rawdata.R1.gz -r2 AA.rawdata.R2.gz -c1 AA.cleandata.R1.gz -c2 AA.cleandata.R2.gz -s AA -o result/

USAGE

die $usage if(!$raw1 || !$raw2 || !$clean1 || !$clean2 || !$sample || !$odir || $help);
(-s $raw1) || die "Error: don't find open r1: $raw1 !\n";
(-s $raw2) || die "Error: don't find open r2: $raw2 !\n";
(-s $clean1) || die "Error: don't find open c1: $clean1 !\n";
(-s $clean2) || die "Error: don't find open c2: $clean2 !\n";
(-d $odir) || system("mkdir -p $odir");
$raw1=File::Spec->rel2abs($raw1);
$raw2=File::Spec->rel2abs($raw2);
$clean1=File::Spec->rel2abs($clean1);
$clean2=File::Spec->rel2abs($clean2);
$odir=File::Spec->rel2abs($odir);

system("$fastp -w 15 -i $raw1 -I $raw2 --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --json $odir/$sample.rawdata.fastp.json --html $odir/$sample.rawdata.fastp.html");
system("$fastp -w 15 -i $clean1 -I $clean2 --disable_adapter_trimming --disable_trim_poly_g --disable_quality_filtering --disable_length_filtering --json $odir/$sample.cleandata.fastp.json --html $odir/$sample.cleandata.fastp.html");

my @infile=("$odir/$sample.rawdata.fastp.json", "$odir/$sample.cleandata.fastp.json");
my @head=("column", "mean_qual", "A_count", "T_count", "C_count", "G_count", "N_count", "GC_count");
$/="\n\t\},\n";
my %summary; my %sum; 
for my $f (@infile){
	my $mark=$1 if($f=~/\.(raw|clean)data\.fastp\.json$/);
	open JSON,"<$f" || die $!;
	while(<JSON>){
		chomp;
		if(/\"summary\": \{/){
			my @b=split /\n\t\t\},\n/,$_;
			for my $t (@b){
				if($t=~/\"before_filtering\": \{/){
					if($t=~/\"total_reads\":(\d+),/){
						my $val=$1;
						$summary{"$mark\_reads"}=$val;
						$sum{"$mark\_reads"}=sprintf("%.2f",$val/1000000)."M";
					}
					if($t=~/\"total_bases\":(\d+),/){
						my $val=$1;
						$summary{"$mark\_bases"}=$val;
						$sum{"$mark\_bases"}=sprintf("%.2f",$val/1000000000)."G";
					}
					if($t=~/\"q30_rate\":([\.\d]+),/ && $mark eq "raw"){
						my $val=$1;
						$summary{"Q30"}=sprintf("%.2f",$val*100);
						$sum{"Q30"}=sprintf("%.2f",$val*100);
					}
					if($t=~/\"gc_content\":([\.\d]+)/ && $mark eq "clean"){
						my $val=$1;
						$summary{"GC"}=sprintf("%.2f",$val*100);
						$sum{"GC"}=sprintf("%.2f",$val*100);
					}
				}else{
					next;
				}
			}
		}elsif(/\"read(1|2)_before_filtering\": \{/){
			my $readid=$1; my %stats;
			my @b=split /\n\t\t\},\n/,$_;
			for my $t (@b){
				if($t=~/\"quality_curves\": \{/){
					my $mean_qual=$1 if($t=~/\"mean\":\[(\S+)\]/);
					my @qual_val=split /,/,$mean_qual;
					for (my $i=0; $i<=$#qual_val; $i++){
						$stats{$i+1}{"mean_qual"}=$qual_val[$i];
					}
				}elsif($t=~/\"content_curves\": \{/){
					while($t=~/\"([ATCGN]+)\":\[(\S+)\]/g){
						my $base=$1; my $base_count=$2;
						my @count=split /,/,$base_count;
						for (my $i=0; $i<=$#count; $i++){
							$stats{$i+1}{"$base\_count"}=$count[$i]*100;
						}
					}
				}else{
					next;
				}
			}
			open STAT,">$odir/$sample\.$mark"."data.R$readid.stats.xls" || die $!;
			print STAT join("\t",@head)."\n";
			for my $k (sort {$a <=> $b} keys %stats){
				print STAT "$k";
				for (my $j=1; $j<=$#head; $j++){
					print STAT "\t$stats{$k}{$head[$j]}";
				}
				print STAT "\n";
			}
			close STAT;
			undef %stats;
		}else{
			next;
		}
	}
	close JSON;
}
undef @head;
$/="\n";

@head=("Sample","raw_reads","raw_bases","clean_reads","clean_bases","valid_bases","Q30","GC");
open SUM,">$odir/$sample.summary.xls" || die $!;
print SUM join("\t",@head)."\n";
my $valid_bases=sprintf("%.2f",$summary{clean_bases}*100/$summary{raw_bases});
print SUM "$sample\t$summary{raw_reads}\t$summary{raw_bases}\t$summary{clean_reads}\t$summary{clean_bases}\t$valid_bases%\t$summary{Q30}%\t$summary{GC}%\n";
print SUM "$sample(report)\t$sum{raw_reads}\t$sum{raw_bases}\t$sum{clean_reads}\t$sum{clean_bases}\t$valid_bases%\t$sum{Q30}%\t$sum{GC}%\n";
close SUM;

print "draw rawdata base/qual picture...\n";
system("$env && Rscript $Bin/ngsqc.r --fq1 $odir/$sample.rawdata.R1.stats.xls --fq2 $odir/$sample.rawdata.R2.stats.xls --len 150 --name $sample --key rawdata --od $odir");
print "draw cleandata base/qual picture...\n";
system("$env && Rscript $Bin/ngsqc.r --fq1 $odir/$sample.cleandata.R1.stats.xls --fq2 $odir/$sample.cleandata.R2.stats.xls --len 150 --name $sample --key cleandata --od $odir");
system("rm $odir/$sample.rawdata.fastp.json $odir/$sample.rawdata.fastp.html $odir/$sample.cleandata.fastp.json $odir/$sample.cleandata.fastp.html");
