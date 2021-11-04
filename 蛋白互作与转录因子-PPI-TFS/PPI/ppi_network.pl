#!/usr/bin/env perl
use strict;
use warnings;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;
use FindBin qw($Bin);

########################### software/env/database #####################################
#######################################################################################

my ($diff_indir, $genome_path,$outdir, $help);
GetOptions(
	"i:s"     => \$diff_indir,
	"s:s"       => \$genome_path,
	"o:s"     => \$outdir,
	"h|help!" => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description:
Options:
	-i      <indir>      The input directory of *-vs-*-diff-p*FC-*.xls               [Required]
	-s      <file>       genome_path  [Required]
	-o      <outdir>     The output directory of result                              [Required]
	-h|help              print help info
Example:
	perl ppi_network.pl -i diff/ -s Populus_trichocarpa_3694 -o out

USAGE

die $usage if(!$diff_indir || !$outdir || !$genome_path || $help);
(-d $diff_indir) || die "Error: don't find indir: $diff_indir !\n";
(-d $outdir) || mkdir $outdir;
$diff_indir=File::Spec->rel2abs($diff_indir);
$outdir=File::Spec->rel2abs($outdir);
#$string_species=File::Spec->rel2abs($string_species);
my $stringdb="$genome_path/gene2gene_network.xls";
(-s $stringdb) || die "Error: $genome_path is error, don't find string database: $stringdb !\n";

my @diff_file=glob "$diff_indir/*-vs-*-diff-p*FC-*.xls";
(@diff_file==0) && die "Error: don't find *-vs-*-diff-p*FC-*.xls in $diff_indir !\n";

for my $diff (@diff_file){
	(-s $diff) || die "Error: don't open file: $diff !\n";
	open DF,"<$diff" || die $!;
	my %diff_gene;
	while(<DF>){
		chomp;
		next if(/^#/ || /^\s*$/ || /\S+\tbaseMean/);
		my @l=split /\t/;
		$diff_gene{$l[0]}=$l[8];
	}
	close DF;

	my $group_name=basename($diff);
	$group_name=(split /-diff-p/,$group_name)[0];
	open OUT,">$outdir/diff-$group_name\_gene2gene_network.txt" || die $!;
	open TOP,">$outdir/top_300_diff-$group_name\_gene2gene_network.txt" || die $!;
	print OUT "gene1\tgene_type\tgene2\tgene_type\tcombined_score\n";
	print TOP "gene1\tgene_type\tgene2\tgene_type\tcombined_score\n";
	my $n=1; my %nodes_num; my %nodes_cfg;
	open IN,"<$stringdb" || die $!;
	while(<IN>){
		chomp;
		next if(/^#/ || /^\s*$/ || /combined_score/);
		my @l=split /\t/;
		if(exists $diff_gene{$l[0]} && exists $diff_gene{$l[1]}){
			print OUT "$l[0]\t$diff_gene{$l[0]}\t$l[1]\t$diff_gene{$l[1]}\t$l[2]\n";
			if($n <= 300){
				print TOP "$l[0]\t$diff_gene{$l[0]}\t$l[1]\t$diff_gene{$l[1]}\t$l[2]\n";
				$nodes_num{$l[0]}++; $nodes_num{$l[1]}++;
				my $col1;my $col2;
				($diff_gene{$l[0]}=~/Up/i) ? ($col1="red") : ($col1="green");
				($diff_gene{$l[1]}=~/Up/i) ? ($col2="red") : ($col2="green");
				$nodes_cfg{$l[0]}="$col1\t$diff_gene{$l[0]}";
				$nodes_cfg{$l[1]}="$col2\t$diff_gene{$l[1]}";
			}
			$n++;
		}
	}
	close IN;
	close OUT;
	close TOP;
}
