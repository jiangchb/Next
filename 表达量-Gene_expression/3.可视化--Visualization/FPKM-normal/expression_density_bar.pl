#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use FindBin qw($Bin);
use Cwd;

###################### software or env ###############################################################
my $env="export PATH=/data/software/gcc/gcc-v6.4.0//bin/:/data/software/Anaconda3/envs/mro-v3.4.3/bin/:\$PATH";
$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
######################################################################################################

my ($infile, $express_type, $group_file, $outdir, $help);
GetOptions(
	"i:s"       => \$infile,
	"t:s"       => \$express_type,
	"g:s"       => \$group_file,
	"o:s"       => \$outdir,
	"h|help!"   => \$help,
);

my $usage=<< "USAGE";
Program: $0
Description: draw expression density picture and expression region picture
Options:
	-i     <file>     The infile(fpkm.xls)             [Required]
	-t     <str>      fpkm or rpm   [default: fpkm]    [Optional]
	-g     <file>     The infile(sample_group.xls)     [Optional]
	                  format: 2 columns
	                  Sample	Group
	                  A1	A
	                  A2	A
	                  B1	B
	-o     <outdir>   The output directory of result   [Required]
	-h|help           print help info
Example:
	perl expression_density_bar.pl -i fpkm.xls -o outdir/
USAGE

die $usage if(!$infile || !$outdir || $help);
(-s $infile) || die "Error: don't find infile: $infile !\n";
(-d $outdir) || mkdir $outdir;
$infile=File::Spec->rel2abs($infile);
$outdir=File::Spec->rel2abs($outdir);
if(defined $group_file){
	(-s $group_file) || die "Error: don't find infile: $group_file !\n";
	$group_file=File::Spec->rel2abs($group_file);
}
$express_type ||= "fpkm";

my %express_region;
open IN,"<$infile" || die $!;
my $head=<IN>;chomp $head;
my @head=split /\t/,$head;
open FD,">$outdir/$express_type\_density.xls" || die $!;
print FD "expression\tSample\n";
while(<IN>){
	chomp;
	my @l=split /\t/;
	for (my $i=1; $i<=$#l; $i++){
		my $express=$l[$i];
		print FD "$express\t$head[$i]\n" if($express>10**-40);
		if($express_type eq "fpkm"){
			if($express < 0.5){
				$express_region{$head[$i]}{"FPKM 0-0.5"}++;
			}elsif($express >= 0.5 && $express < 1){
				$express_region{$head[$i]}{"FPKM 0.5-1"}++;
			}elsif($express >= 1 && $express < 10){
				$express_region{$head[$i]}{"FPKM 1-10"}++;
			}else{
				$express_region{$head[$i]}{"FPKM >=10"}++;
			}
		}
		if($express_type eq "rpm"){
			if($express < 0.25){
				$express_region{$head[$i]}{"RPM 0-0.25"}++;
			}elsif($express >= 0.25 && $express < 1.5){
				$express_region{$head[$i]}{"RPM 0.25-1.5"}++;
			}elsif($express >= 1.5 && $express < 2.5){
				$express_region{$head[$i]}{"RPM 1.5-2.5"}++;
			}elsif($express >= 2.5){
				$express_region{$head[$i]}{"RPM >=2.5"}++;
			}
		}
	}
}
close IN;
close FD;

shift(@head);
my @region;
if($express_type eq "fpkm"){
	@region=("FPKM 0-0.5", "FPKM 0.5-1", "FPKM 1-10", "FPKM >=10");
}
if($express_type eq "rpm"){
	@region=("RPM 0-0.25", "RPM 0.25-1.5", "RPM 1.5-2.5", "RPM >=2.5");
}
open FR,">$outdir/$express_type\_region.xls" || die $!;
print FR "Sample\texpression_region\tnumber\n";
for my $s(@head){
	for my $r (@region){
		print FR "$s\t$r\t$express_region{$s}{$r}\n" if(exists $express_region{$s}{$r});
	}
}
close FR;

my $pwd=getcwd;
chdir $outdir;
if(defined $group_file){
	system("$env && Rscript $Bin/expression_density.r -i $outdir/$express_type\_density.xls -u $express_type -d $group_file && rm $outdir/$express_type\_density.xls");
	system("$env && Rscript $Bin/expression_region.r -i $outdir/$express_type\_region.xls -u $express_type -d $group_file -w 7 -t 6 && rm $outdir/$express_type\_region.xls");
}else{
	system("$env && Rscript $Bin/expression_density.r -i $outdir/$express_type\_density.xls -u $express_type && rm $outdir/$express_type\_density.xls");
	system("$env && Rscript $Bin/expression_region.r -i $outdir/$express_type\_region.xls -u $express_type -w 7 -t 6 && rm $outdir/$express_type\_region.xls");
}
chdir $pwd;
