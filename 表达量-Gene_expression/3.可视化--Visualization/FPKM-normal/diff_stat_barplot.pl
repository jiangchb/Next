#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Spec;

############################# Initialization of software ################################
my $env="export PATH=/data/software/gcc/gcc-v6.4.0//bin/:/data/software/Anaconda3/envs/mro-v3.4.3/bin/:\$PATH";
$env.=" && export LD_LIBRARY_PATH=/data/software/gcc/gcc-v6.4.0/lib64/:\$LD_LIBRARY_PATH";
#########################################################################################

my ($tab,$outdir,$height,$width,$help);
GetOptions(
	"tab|i=s"    => \$tab,
	"out|o=s"    => \$outdir,
	"hgt|g=i"    => \$height,
	"wid|w=i"    => \$width,
	"help|h!"    => \$help,
);

my $usage =<< "USAGE";
Program: $0
Description: draw up,down diffgene number barplot
Version: 1.0
Options:
	-i  --tab      The input diff_stat.xls file
	-o  --out      The output file directory
	-g  --hgt      The pucture height, default=6
	-w  --wid      The pucture width, default=7
	-h  --help     print this help info
Example:
	perl diff_stat_barplot.pl -i diff_stat.xls -o Quantification/

USAGE

die $usage if(!$tab || !$outdir || $help);
$height ||= 6;
$width ||= 7;
(-s $tab) || die "Error: don't open infile: $tab !\n";
(-d $outdir) || mkdir $outdir;
$tab=File::Spec->rel2abs($tab);
$outdir=File::Spec->rel2abs($outdir);$outdir =~ s/\/$//;

open IN,"<$tab" || die $!;
open OUT,">$outdir/diff_stat.plot.xls" || die $!;
print OUT "Group\tType\tGene_number\n";
while(<IN>){
	chomp;
	next if(/^#/ || /^\s*$/);
	my @l=split /\t/;
	if($.==1){
		$l[0]=~/Case/i || die "Error: infile header is error!\n";
		$l[1]=~/Control/i || die "Error: infile header is error!\n";
		$l[2]=~/Up/i || die "Error: infile header is error!\n";
		$l[3]=~/Down/i || die "Error: infile header is error!\n";
	}else{
		my $case_name=$l[0]; my $control_name=$l[1];
		$case_name=$1 if($l[0]=~/\((\S+)\)/);
		$control_name=$1 if($l[1]=~/\((\S+)\)/);
		print OUT "$case_name\_vs_$control_name\tUp\t$l[2]\n$case_name\_vs_$control_name\tDown\t$l[3]\n";
	}
}
close IN;
close OUT;

open R,">$outdir/diff_stat_barplot.r" || die $!;
print R "#!/usr/bin/env Rscript
library(ggplot2)
data <- read.table(\"$outdir/diff_stat.plot.xls\", header=T, sep=\"\\t\", quote=\"\", check.names=F);
uniq<-data\$Group[!duplicated(data\$Group)]
data\$Group <- factor(data\$Group, levels=uniq)
data\$Type <- factor(data\$Type, levels=c(\"Up\", \"Down\"));
p=ggplot(data=data,aes(x=Group, y=Gene_number, fill=Type, width=0.7, space=0))+
geom_bar(stat=\"identity\", position=\"dodge\")+
geom_text(aes(label=data\$Gene_number),hjust=0.5,vjust=-0.5,size=2.7,position = position_dodge(0.7))+
xlab(\"\")+ylab(\"DEG number\")+
labs(title=\"Statistic of Differently Expressed Gene\") +
theme(plot.title = element_text(hjust = 0.5,size=13)) +
theme(panel.border=element_rect(fill=NA,colour=\"black\"),
panel.background = element_rect(fill=\"transparent\",colour=NA),
plot.background = element_rect(fill=\"transparent\",colour=NA))+
theme(axis.text.y=element_text(size=10,color=\"black\")) + 
theme(axis.text.x=element_text(angle = 45, hjust=1, vjust = 1, size=10,color=\"black\")) +
theme(legend.text=element_text(size=10))+
theme(panel.grid =element_blank())\n
ggsave(\"$outdir/gene_diff_stat_barplot.pdf\", height=$height, width=$width, plot=p)
ggsave(\"$outdir/gene_diff_stat_barplot.png\", type=\"cairo-png\", height=$height, width=$width, plot=p)\n";

system("$env && R --restore --no-save < $outdir/diff_stat_barplot.r && rm $outdir/diff_stat_barplot.r $outdir/diff_stat.plot.xls");
