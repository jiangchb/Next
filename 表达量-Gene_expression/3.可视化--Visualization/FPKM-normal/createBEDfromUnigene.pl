#!/usr/bin/env perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd;
my $cwd = getcwd();

### Usage
sub Usage {
  print STDERR << "USAGE";
Description:
Author:
Version:
=========================================================================
Options:
		-i  --in   input file
		-o  --out  output file
		-h  --help print this help info

=========================================================================

USAGE

  exit();
}

my ($in,$out,$help);

GetOptions("in|i=s" => \$in,
	   "out|o=s" => \$out,
	   "help|h" => \$help);

if( !defined($in) || !defined($out) ||
    $help ) {
  &Usage();
}

### core

&createBEDfromUnigene($in,$out);


### subroutine for createBEDfromUnigene
sub createBEDfromUnigene {
  my $fasta = shift;
  my $bed = shift;

  my $seq = &readFASTA($fasta);  
  open BED,">$bed" || die "Cannot write to this file!$!";
  foreach my $i(0 .. $#$seq) {
    my $len = $seq->[$i]{'len'};
    my $name = $seq->[$i]{'head'};
    my $chromStart = 0;
    my $chromEnd = $len - 1;
    my $score = 0;
    my $strand = '+';
    my $thickStart = 0;
    my $thickEnd = 0;
    my $itemRGB = 0;
    my $blockCount = 1;
    my $blockSize = "$len," ;
    my $blockStarts = "0,";
    print BED "$name\t$chromStart\t$chromEnd\t$name\t$score\t$strand\t$thickStart\t$thickEnd\t$itemRGB\t$blockCount\t$blockSize\t$blockStarts\n";
  }
  

  close BED;
}

###read fasta
sub readFASTA {
  my $fasta = shift;
  my @sequence = ();
  my $cnt = -1;

  open FASTA,"<$fasta" || die "Cannot open this file!$!";
  while(<FASTA>) {
    chomp;
    if(/^>(\S+)/) {
      $cnt++;
      $sequence[$cnt]{'head'} = $1;
    }
    if(/^\w+/) {
      $sequence[$cnt]{'seq'} .= $_;
    }
  }

  close FASTA;

  foreach my $i(0 .. $cnt) {
    my $len = length $sequence[$i]{'seq'};
    $sequence[$i]{'len'} = $len;
  }

  return \@sequence;
}
