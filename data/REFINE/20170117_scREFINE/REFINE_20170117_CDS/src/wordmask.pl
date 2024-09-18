#!/usr/bin/perl
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# wordmask.pl
# 
# usage: wordmask.pl [-ss|-ds] [-x | -lc] [-wds wdlist.txt] seqfile.fa
#
# this program reads in a fasta sequence file and a list of identical-length
# k-mer words in the text file wdlist.txt.  it then finds all instances of
# words from the list that exist in the sequences, and masks all parts of the
# sequence that are not tiled by a word from the list as either lowercase (-lc)
# or X bases (-x).
#
#

use strict;
use Bio::SeqIO;

my $usage = "usage:\twordmask.pl [-ss|-ds] [-x|-lc] [-flank N] [-gap N] "
  . "[-ct N] [-max N] -wds wordlist -i seqs.fa -o outfile.fa\n\n";

my $K      =  0;   # word length for k-mers
my $DS     =  0;   # consider only fwd strand words
my $MASK   =  1;   # 1 => X-masking, 0 => lower-case masking
my $wdfile = "";   # filename for list of words when not using ALL k-mers
my $infile = "";   # filename for fasta sequences
my $outfile = "";  # filename for sequence output file
my %words  = ();   # hash for words in wdfile
my $GAPLEN =  0;   # minimum gap length for masking btwn listed words
my $VERB   =  0;   # verbose mode
my $flank  =  0;   # automatically unmask bases adjacent to listed words
my $COUNTFACT = 0; # sets limit for number of matching words req'd for inclusion
my $MAXMASK = 0;   # maximum length of x's in a row
my $CLIPX = 0;     # 1 => remove beginning and trailing x's after masking

#process input arguments
for (my $i=0; $i<$#ARGV; $i++) {
  my $arg = $ARGV[$i];
  if ($arg eq "-ss")   {$DS = 0;}
  if ($arg eq "-ds")   {$DS = 1;}
  if ($arg eq "-x")    {$MASK = 1;}
  if ($arg eq "-lc")   {$MASK = 0;}
  if ($arg eq "-i")    {$infile  = $ARGV[($i+1)];}
  if ($arg eq "-o")    {$outfile = $ARGV[($i+1)];}
  if ($arg eq "-wds")  {$wdfile  = $ARGV[($i+1)];}
  if ($arg eq "-gap")  {$GAPLEN  = $ARGV[($i+1)];}
  if ($arg eq "-flank")  {$flank  = $ARGV[($i+1)];}
  if ($arg eq "-ct")  {$COUNTFACT  = $ARGV[($i+1)];}
  if ($arg eq "-max") {$MAXMASK = $ARGV[($i+1)];}
}
if ($infile eq "" || $wdfile eq "" || $outfile eq "") {die "\n$usage";}

my $numWords = 0;
open(WDS,$wdfile) or die "cannot open $wdfile\n";
while (<WDS>) {
  chomp;
  if ($_=~/^([ACGT]+)$/) {
    $words{uc($1)} = 1;
    $numWords++;
    if ($DS>0) {
      my $rc = uc(reverse($1));
      $rc =~ tr/ACGT/TGCA/;
      $words{$rc} = 1;
      $numWords++;
    }
    if ($K>0 && length($1)>0 && length($1)!=$K) {
      die "word size error:\t$1\n";
    }
    $K = length($1);
  }
}
close(WDS);

# Process sequences
my $totlen = 0;  # number of total bases in all input seqs
my $unmlen = 0;  # number of unmaskes bases in all seqs
my $totseq = 0;
my $unmseq = 0;

my $out   = Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');
my $seqio = Bio::SeqIO->new(-file => "$infile");
while (my $seqobj = $seqio->next_seq) {
  
  my $name = $seqobj->display_name;
  my $seq = lc($seqobj->seq);
  my $len = length($seq);
  my $count = 0;
  my $minCount = ($COUNTFACT)*($numWords/pow(4,$K))*($len-$K+1);
  if ($minCount > ($len-$K+1)) {$minCount = 0;}
  for (my $i=0; $i<=$len-$K; $i++) {
    my $wd = uc(substr($seq, $i, $K));
    if (exists($words{$wd})) {
      
      my $left  = substr($seq, 0, $i);
      my $right = substr($seq, $i+$K);
      if ($flank>0) {
	$left  =~ s/(.{0,$flank}$)/uc($1)/e;
	$right =~ s/(^.{0,$flank})/uc($1)/e;
      }
      $seq = $left . $wd . $right;
      $count++;
    }
  }
  if ($GAPLEN>0)  {
    while ($seq =~ m/([A-Z][a-z]{1,$GAPLEN}[A-Z])/) {
      $seq =~ s/([A-Z][a-z]{1,$GAPLEN}[A-Z])/uc($1)/eg;
    }
  }
  if ($MASK==1) {$seq =~ s/[a-z]/x/g;}
  if ($len != length($seq)) {die "Length mismatch error\!\n\n";}
  if ($MASK==1 && $MAXMASK>0) {
    $seq =~ s/(x{$MAXMASK})x+/$1/g;
    if ($CLIPX>0) {
      $seq =~ s/^x+//;  # remove starting x's
      $seq =~ s/x+$//;  # remove trailing x's
    }
  }
  $seqobj->seq($seq);
  if ($COUNTFACT>0) {$name = "$name:" . $count; $seqobj->display_name($name);}
  if ($count>0 && $count>=$minCount) {
    $out->write_seq($seqobj);
    $unmseq++;
  }
  $totseq++;
  $totlen += length($seq);
  $unmlen += ($seq =~ tr/ACGT/ACGT/);
}
if ($VERB==1) {
  printf STDERR "%1.3f of %d bases masked\n", (1-($unmlen/$totlen)), $totlen;
  printf STDERR "$unmseq of $totseq sequences retained\n";
}

sub pow {
  my ($a, $b) = @_;
  my $x = 1;
  for (my $i=0; $i<$b; $i++) {
    $x *= $a;
  }
  return $x;
}
