#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# getmaxlodsite.pl
#

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $usage = "usage:  getmaxlodsite.pl -lod lodfile -i seqfile -o outfile\n";

my $infile = "";
my $outfile = "";
my $lodfile = "";
my $VERB=0;
my $motlen = 0;
my $motif;

# Parse Input Arguments
for (my $i=0; $i<=$#ARGV; $i++) {
    
    my $arg = $ARGV[$i];
    if ($arg eq "-i") {$infile  = $ARGV[($i+1)];}
    if ($arg eq "-o") {$outfile = $ARGV[($i+1)];}
    if ($arg eq "-lod") {$lodfile = $ARGV[($i+1)];}
}
die "Input Error\!\n$usage" if ($outfile eq "" || $infile eq "" || $lodfile eq "");

# Read Log-Odds Matrix File (lodfile) for Motif
open(IN,"<$lodfile") or die "cannot open lodfile $lodfile\n";
while(<IN>) {
  
  chomp;
  if ($_ =~ /^\#/ || /ALPHABET/ || /log-odds matrix: alength/ || /CUTOFF/) {}
  #elsif ($_=~ /(\d+\.?\d*)\t(\d+\.?\d*)\t(\d+\.?\d*)\t(\d+\.?\d*)/) {
  
  else {
    
    my ($a, $c, $g, $t) = split /\t/, $_;
    
    if ($VERB==1) {print STDERR "Motif $motlen: A $a\tC $c\tG $g\tT $t\n";}
    
    $motif->{"A"}->[$motlen] = $a;
    $motif->{"C"}->[$motlen] = $c;
    $motif->{"G"}->[$motlen] = $g;
    $motif->{"T"}->[$motlen] = $t;
    $motlen++;
  }
}
if ($motlen>0) {
  
  # Read Input File & Find Max-Score Site For Each Sequence
  open(OUT,">$outfile") or die "cannot write to outfile $outfile\n";
  my $seqin  = Bio::SeqIO->new( -file => "<$infile");
  while (my $seqobj = $seqin->next_seq) {
    
    my $seq = uc($seqobj->seq);
    my $seqlen = length($seq);
    
    my ($maxscore, $maxi, $maxword) = (-1000, -1, "maxword");
    for (my $i=0; $i<=$seqlen-$motlen; $i++) {
      
      my $word = substr($seq, $i, $motlen);
      my $score = getlodscore($motif, $word, $motlen);
      if ($score>$maxscore) {
	$maxscore = $score;
	$maxi     = $i;
	$maxword  = $word;
      }
    }
    printf OUT "%s\t%1.6f\t%d\t%d\t%s\n", ($seqobj->display_name, $maxscore, ($maxi+1), ($maxi+$motlen), $maxword);
  }
  close(OUT);
}


sub getlodscore
  {
    my ($mot, $wd, $k) = @_;
    my $s = 0;
    if (length($wd) != $k) {die "getlodscore word length error\n";}
    for (my $i=0; $i<$k; $i++) {
      my $char = substr($wd, $i, 1);
      $s += $mot->{$char}->[$i];
    }
    return $s;
  }
