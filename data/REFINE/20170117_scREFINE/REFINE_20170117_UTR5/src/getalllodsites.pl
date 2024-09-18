#!/usr/bin/perl
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# getalllodsites.pl
#
# this program reads in a fasta seq library file and finds all non-overlapping
# sites with log-odds scores above the specified cutoff score.
# it prints data on the sites in the same format as getmaxlodsite.pl
# if the cutoff score is not specified as an input parameter, it
# is calculated as the midpoint between the avg score and the minimum score
# based on the sequences in the cutfile
#

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $usage = "usage:  getalllodsites.pl [-c C | -cut cutfile] -lod lodfile -i seqfile -o outfile\n";

my $infile  = "";
my $outfile = "";
my $lodfile = "";
my $cutfile = "";
my $cutoff  = 0;
my $mincut  = 0;
my $maxcut  = 0;
my $avgcut  = 0;
my $numcut  = 0;
my $VERB    = 0;
my $motlen  = 0;
my $motif;
my $LODCUT  = 1;  # 1 => read cutoff lod-score from lodfile matrix
my $flank   = 0;  # number of adjacent bases to include
my $OVERLAP = 0;  # 1 allows over-lapping motif sites

sub getCutoff {
  
  my ($min, $avg, $max) = @_;
  my $cut = ($avg+$min)/2;
  return $cut;
}

# Parse Input Arguments
for (my $i=0; $i<=$#ARGV; $i++) {
    
    my $arg = $ARGV[$i];
    if ($arg eq "-c") {$cutoff  = $ARGV[($i+1)]; $LODCUT=0;}
    if ($arg eq "-i")   {$infile  = $ARGV[($i+1)];}
    if ($arg eq "-o")   {$outfile = $ARGV[($i+1)];}
    if ($arg eq "-lod") {$lodfile = $ARGV[($i+1)];}
    if ($arg eq "-cut") {$cutfile = $ARGV[($i+1)]; $LODCUT=0;}
    if ($arg eq "-flank") {$flank = $ARGV[($i+1)];}
    if ($arg eq "-overlap") {$OVERLAP = 1;}
  }
die "Input Error\!\n$usage" if ($outfile eq "" || $infile eq "" || $lodfile eq "");

# Read Log-Odds Matrix File (lodfile) for Motif
open(IN,"<$lodfile") or die "cannot open lodfile $lodfile\n";
while(<IN>) {
  
  chomp;
  if (m/ALPHABET|log-odds matrix:/) {}
  elsif (m/CUTOFF= (.+)/ && $LODCUT==1) {
    $cutoff =  sprintf("%1.6f",$1);
    if ($VERB==1) {
      print STDERR "Cutoff lodscore = $cutoff from $lodfile\n";
    }
  }
  else {
    
    my ($a, $c, $g, $t) = split /\t/, $_;
    
    if ($VERB==1) {print STDERR "Motif $motlen: A $a\tC $c\tG $g\tT $t\n";}
    
    $motif->{"A"}->[$motlen] = $a;
    $motif->{"C"}->[$motlen] = $c;
    $motif->{"G"}->[$motlen] = $g;
    $motif->{"T"}->[$motlen] = $t;
    $motlen++;
  }
  #else {die "lodfile format error: $_\n";}
}
if ($motlen<=0) {die "motif length zero from $lodfile\n";}
close(IN);


# Read cutfile if specified
if ($cutfile ne "") {
  open(IN,"<$cutfile") or die "cannot open cutfile $cutfile\n\n";
  while(<IN>) {
    chomp;
    next if (m/^>/);
    unless ($_ =~ /x|n/i) {
      my $s = getlodscore($motif, $_, $motlen);
      #if ($cutoff==0 || $s<$cutoff) {$cutoff=$s;}
      if ($mincut==0 || $mincut>$s) {$mincut = $s;}
      if ($maxcut<$s) {$maxcut = $s;}
      $avgcut += $s;
      $numcut++;
      print STDERR "lodscore of $_ = $s\n";
    }
  }
  close(IN);
  $avgcut = ($avgcut/$numcut);
  $cutoff = getCutoff($mincut, $avgcut, $maxcut);
  printf "score cutoff = %1.3f (%1.3f, %1.3f-%1.3f) %s\n",
    ($cutoff, $avgcut, $mincut, $maxcut, $cutfile);
}

# Read Input File & Find Above-Cutoff-Score Sites For Each Sequence
open(OUT,">$outfile") or die "cannot write to outfile $outfile\n";
my $seqin  = Bio::SeqIO->new( -file => "<$infile");
while (my $seqobj = $seqin->next_seq) {
  
  my $seq = uc($seqobj->seq);
  my $seqlen = length($seq);
  
  my ($maxscore, $maxi, $maxword) = (-1000, -1, "maxword");
  for (my $i=0; $i<=$seqlen-$motlen; $i++) {
    
    my $word = substr($seq, $i, $motlen);
    my $score = getlodscore($motif, $word, $motlen);
    if ($score>=$cutoff) {
      
      my ($lflank, $rflank) = ("","");
      if ($flank>0) {
	if ($i-$flank>=0) {
	  $lflank = lc(substr($seq, $i-$flank, $flank));
	}
	else {$lflank = "x" x ($flank-$i) . lc(substr($seq, 0, $i));}
	
	if ($i+$flank+$motlen<=$seqlen) {
	  $rflank = lc(substr($seq, $i+$motlen, $flank));
	}
	else {$rflank = lc(substr($seq, $i+$motlen)) . "x" x ($i+$flank+$motlen-$seqlen);}
      }
      
      printf OUT "%s\t%1.6f\t%d\t%d\t%s\n", 
	($seqobj->display_name, $score, ($i+1), ($i+$motlen), ($lflank . $word . $rflank));
      if ($OVERLAP==0) {$i = $i + $motlen;}
    }
  }
}
close(OUT);

sub getlodscore {
  my ($mot, $wd, $k) = @_;
  my $s = 0;
  if (length($wd) != $k) {die "getlodscore word length error\n";}
  for (my $i=0; $i<$k; $i++) {
    my $char = substr($wd, $i, 1);
    $s += $mot->{$char}->[$i];
  }
  return sprintf("%1.6f",$s);
  #return $s;
}
