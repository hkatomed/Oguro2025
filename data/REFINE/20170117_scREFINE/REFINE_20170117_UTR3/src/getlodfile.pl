#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# getlodfile.pl
#

# This program takes a set of input sequences from infile and outputs 
# a log-odds scoring matrix to output (lodfile) based on a 
# user-specified background distribution.

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $usage = "usage:  getlodfile.pl -i infile [-bfile nt.freqs] [-base 10] [-pseud N]\n";

sub getCutoff {
  
  my ($min, $avg, $max) = @_;
  my $cut = ($avg+$min)/2;
  return $cut;
}

my $infile = "";
my $outfile = "";
my $lodfile = "";
my $bakfile = "";
my $VERB=0;
my $base = 10;
my %f;
$f{"A"} = 0.25;
$f{"C"} = 0.25;
$f{"G"} = 0.25;
$f{"T"} = 0.25;
my $motlen = 0;
my $motif;
my $PSEUD = 0;


# Parse Input Arguments
for (my $i=0; $i<=$#ARGV; $i++) {
  
  my $arg = $ARGV[$i];
  if ($arg eq "-i") {$infile  = $ARGV[($i+1)];}
  if ($arg eq "-o") {$outfile = $ARGV[($i+1)];}
  if ($arg eq "-lod") {$lodfile = $ARGV[($i+1)];}
  if ($arg eq "-bfile") {$bakfile = $ARGV[($i+1)];}
  if ($arg eq "-len") {$motlen  = $ARGV[($i+1)];}
  if ($arg eq "-a") {$f{"A"}  = $ARGV[($i+1)];}
  if ($arg eq "-c") {$f{"C"}  = $ARGV[($i+1)];}
  if ($arg eq "-g") {$f{"G"}  = $ARGV[($i+1)];}
  if ($arg eq "-t") {$f{"T"}  = $ARGV[($i+1)];}
  if ($arg eq "-pseud") {$PSEUD  = $ARGV[($i+1)];}
  if ($arg eq "-base") {$base  = $ARGV[($i+1)];}
}
die "Input Error\!\n$usage" if ($infile eq "" || $base<=1 || $motlen<0);

if ($bakfile ne "") {
  
  open(BAK,$bakfile) or die "cannot open bfile $bakfile\n";
  while(<BAK>) {
    chomp;
    if (m/([ACGT])+\t(\d+\.*\d*)/) {
      $f{uc($1)} = $2;
    }
  }
}

# Normalize Bkgrd Distribution
my $sum = $f{"A"} + $f{"C"} + $f{"G"} + $f{"T"};
$f{"A"} = $f{"A"}/$sum;
$f{"C"} = $f{"C"}/$sum;
$f{"G"} = $f{"G"}/$sum;
$f{"T"} = $f{"T"}/$sum;


# Initializations
my @sums;
my @seqs = ();
# Read Input Sequences
open(IN,"$infile") or die "cannot open infile $infile\n";
while(<IN>) {
  
  chomp;
  
  next if (m/^>/);
  if ($_=~ /([\.\-acgtuxnACGTUXN]+)/) {
    my $seq = uc($1);
    $seq =~ s/U/T/g;
    push @seqs, $seq;
    
    if ($motlen==0) {
      $motlen = length($seq); 
      for (my $i=0; $i<$motlen; $i++) {
	foreach my $char ("A", "C", "G", "T") {
	  $motif->{$char}->[$i] = $PSEUD;
	}
	$sums[$i]=4*$PSEUD;
      }
    }
    
    if (length($seq) != $motlen || $motlen<=0) {
      die "";
    }
    else {
      for (my $i=0; $i<$motlen; $i++) {
	my $char = substr($seq, $i, 1);	
	if ($char =~ /[ACGT]/) {
	  my $tmp = $motif->{$char}->[$i];
	  $motif->{$char}->[$i] = $tmp + 1;
	  $sums[$i] = $sums[$i]+1;
	}
      }
    }
  }
}
close(IN);

if ($motlen>0) {
  
  my @relent = ();
  
  # Normalize Motif Model freqs and Calc Log-Odds Scores
  for (my $i=0; $i<$motlen; $i++) {
    
    $relent[$i] = 0;
    foreach my $char ("A","C","G","T") {
      
      my $tmp = $motif->{$char}->[$i];
      if ($sums[$i]>4*$PSEUD) {$tmp = $tmp/$sums[$i];}
      else { die "FATAL ERROR:  abnormally low sum for col $i:  $sums[$i]\n";}
      
      my $tmp2 = $tmp/$f{$char};
      my $tmp3 = log($tmp2)/log($base);
      $motif->{$char}->[$i] = $tmp3;
      $relent[$i] += ($tmp*$tmp3);	
    }
  }
  
  # Calculate Cutoff Score
  my ($mincut, $maxcut, $avgcut, $numcut) = (0,0,0,0);
  foreach my $seq (@seqs) {
    unless ($seq =~ /x|n/i) {
      my $s = getlodscore($motif, $seq, $motlen);
      if ($mincut==0 || $mincut>$s) {$mincut = $s;}
      if ($maxcut<$s) {$maxcut = $s;}
      $avgcut += $s;
      $numcut++;
      #print STDERR "$numcut:\tlodscore of $seq = $s\n";
    }
  }
  $avgcut = ($avgcut/$numcut);
  my $cutoff = getCutoff($mincut, $avgcut, $maxcut);
  #printf STDERR "cut=%1.4f\tavg=%1.4f\tmin=%1.4f\tmax=%1.4f\n", ($cutoff, $avgcut, $mincut, $maxcut);
  
  # Print Motif Model to Output
  
  $cutoff = sprintf("%1.6f",$cutoff);
  
  print "CUTOFF= $cutoff\n";
  print "ALPHABET= ACGT\n";
  print "log-odds matrix: alength= 4 w= " . $motlen . "\n";
  
  #print STDERR "RelEnt\n";
  for (my $i=0; $i<$motlen; $i++) {
    foreach my $char ("A", "C", "G", "T") {
      printf "%1.6f\t", $motif->{$char}->[$i];
    }
    print "\n";
    #printf STDERR "$i\t%1.3f\n", $relent[$i];
  }
}

sub getlodscore {
  my ($mot, $wd, $k) = @_;
  my $s = 0;
  if (length($wd) != $k) {die "getlodscore word length error\n";}
  for (my $i=0; $i<$k; $i++) {
    my $char = substr($wd, $i, 1);
    $s += $mot->{$char}->[$i];
  }
  return $s;
}








