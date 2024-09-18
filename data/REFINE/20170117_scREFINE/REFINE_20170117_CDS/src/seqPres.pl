#!/usr/bin/perl
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# seqPres.pl
#
# usage: seqPres.pl [-ss|-ds] [-k k | -size k] [-wds wdlist.txt] seqfile
#
#

use strict;
use Bio::SeqIO;

my $usage = "seqCons.pl [-ss | -ds] [-size k] [-wds wdlist.txt] seqfile\n\n";
my $K      =  6;   # word length for k-mers
my $DS     =  0;   # consider only fwd strand words
my $CTS    =  1;   # report counts (0 => frequencies)
my $total  =  0;   # total number of k-mers in sequence set
my %freq   =  ();  # hash of k-mer frequencies/counts
my %pwfreq = ();   # hash for per-species full-conservation counts
my $wdfile = "";   # filename for list of words when not using ALL k-mers
my $maxnseqs = 7;
my $total  =  0;
my @alnbeg = ();
my @alnend = ();
my $GAPPCT = .5;

#process input arguments
for (my $i=0; $i<$#ARGV; $i++) {
  my $arg = $ARGV[$i];
  if ($arg eq "-size") {$K  = $ARGV[($i+1)];}
  if ($arg eq "-ss") {$DS = 0;}
  if ($arg eq "-ds") {$DS = 1;}
  if ($arg eq "-cts") {$CTS = 1;}
  if ($arg eq "-nocts") {$CTS = 0;}
  if ($arg eq "-wds")  {$wdfile = $ARGV[($i+1)];}
  if ($arg eq "-maxn") {$maxnseqs = $ARGV[($i+1)];}
}
my $infile = $ARGV[$#ARGV];
if ($infile eq "" || $K<1) {die "$usage";}

my @Kmers = initKmers($wdfile, $K);
foreach my $wd (@Kmers) {
  
  $freq{$wd}->{"con"} = 0; # 
  $freq{$wd}->{"non"} = 0; # 
}
undef @Kmers; # free up memory

my $seqio  = Bio::SeqIO->new(-file => "$infile");
while (my $seqobj = $seqio->next_seq) {
  
  my %words = ();
  
  my $seq = uc($seqobj->seq);
  my $len = length($seq);
  for (my $i=0; $i<$len-$K+1; $i++) {
    
    my $wd = substr($seq, $i, $K);
    if (exists($freq{$wd})) {
      
      $words{$wd} = 1;
      
      if ($DS>0) {
	my $rc = reverse($wd);
	$rc =~ tr/TGCA/ACGT/;
	$words{$rc} = 1;
      }
    }
  }
  foreach my $wd (keys %freq) {
    if (exists($words{$wd})) {
      $freq{$wd}->{"con"} += 1;
    }
    else {
      $freq{$wd}->{"non"} += 1;
    }
  }
}

# Print Output Data Rows
foreach my $word (sort keys %freq) {    
  
  my $cons = $freq{$word}->{"con"};
  my $nonc = $freq{$word}->{"non"};
  printf STDOUT "%s\t%d\t%d\n", ($word, $cons, ($cons+$nonc));
}


#####  SUBROUTINES  #####

sub getallnmers {
  
  my $n     = pop(@_);
  my @alph  = @{shift(@_)};
  #my @nucs = qw(A C G T);
  my @nmers = @alph;
  
  for (my $i=1; $i<$n; $i++) {
    my @old = @nmers;
    @nmers = ();
    
    foreach my $base (@alph) {
      foreach my $word (@old) {
        my $tmp = "$word" . "$base";
        push @nmers, $tmp;
      }
    }
  }
  return @nmers;
}

sub initKmers {
  
  my ($wdfile, $K) = @_;
  my @Kmers = ();
  if ($wdfile ne "") {
    open(WDS,$wdfile) or die "cannot open $wdfile\n";
    while (<WDS>) {
      chomp;
      if ($_=~/^([ACGT]{$K})/) {push @Kmers, $1;}
    }
    close(WDS);
  }
  else {
    my @nucs = qw(A C G T);
    @Kmers = getallnmers(\@nucs, $K);
  }
  return @Kmers;
}
