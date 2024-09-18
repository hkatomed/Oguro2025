#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# meme2blocks.pl
#
# this program reads in a meme output file and prints out the motif sites
# (given in blocks format table) out in fasta format
#
#

my $FLAG=0;      # 1 => reading blocks format table
my $N = 0;       # Motif number
my $MOTIFN = 1;  # print only this motif
my $infile = pop(@ARGV);

for (my $i=0; $i<$#ARGV; $i++) {
  my $arg = $ARGV[$i];
  if ($arg eq "-n") {$MOTIFN = $ARGV[($i+1)];}
}

open(IN,"$infile") or die "cannot open $infile\n";
while (<IN>) {
  
  if (m/^BL/) {$FLAG=1; $N++;}
  if ($FLAG==1 && m/^\/\//) {$FLAG=0;}
  if ($FLAG==1 && m/\(\s*\d+\)\s*[ACGTNX]+/) {
    s/\s+\(\s*/\_/;
    s/\s*\)\s*/\n/;
    s/\n([ACGTXN]+).+/\n$1/;
    if ($N==$MOTIFN) {
      s/^.+\n//;
      print $_;
    }	
  }
}
close(IN);
