#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# fa2longfa.pl
#
# reads in fasta file(s), prints them out with the entire sequence on a single line

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $infile  = "";
my $format  = 'fasta';
my $RCflag  = 0;  ###  if 1, use reverse complement
#my $pat = "";
#my $flanklen = 0;

my $usage = "input error\!\n\nfa2longfa.pl -i infile.fa\n";

for (my $i=0; $i<=$#ARGV; $i++) {
    
    my $arg = $ARGV[$i];
    if ($arg eq "-i") {$infile  = $ARGV[($i+1)];}
    #if ($arg eq "-o") {$outfile = $ARGV[($i+1)];}
    if ($arg eq "-f") {$format  = $ARGV[($i+1)];}
    if ($arg eq "-rc") {$RCflag = 1; print STDERR "using reverse complement\n\n";}
    #if ($arg eq "-flank") {$flanklen = $ARGV[($i+1)];}
  }

die "$usage\n" if ($infile eq "");

my $seqio  = Bio::SeqIO->new(-format => "$format", -file => "<$infile");
#my $seqout = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
while (my $seqobj = $seqio->next_seq) {
  
  my $seq = $seqobj->seq;
  my $name = $seqobj->display_name;
  if ($seq ne "") {
  
  if ($RCflag==1) {
    $seq = $seqobj->revcom->seq;
    $name = $name . "_RC";
  }
}
  print STDOUT ">" . $name . "\n";
  print STDOUT $seq . "\n";
}


