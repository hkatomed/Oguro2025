#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# selectfaseqs.pl
#

use strict;
use Bio::SeqIO;

my $infile = "";
my $outfile = "";
my $listfile = "";
my $MINLEN = 0;
my $MAXLEN = 0;
my $PRINTLEN = 0;

my $usage = "input error\!\n\nselectfaseqs.pl -list listfile -i infile.fa -o outfile.fa [-minlen X] [-maxlen X] [-printlen]\n";
my %list = ();

for (my $i=0; $i<=$#ARGV; $i++) {
  
  my $arg = $ARGV[$i];
  if ($arg eq "-maxlen") {$MAXLEN = $ARGV[($i+1)];}
  if ($arg eq "-minlen") {$MINLEN = $ARGV[($i+1)];}
  if ($arg eq "-printlen") {$PRINTLEN = 1;}
  if ($arg eq "-i") {$infile = $ARGV[($i+1)];}
  if ($arg eq "-o") {$outfile = $ARGV[($i+1)];}
  if ($arg eq "-list") {$listfile = $ARGV[($i+1)];}
}

die "$usage\n" if ($infile eq "" || $listfile eq "" || $outfile eq "" || $MAXLEN<0);

my $rank = 1;
open(IN,$listfile) or die "cannot open $listfile\n";
while(<IN>) {
  chomp;
  $list{$_} = $rank++;
}
close(IN);

#my $seqio  = Bio::SeqIO->new(-format => 'fasta', -file => "<$infile");
my $seqio  = Bio::SeqIO->new(-file => "<$infile");
my $seqout = Bio::SeqIO->new(-format => 'fasta', -file => ">$outfile");
while (my $seqobj = $seqio->next_seq) {
  if (exists($list{$seqobj->display_name})) {
    if (length($seqobj->seq)>=$MINLEN && (length($seqobj->seq)<=$MAXLEN || $MAXLEN<=0)) {
      if ($PRINTLEN==1) {
	my $tmp = $seqobj->display_name . "_" . $seqobj->length . "bp";
	$seqobj->display_name($tmp);
      }
      $seqout->write_seq($seqobj);
    }
    else {
      print STDERR $seqobj->display_name . " has length outside range $MINLEN to $MAXLEN\n";
    }
  }
}



