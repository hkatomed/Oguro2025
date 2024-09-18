#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# updatelodfile.pl
#
# this program reads in a ss.* file and a lod file
# a cutoff score is calculated from the ssfile as the
# score (column 1) with the minimum finite log-p-value
# (column 2) from the ssfile.
# the lodfile is read and its previous cutoff value
# is replaced with the new cutoff score from ssfile
#

my ($ssfile, $lodfile) = @ARGV;
my $cutscore = `grep -v Inf $ssfile | sort -n -k 2 | head -1 |cut -f1`;
chomp($cutscore);
#print STDERR "cutoff=\t$cutscore\n";

my $DONE = 0;
open(IN,$lodfile) or die "cannot open $lodfile\n";
while(<IN>) {
  
  if (m/CUTOFF=/) {
    print "CUTOFF= $cutscore\n";
    $DONE = 1;
  }
  else {print $_;}
}
#if ($DONE==0) {print STDERR "WARNING:\tNO CUTOFF present in $lodfile\n\n";}
