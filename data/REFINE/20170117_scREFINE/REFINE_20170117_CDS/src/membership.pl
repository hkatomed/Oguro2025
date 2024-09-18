#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# membership.pl
#
# this program reads in a listfile and a text file and prints out lines
# from the text file that have an entry in the first field that is present
# in the list file (tab-delimited).
#

my $usage = "membership.pl -list listfile -i infile\n";
my %list = ();
my $listfile = "";
my $infile   = "";

for (my $i=0; $i<$#ARGV; $i++) {
  my $arg = $ARGV[$i];
  if ($arg eq "-list") {$listfile = $ARGV[($i+1)];}
  if ($arg eq "-i")    {$infile = $ARGV[($i+1)];}
}
if ($infile eq "" || $listfile eq "") {die "$usage\n";}

open(LIST,$listfile) or die "cannot open list $listfile\n";
while(<LIST>) {
  chomp;
  s/\s+//g;
  $list{$_} = 1;
}
close(LIST);

#foreach my $k (sort keys %list) {print ":$k:\n";}

open(IN,$infile) or die "cannot open input file $infile\n";
while(<IN>) {
  
  chomp($_);
  s/\s+/\t/g;
  my @data = split /\t/, $_;
  my $name = $data[0];
  if (exists($list{$name})) {
    print "$_" . "\t" . "1\n";
  }
  else {
    print "$_\t0\n";
  }
}
close(IN);
