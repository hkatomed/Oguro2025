#!/usr/bin/perl
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# cutoff.pl
#

my $file  = pop(@ARGV);
my $field = 5;
my $locut = 0;
my $hicut = 0;
my $LO = 0;
my $HI = 0;

my $MODE = 1; # 1 => print data outside of [LO,HI] interval, 0 => opposite

for (my $i=0; $i<=$#ARGV; $i++) {
  my $arg = $ARGV[$i];
  if ($arg eq "-incl") {$MODE=0;}
  if ($arg eq "-lo") {$LO=1; $locut = $ARGV[($i+1)];}
  if ($arg eq "-hi") {$HI=1; $hicut = $ARGV[($i+1)];}
  if ($arg =~ /^-f(\d+)$/) {$field = ($1-1);}
  #if ($arg =~ /^-c(\d+\.?\d*)$/) {$cutoff = $1;}
}

open(IN,$file) or die, "cannot open $file: $!\n\n";
while(<IN>) {
  
  chomp();
  my @data = split /\t/, $_;
  my $val = $data[$field];
  if ($MODE>0) {
      
      if ($val<=$locut && $LO==1) {print $_, "\n";}
      if ($val>=$hicut && $HI==1) {print $_, "\n";}
  }
  else {
      
      if (($val>=$locut || $LO==0) && ($val<=$hicut || $HI==0)) {print $_, "\n";}
  }
}
close(IN);
