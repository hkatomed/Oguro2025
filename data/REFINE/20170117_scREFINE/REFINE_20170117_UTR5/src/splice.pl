#!/usr/bin/perl -w
#
# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# splice.pl
#
# usage:  splice.pl f1 ... fN
#


my @lines;
my $Nlines=0;
my $Nfiles=0;

#make sure at least two input files given
if ($#ARGV < 1) {die "input error!\n";}
else {
  
  foreach $file (@ARGV) {

    $Nfiles++;
    my $n=0;
    open (FILE, "$file") or die "cannot open $file\n $!\n";
    
    while (<FILE>) {
      chomp($_);
      
      if ($Nfiles<=1) {
	push @lines, $_;
      }
      else {
	
	#add nth line from fileN to @lines
	if ($n<=$Nlines) {
	  my $tmp = $lines[$n];
	  $lines[$n] = $tmp . "\t" . $_;
	}
	#add (N-1) tabs to @lines before adding nth line from fileN
	else {
	  my $tmp = "";
	  for (my $i=1; $i<$Nfiles; $i++) {
	    $tmp .= "\t";
	  }
	  $lines[$n] = $tmp . $_;
	}
      }
      $n++;
    }
    $Nlines = $#lines;
    close FILE;
  }
}

foreach $line (@lines) {print STDOUT $line . "\n";}
