#!/usr/bin/perl

use strict;use warnings;
use lib "/home/yhmao/script/all_codes_backup/Script/perl_module/";
use START_POS;

my $start_pos=shift;
my $cds=shift;
my $align=shift;

my $flag=shift;

my %start_pos;my %cds;
START_POS($start_pos,\%start_pos);

open(F,$cds);
while(my $line=<F>){
	my $line1=<F>;
	$line=~s/\s+$//;$line=~s/>//;
	$line1=~s/\s+$//;
	$cds{$line}=$line1;
}
close(F);

open(F,"samtools view $align|");
my %hash;my %hash1;my %all;my $all1=0;
while(my $line=<F>){
	$line=~s/\s+$//;
	my @tps=split(/\t/,$line);

	my $name=$tps[2];
	my $pos=$tps[3]-1+12;
	next unless(exists $start_pos{$name});
	next if($pos<$start_pos{$name}{'start'}-1 or $pos>$start_pos{$name}{'start'}+1);
	
	my $CIGAR=$tps[5];
	if($CIGAR=~/^\d{0,}S{0,}\d+M\d{0,}S{0,}$/){
		my($match_len)=$CIGAR=~/(\d+)M/;
		#		print $CIGAR,"\t",$match_len;<STDIN>;
		if($CIGAR=~/\d+S$/){
			$hash{$name}{$match_len}++;$hash1{$name}++;

			$all{$match_len}++;$all1++;
		}
		elsif($CIGAR=~/\d+M$/){
			my $pos1=$tps[3]-1+$match_len;
			my $base=substr($cds{$name},$pos1,15);
			my ($string_a)=$base=~/^(A{0,})/;
			my $length_a=length($string_a);
			
			for(my $i=$match_len;$i<=$match_len+$length_a;$i++){
				my $frac=1/($length_a+1);
				$hash{$name}{$i}+=$frac;
				$all{$i}+=$frac;
			}
			
			$hash1{$name}++;$all1++;
		}
	}
}
close(F);

my %mean;my %mean_count;
foreach my $name(keys %hash){
	foreach my $len(keys %{$hash{$name}}){
		$mean{$name}+=($len*$hash{$name}{$len});
		$mean_count{$name}+=$hash{$name}{$len};
	}
}

my $mean_all=0;my $mean_all_count=0;
foreach my $len(keys %all){
	$mean_all+=($len*$all{$len});
	$mean_all_count+=$all{$len};
}
$mean_all/=$mean_all_count;

open(F,"common_name_top50.txt");
open(F1,">start_length_top50_$flag");
while(my $line=<F>){
	$line=~s/\s+$//;
	my($name1,$name2)=split(/\t/,$line);
	if(exists $hash{$name1}){
		$mean{$name1}/=$mean_count{$name1};
		for(my $i=20;$i<35;$i++){
			unless(exists $hash{$name1}{$i}){
				print F1 "$name1\_$name2\t$i\t0\t$flag\t$mean{$name1}\n";
			}
			else{
				my $frac=$hash{$name1}{$i}/$hash1{$name1};
				print F1 "$name1\_$name2\t$i\t$frac\t$flag\t$mean{$name1}\n";
			}
		}
	}
}
close(F);


for(my $i=20;$i<35;$i++){
	unless(exists $all{$i}){
		print F1 "all\t$i\t0\t$flag\t$mean_all\n";
	}
	else{
		my $frac=$all{$i}/$all1;
		print F1 "all\t$i\t$frac\t$flag\t$mean_all\n";
	}
}
close(F1);
