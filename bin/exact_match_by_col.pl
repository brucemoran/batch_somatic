#! /usr/bin/env perl5

##exact match two files based on selected columns
use strict;
use warnings;

##input
my @F1=split(/\,/,$ARGV[0]);
my @F2=split(/\,/,$ARGV[1]);

$F1[2]="\t" unless defined $F1[2];
$F2[2]="\t" unless defined $F2[2];

##iterate over fai, store in hash
open(FIL, $F1[0]);
my %FILE;
while(<FIL>){
  chomp;
  my @sp=split($F1[2]);
  $FILE{$sp[$F1[1]]}=$sp[$F1[1]];
}
close FIL;

##iterate over vcf, print corrected contigs and all other lines
open(FIL, $F2[0]);
while(<FIL>){
  chomp;
  my @sp=split($F2[2]);
  if(exists($FILE{$sp[$F2[1]]})){
      print $_ . "\n";
  }
  else{
    next;
  }
}
close FIL;
