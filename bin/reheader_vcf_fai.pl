#! /usr/bin/env perl5

##reheader of VCF based on an fai
##to include lengths for GATK
use strict;
use warnings;

##input
my $VCF=$ARGV[0];
my $FAI=$ARGV[1];

##iterate over fai, store in hash
open(FAI, $FAI);
my %FAIL;
while(<FAI>){
  my @sp=split(/\t/);
  $FAIL{"##contig=<ID=$sp[0]>"}="##contig=<ID=$sp[0],length=$sp[1]>";
}
close FAI;

##iterate over vcf, print corrected contigs and all other lines
open(VCF, $VCF);
while(<VCF>){
  if($_=~m/^\#\#contig=/){
    chomp;
    if(exists($FAIL{$_})){
      print $FAIL{$_} . "\n";
    }
  }
  else{
    print $_;
  }
}
close VCF;
