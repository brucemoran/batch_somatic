#! perl5

use warnings;
use strict;
use Sort::Naturally;

##script to take in multiple single-sample VCFs
##first files header is retained as header.txt, rest are lost
##take all extant variants (no # prefix; $ID_$f[0]_$f[1] as hash key -> $_;)
##if found, report GT string; if not report '-'
##report values for all extant postions in a table
##report full VEP string annotation for each row

##set inputs; runs on a dir (currdir) of VCFs
my @vcfs = `ls *vep.vcf`;
my $hcount = 0;

##capture variants (rows not prefixed #)
my %vars;
my %chrsort1;
my %chrsort2;
my @IDs;
my @ID_vars;
my $header;
my $chrom;

##iterate over VCFs in dir
foreach my $f (@vcfs){

  print "Reading in: " . $f;
  chomp $f;

  open(FIL, $f);

  while(<FIL>){
    if($hcount == 0){
      if($_=~m/^##/){
        $header .= $_;
      }
    }
    if($_=~m/^##/){
      next;
    }
    chomp $_;
    my @sp = split(/\t/);

    ##capture ID into array
    if($_=~m/^#CHROM/){
      push(@IDs, $sp[-1]);
      $chrom = join("\t", $sp[2], $sp[3], $sp[4], $sp[5], $sp[6]);
    }

    ##parse variants
    else{
      $vars{$sp[0] . "_" . $sp[1]} = join("\t", $sp[2], $sp[3], $sp[4], $sp[5], $sp[6], $sp[7]);
      $chrsort1{$sp[0]} = $sp[1];
      $chrsort2{$sp[0] . "_" . $sp[1]} = $sp[0] . "_" . $sp[1];
      $ID_vars[$hcount]{$IDs[$hcount] . "_" . $sp[0] . "_" . $sp[1]} =  $sp[-1];
    }
  }

  $hcount++;
  close FIL;
}

##print header
open(OUTH, ">header.txt");
print OUTH "$header\n";
close OUTH;

##now we have all vars, and all vars per sample
##we parse this and print output
##first header line
open(OUT, "> " . $ARGV[0] . ".tabvcf.tsv");
print OUT "variant\t" . join("\t", @IDs) . "\t$chrom\tVEP\n";

foreach my $ky (nsort keys %vars){

  ##therefore key of %vars, values of which are the field info
  ##setup output line
  my $oline = $ky . "\t";

  ##iterate over IDs, each value into output line
  for(my $i=0; $i<@IDs; $i++){
    ##therefore use $IDs[$i] . "_" . $k as key of $id_vars to get FORMAT data to populate table

    if($ID_vars[$i]{$IDs[$i] . "_" . $ky}){
      $oline .=  $ID_vars[$i]{$IDs[$i] . "_" . $ky} . "\t";
    }
    else {
      $oline .= "-\t";
    }
  }
  $oline .= $vars{$ky} . "\n";
  print OUT $oline;
}
close OUT;
