#! perl5

use warnings;
use strict;

##script to take in two files, match on first two elements, return cat'd line
##file_1 is master, return cat'd to it's lines
##any double hash lines returned without change
##search file_2 for lines with first 2 elements starting same as file_1
##cat last element of file_2 lines to matching file_1

##set inputs; runs on a dir (currdir) of VCFs
my $file_1 = $ARGV[0];
my $file_2 = $ARGV[1];
my $outf = $ARGV[2];

##capture variants (rows not prefixed #)
my %keyd;
my $outs;

##iterate over VCFs in dir
if($file_1 =~ m/".gz$"/){
  open(FIL1, "gunzip -c $file_1 |");
}
else{
  open(FIL1, $file_1);
}

while(<FIL1>){
  if($_=~m/^##/){
      $outs .= $_;
  }

  else{
    chomp $_;
    my @sp = split(/\t/);
    $keyd{$sp[0] . "\t" . $sp[1]} = $_;
  }
}

close FIL1;

open(FIL2, $file_2);

while(<FIL2>){
  if($_=~m/^##/){
      next;
  }

  else{
    chomp $_;
    my @sp = split(/\t/);
    if($keyd{$sp[0] . "\t" . $sp[1]}){
      $outs .= $keyd{$sp[0] . "\t" . $sp[1]} . "\t" . $sp[-1] . "\n";
    }
  }
}

close FIL2;

##print output
open(OUTF, ">$outf");
print OUTF $outs;
close OUTF;
