#! perl

use strict;
use warnings;
use Cwd;
use File::Basename;

##Lancet VCF format does not contain FORMAT AF
##inputs:
##-ID = Tumour ID (VCF must be somatic, 2 sample, Tumour ID -> other is Germline)
##-DP = Required depth Tumour and Germline
##-MD = Required min depth of ALT in Tumour
##-VCF = VCF input file (in current wd or full path)
##NB hardcoded for 0 ALT in Germline

my ($no,$id,$dp,$md,$vcf)="";
if(scalar(@ARGV)!=4){
  print "Require 4 inputs:\n-ID=Tumour ID (VCF must be somatic, 2 sample, Tumour ID -> other is Germline)\n-DP=Required depth Tumour and Germline (>=)\n-MD=Required min depth of ALT in Tumour (>=)\n-VCF=VCF input file (in current dir or full path)\n";
  exit;
}

for(my $i=0;$i<@ARGV;$i++){
  if($ARGV[$i]=~m/ID/){($no,$id)=split(/\=/,$ARGV[$i]);}
  if($ARGV[$i]=~m/DP/){($no,$dp)=split(/\=/,$ARGV[$i]);}
  if($ARGV[$i]=~m/MD/){($no,$md)=split(/\=/,$ARGV[$i]);}
  if($ARGV[$i]=~m/VCF/){($no,$vcf)=split(/\=/,$ARGV[$i]);}
}

if($vcf!~m/\//){
  my $dir=getcwd();
  my $vcf1=$dir . "/" . $vcf;
  $vcf=$vcf1;
}

##open VCF, run filter
my $tcol=-1;
my $gcol=-2;
my $filtflag=0;
my $header="";
my $snv="";
my $indel="";
my $raw="";
my $outDir="";
my $outName="";

############
## Lancet ##
############

if($vcf=~m/lancet/){
  open(VCF,$vcf);
  while(<VCF>){
    chomp;
    my @sp=split(/\t/);

    if($_=~m/^##/){
      if($_=~m/FORMAT/){
        if($filtflag==0){
          $filtflag++;
          $header.="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=AF,Number=1,Type=String,Description=\"Allele Frequency\">\n##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"allele depth: # of supporting ref,alt reads at the site\">\n##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"strand counts for ref: # of supporting forward,reverse reads for reference allele\">\n##FORMAT=<ID=SA,Number=.,Type=Integer,Description=\"strand counts for alt: # of supporting forward,reverse reads for alterantive allele\">\n##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Total Depth\">\n";
          next;
        }
        else{
          next;
        }
      }
      $header.=$_ . "\n";
      next;
    }

    if($_=~m/^#CHROM/){
      if($sp[-1] eq $id){
        $header.=$_ . "\n";
        next;
      }
      else{
        $tcol=-2;$gcol=-1;
        $header.=$_ . "\n";
        next;
      }
    }

    ##down to pay-dirt
    else{
      my $lastthree=&lancetAF(@sp);
      my $out="";
      if($lastthree eq 0){
        next;
      }
      else{
        my @popd=pop(@sp);
        @popd=pop(@sp);
        @popd=pop(@sp);
        $out.=join("\t",@sp[0..$#sp]);
        $out.="\t$lastthree\n";
        @sp=split(/\t/,$out);
        $raw.=$out;

        if($sp[6] eq "PASS"){
          my $tflag=&lancetTumourFilter(@sp);
          my $gflag=&lancetGermlineFilter(@sp);
          my $iflag=&lancetIndelgermlineFilter(@sp);
          my $siflag=&indelOrSNV(@sp);
          if($siflag==0){
            if(($tflag == 1) && ($gflag == 1)){
              $snv.=$out;
              next;
            }
          }
          else{
            if(($tflag == 1) && ($iflag == 1)){
              $indel.=$out;
              next;
            }
          }
        }
      }
    }
  }
  $outDir=dirname($vcf);
  $outName=$outDir . "/" . $id . ".lancet";
}

#############
## MuTect2 ##
#############
if($vcf=~m/mutect2/){
  open(VCF,$vcf);
  while(<VCF>){
    chomp;
    my @sp=split(/\t/);

    ##header
    if($_=~m/^##/){
      $header.=$_ . "\n";
      next;
    }

    if($_=~m/^#CHROM/){
      if($sp[-1] eq $id){
        $header.=$_ . "\n";
        next;
      }
      else{
        $tcol=-2;$gcol=-1;
        $header.=$_ . "\n";
        next;
      }
    }

    ##variants
    else{
      $raw.=$_ . "\n";
      if(($sp[6] eq "PASS") || ($sp[6] eq ".")){
        my $tflag=&mutect2TumourFilter(@sp);
        my $gflag=&mutect2GermlineFilter(@sp);
        my $iflag=&mutect2IndelgermlineFilter(@sp);
        my $siflag=&indelOrSNV(@sp);
        if($siflag==0){
          if($tflag == 1){
            $snv.=$_ . "\n";
            next;
          }
        }
        else{
          if($tflag == 1){
            $indel.=$_ . "\n";
            next;
          }
        }
      }
    }
  }
  $outDir=dirname($vcf);
  $outName=$outDir . "/" . $id . ".mutect2";
}

####################
## Strelka2 joint ##
####################
if($vcf=~m/strelka2/){
  open(VCF,$vcf);
  while(<VCF>){
    chomp;
    my @sp=split(/\t/);

    if($_=~m/^##/){
      if($_=~m/FORMAT/){
        if($filtflag==0){
          $filtflag++;
          $header.="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GENOTYPE\">\n";
          $header.="##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ALLELE DEPTHS\">\n";
          $header.="##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"ALLELE FREQUENCY (ADalt/DP)\">\n";
          $header.=$_ . "\n";
          next;
        }
      }
      $header.=$_ . "\n";
      next;
    }

    if($_=~m/^#CHROM/){
      if($sp[-1] eq $id){
        $header.=$_ . "\n";
        next;
      }
      else{
        $tcol=-2;$gcol=-1;
        $header.=$_ . "\n";
        next;
      }
    }

    ##down to pay-dirt
    else{
      my ($line,$lasttwo)="";
      my $siflag=&indelOrSNV(@sp);
      if($siflag == 0){ $lasttwo=&Strelka2GTADAF(@sp); }
      if($siflag == 1){ $lasttwo=&Strelka2IndelGTADAF(@sp); }
      if($lasttwo eq 0){ next; }
      else{
        # print "@sp\n";
        # print $siflag . "\n";
        my @popd=pop(@sp);
        @popd=pop(@sp);
        my $pref="GT:AD:AF:";
        my $sp1=$pref . $sp[-1];
        $sp[-1]=$sp1;
        $line.=join("\t",@sp[0..$#sp]);
        $line.="\t$lasttwo\n";
        @sp=split(/\t/,$line);
        $raw.=$line;
        if($sp[6] eq "PASS"){
          # print $line . "\n";
          if($siflag == 1){
            my $itflag=&Strelka2IndelTumourFilter(@sp);
            my $igflag=&Strelka2IndelGermlineFilter(@sp);
            # print "$itflag\t$igflag\t$siflag\n";
            if(($siflag == 1) && ($itflag == 1) && ($igflag == 0)){
                $indel.=$line;
                next;
            }
            else{
              next;
            }
          }
          else{
            my $tflag=&Strelka2SNVTumourFilter(@sp);
            my $gflag=&Strelka2SNVGermlineFilter(@sp);
            # print "$tflag\t$gflag\t$siflag\n";
            if(($tflag == 1) && ($gflag == 1)){
              $snv.=$line;
              next;
            }
            else{
              next;
            }
          }
        }
      }
    }
  }
  $outDir=dirname($vcf);
  $outName=$outDir . "/" . $id . ".strelka2";
}

my $snvindelName=$outName . ".snv_indel.pass.vcf";
open(OUT,">$snvindelName");
print OUT $header;
print OUT $snv;
print OUT $indel;
close OUT;

my $rawName=$outName . ".raw.vcf";
open(OUT,">$rawName");
print OUT $header;
print OUT $raw;
close OUT;

###############################
###############################
####                       ####
#### S U B R O U T I N E S ####
####                       ####
###############################
###############################

############
## Shared ##
############

sub indelOrSNV {
  my @ref=split(//,$_[3]);
  my @alt=split(//,$_[4]);
  if((scalar(@ref) > 1) || (scalar(@alt) > 1)){
      return(1);
  }
  else{
      return(0);
  }
}

############
## Lancet ##
############

sub lancetAF {
  my $gaf="";
  my $taf="";
  my $out="";
  my @tf=split(/\:/,$_[$tcol]);
  my @gf=split(/\:/,$_[$gcol]);
  if(($tf[0] eq "0/0") || ($tf[0] eq ".") || ($gf[0] eq ".")){
    return(0);
  }
  else{

    ##REF, ALT for T, G
    my @tra=split(/\,/,$tf[1]);
    my @gra=split(/\,/,$gf[1]);
    my $tt=$tra[0]+$tra[1];
    my $taff=$tra[1]/$tt;
    my $taf=sprintf("%.3f",$taff);
    my $gt=$gra[0]+$gra[1];
    my $gaff=$gra[1]/$gt;
    my $gaf=sprintf("%.3f",$gaff);

    $out="GT:AD:AF:SR:SA:DP\t";
    if($tcol>$gcol){
      $out.="$gf[0]:$gf[1]:$gaf:" . join(":",@gf[2..$#gf]);
      $out.="\t$tf[0]:$tf[1]:$taf:" . join(":",@tf[2..$#tf]);
    }
    else{
      $out.="$tf[0]:$tf[1]:$taf:" . join(":",@tf[2..$#tf]);
      $out.="\t$gf[0]:$gf[1]:$gaf:" . join(":",@gf[2..$#gf]);
    }
  }
  return($out);
}

sub lancetTumourFilter {
  my @info=split(/\:/,$_[$tcol]);
  if($info[0] eq "."){return(0);}
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] > $md)){
    return(1);
  }
  else{
    return(0);
  }
}

sub lancetGermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  if($info[0] eq "."){return(0);}
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] == 0)){
    return(1);
  }
  else{
    return(0);
  }
}

sub lancetIndelgermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] < $md)){
    return(1);
  }
  else{
    return(0);
  }
}

#############
## MuTect2 ##
#############

sub mutect2TumourFilter {
  my @info=split(/\:/,$_[$tcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] > $md)){
    return(1);
  }
  else{
    return(0);
  }
}

sub mutect2GermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] == 0)){
    return(1);
  }
  else{
    return(0);
  }
}

sub mutect2IndelgermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] < $md)){
    return(1);
  }
  else{
    return(0);
  }
}

##############
## Strelka2 ##
##############
sub Strelka2SNVTumourFilter {
  my @info=split(/\:/,$_[$tcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($refalt[0] eq ".") || ($refalt[1] eq ".")){
    return(0);
  }
  if(($tot > $dp) && ($refalt[1] > $md)){
    return(1);
  }
  else{
    return(0);
  }
}

##total reads above depth threshold && no reads at alt
sub Strelka2SNVGermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] == 0)){
    return(1);
  }
  else{
    return(0);
  }
}

sub Strelka2SNVFindBaseInfo {
  ##which index for base in info?
  my $index="";
  if($_[0] eq "A"){$index=-4}
  if($_[0] eq "C"){$index=-3}
  if($_[0] eq "G"){$index=-2}
  if($_[0] eq "T"){$index=-1}
  return($index);
}

sub Strelka2GTADAF {
  my $ggt="";
  my $tgt="";
  my @tf=split(/\:/,$_[$tcol]);
  if($tf[0]==0){
    return(0);
  }
  else{
    my @gf=split(/\:/,$_[$gcol]);
    my $ref=&Strelka2SNVFindBaseInfo($_[3]);
    my $alt=&Strelka2SNVFindBaseInfo($_[4]);

    ##REF, ALT for T, G
    my @tr=split(/\,/,$tf[$ref]);
    my @ta=split(/\,/,$tf[$alt]);
    my @gr=split(/\,/,$gf[$ref]);
    my @ga=split(/\,/,$gf[$alt]);

    ##make GTs
    ##"A,T,11,2" -> 0/1:11,2:0.154:13
    ##"C,G,11,0" -> 0/0:11,0:0:11

    ##germline 0/0
    if($ga[0]==0){
      $ggt="0/0:$gr[0],0:0:";
    }
    ##tumour 0/1
    if(($tr[0]>0) && ($ta[0]>0)){
      my $tot=$tr[0]+$ta[0];
      my $fq=$ta[0]/$tot;
      my $fqr=sprintf("%.3f", $fq);
      $tgt="0/1:$tr[0],$ta[0]:$fqr:";
    }
    ##germline 0/1
    if(($gr[0]>0) && ($ga[0]>0)){
      my $tot=$gr[0]+$ga[0];
      my $fq=$ga[0]/$tot;
      my $fqr=sprintf("%.3f", $fq);
      $ggt="0/1:$gr[0],$ga[0]:$fqr:";
    }
  }
  if(($tgt eq "") || ($ggt eq "")){
    return(0);
  }
  elsif($tcol eq -2){
    my $out="$tgt$_[$tcol]\t$ggt$_[$gcol]";
    return($out);
  }
  else{
    my $out="$ggt$_[$gcol]\t$tgt$_[$tcol]";
    return($out);
  }
}

sub Strelka2IndelTumourFilter {
  my @info=split(/\:/,$_[$tcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($refalt[0] eq ".") || ($refalt[1] eq ".")){
    return(0);
  }
  if(($tot > $dp) && ($refalt[1] > $md)){
    return(1);
  }
  else{
    return(0);
  }
}

sub Strelka2IndelGermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] < $md)){
    return(1);
  }
  else{
    return(0);
  }
}

sub Strelka2IndelGTADAF {

  ##make GT:AD:AF
  my $ggt="";
  my $tgt="";

  my @ti=split(/\:/,$_[$tcol]);
  my @gi=split(/\:/,$_[$gcol]);
  my @tref=split(/\,/,$ti[2]);
  my @talt=split(/\,/,$ti[3]);
  my @gref=split(/\,/,$gi[2]);
  my @galt=split(/\,/,$gi[3]);

  ##make GTs
  ##"A,T,11,2" -> 0/1:11,2:0.154:13
  ##"C,G,11,0" -> 0/0:11,0:0:11

  ##germline 0/0
  if($galt[0]==0){
    my $ggt="0/0:$gref[0],0:0:";
  }
  ##germline 0/1
  if(($gref[0]>0) && ($galt[0]>0)){
    my $tot=$gref[0]+$galt[0];
    my $fq=$galt[0]/$tot;
    my $fqr=sprintf("%.3f", $fq);
    $ggt="0/1:$gref[0],$galt[0]:$fqr:";
  }
  ##germline 1/1
  if(($gref[0]==0) && ($galt[0]!=0)){
    $ggt="1/1:0,$galt[0]:1:";
  }
  ##tumour 0/0
  if($talt[0]==0){
    my $tgt="0/0:$tref[0],0:0:";
  }
  ##tumour 0/1
  if(($tref[0]>0) && ($talt[0]>0)){
    my $tot=$tref[0]+$talt[0];
    my $fq=$talt[0]/$tot;
    my $fqr=sprintf("%.3f", $fq);
    $tgt="0/1:$tref[0],$talt[0]:$fqr:";
  }
  ##tumour 1/1
  if(($tref[0]==0) && ($talt[0]!=0)){
    $tgt="1/1:0,$talt[0]:0:";
  }

  if(($tgt eq "") || ($ggt eq "")){
    return(0);
  }
  if($tcol == -2){
    my $out="$tgt$_[$tcol]\t$ggt$_[$gcol]";
    return($out);
  }
  else{
    my $out="$ggt$_[$gcol]\t$tgt$_[$tcol]";
    return($out);
  }
}
