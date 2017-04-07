#!/usr/bin/perl -w

###########################################################################
# Generate the input file for HappyMap experiment 
############################################i###############################
# Input : fasta   (de novo)
# Output: pair-end file
###########################################################################


use strict;
use warnings;
use Cwd 'abs_path';
use Time::localtime;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use Statistics::Regression;
use ParseFasta;
use Foundation;




# Declare default values for configure

my $tf = new Foundation                ;
my $PRG = $tf->getProgramInfo('name')  ;
my $REV="1.0"                          ;
my @DEPENDS=("Foundation","ParseFasta","Regression");
my $VERSION = "Version 1.00"           ;

 
my $INPUT_MAP;
my $INPUT_GENO;
my $OUTPUT;
my $INPUT_LINK;
my %dis=();


setDefaults()                          ;
setParametersFromCommandLine()         ;


printPrnLog("Loading the parameters succesfully\n");

open (IN1, "<$INPUT_MAP" )      ||  $tf->bail("Could't open $INPUT_MAP ($!)");
open (O,  ">$OUTPUT" )     ||  $tf->bail("Could't open $OUTPUT ($!)");



my $tmp;
my %markerMap=();
my $lg;
my @f;
while(<IN1>){
      chomp;
      @f=split(/\t|\s+/,$_);
      if(@f<2){
                   $lg=$f[0];
                   $markerMap{$lg}{cnt}=0;
                   $markerMap{$lg}{name}=$f[0];

      }
      else{

                  $markerMap{$lg}{cnt}++;
                  my $tag=$markerMap{$lg}{cnt};
                  $markerMap{$lg}{$tag}{name}=$f[0];
                  
                  $markerMap{$lg}{$tag}{phydis}=$f[1];
      }
}
close(IN1);



open (IN2, "<$INPUT_GENO" )  ||  $tf->bail("Could't open $INPUT_GENO ($!)");
my %geno=();
while(<IN2>){
      chomp;
      my @f=split(/\t|\s+/,$_);
      $geno{$f[0]}=$_;
}
close(IN2);


foreach my $w(keys %markerMap){
      if($markerMap{$w}{cnt}==1){print O "$w\t0\t0\n";}
      else{
        for(my $i=1;$i<$markerMap{$w}{cnt};$i++){
			for(my $j=$i+1;$j<=$markerMap{$w}{cnt};$j++){
				 my $m1=$markerMap{$w}{$i}{name};
				 my $m2=$markerMap{$w}{$j}{name};
				 if(exists $geno{$m1} & exists $geno{$m2}){
						 my $gendis1=$markerMap{$w}{$i}{phydis};
						 my $gendis2=$markerMap{$w}{$j}{phydis};
						 my $tmp=abs($gendis2-$gendis1);
						 #print "$gendis1  $gendis2 $tmp  $tmp1\n";
						 my $tmp1=pairWise($geno{$m1},$geno{$m2});
						 print O "$tmp1\t$tmp\t0\n";
				 }
            }
        }

      }
}
open (IN2, "<$INPUT_LINK" )  ||  $tf->bail("Could't open $INPUT_LINK ($!)");
#open (O,  ">test.mat" )     ||  $tf->bail("Could't open $OUTPUT ($!)");
my $ctg1;
my $ctg2;
my $line=0;
while(<IN2>){
      chomp;
      @f=split(/\t|\s+/,$_);
      my $m1=$f[0];
      my @tmp=split(/-/,$m1);
      my $dis1=$tmp[1];
      my $m2=$f[1];
      @tmp=split(/-/,$m2);
      my $dis2=$tmp[1];
      my $phydis=abs($dis2-$dis1);
      my $genedis=pairWise($geno{$m1},$geno{$m2});
	  print O "$genedis\t$phydis\t1\n";


}
close(IN1);
close(O);



sub pairWise{
    my $tmp=shift;
    my @m1=split(/\t|\s+/,$tmp);
    $tmp=shift;
    my @m2=split(/\t|\s+/,$tmp);
    my $total=0;
    my $diff=0;
    my $myval=0;
    for(my $i=1;$i<@m1;$i++){
      $total++;
      if($m1[$i] ne $m2[$i]){$diff++;}
    }
    $diff=$diff/$total;
    if($diff>0.00001){
	#$myval=-$diff*log($diff)-(1-$diff)*log(1-$diff);
    }
    return $diff;   
}

  





sub setParametersFromCommandLine{

    my @specOpts=();

    my $result = $tf->TIGR_GetOptions(
           "m=s"                      => \$INPUT_MAP,
           "l=s"                      => \$INPUT_LINK,
           "g=s"                      => \$INPUT_GENO,
           "o=s"                      => \$OUTPUT,
    );

    if (@specOpts) {
        print "$0: Argument required.\n";
        usage();
    }
    usage() if  (!$result);
    $tf  -> printUsageInfoAndExit() if (!$result);
}


sub setDefaults{
    my $HELPTEXT = qq~
# Generate the input file for HappyMap experiment 
Author : Jinzhuang Dou
Usage:   $PRG file -i  -o 


Main parameters:

  -i    <s>     The input file with *.map format . 
  -o    <s>     The output file imported for the Matlab.
  -h|help       print this help and exit;
  -V|version    print the version and exit;
  -depend       print the program and database dependency list;
  -debug        set the debug <level> (0, non-debug by default); 
        
OUTPUT:  
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;


    $tf->setHelpInfo($HELPTEXT.$MOREHELP);
    $tf->setUsageInfo($HELPTEXT);
    $tf->setVersionInfo($REV);
    $tf->addDependInfo(@DEPENDS);


}



sub printPrnLog
{
    my $msg = shift ;
    $tf -> logLocal($msg,1);
    print STDERR $msg; 

}
