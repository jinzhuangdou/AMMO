#!/usr/bin/perl -w

###########################################################################
# Group  markers into different linkage groups using the ultra-fast version 
############################################i###############################
# Input : fasta   (de novo)
# Output: pair-end file
###########################################################################

use strict;
use warnings;
use Cwd 'abs_path';
use Time::localtime;
use DBM::Deep;
use Graph::Undirected;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use Foundation;

# Declare default values for configure

my $tf = new Foundation                ;
my $PRG = $tf->getProgramInfo('name')  ;
my $REV="1.0"                          ;
my @DEPENDS=("Foundation","ParseFasta");
my $VERSION = "Version 1.00"           ;

 
my $inputFile="NA";
my $outputDir="NA";
my $samplingNum;
my $groupThreshold;
my $clusterNum=100;



setDefaults()                          ;
setParametersFromCommandLine()         ;
my $threshold=$samplingNum*$groupThreshold;

printPrnLog("Loading the parameters succesfully\n");
system("mkdir -p $outputDir");
### A|B to 0|1
open (IN, "<$inputFile" )      ||  $tf->bail("Could't open $inputFile ($!)");
open (O, ">$outputDir/AB_01" )      ||  $tf->bail("Could't open $inputFile ($!)");
my %strToNum=(A=>"0",B=>"1",0=>"0",1=>"1");  
my %numToStr=(0=>"A",1=>"B",A=>"A",B=>"B");  
my %ID=();
my $cnt=0;
my $sampleSize=0;
my $g=Graph::Undirected->new;

while(<IN>){
  chomp;
  if($_=~/-/){
      my @line=split(/\s+|\t/,$_);
      $cnt++;
      $ID{$cnt}=$_;
      print O "$line[0]";
      for(my $i=1;$i<@line;$i++){
      	  $sampleSize=@line-1;
          print O "\t$strToNum{$line[$i]}";
      }
      print O "\n";
  }
}
close (IN);
close (O);


### Get the cluster information for all markers 

system("Rscript $Bin/markerSort.r  $outputDir/AB_01  $samplingNum  $groupThreshold 0.8 1 1  $clusterNum  $outputDir/cluster");
system("rm -rf $outputDir/*.map");
system("rm -rf $outputDir/*.genotype");
### Group markers into different markers 
open IN, "<$outputDir/cluster";
my %hash=();


my $ tmp_name="";
my %link=();
while(<IN>){
  chomp;
  print "High memory usage now...\n";
  my @line=split(/\s+|\t|,/,$_);
  for(my $i=0;$i<@line;$i++){
  	 $hash{$line[$i]}++;
  	 $g->add_vertex($line[$i]);
     for(my $j=$i+1;$j<@line;$j++){  
        if($line[$j]<$line[$i]){
           $tmp_name="$line[$j]"."_"."$line[$i]";
        }
        else{
           $tmp_name="$line[$i]"."_"."$line[$j]";
        }
        $link{$tmp_name}++;
        if($link{$tmp_name}>$threshold){
        	$g->add_edge($line[$j],$line[$i]);
        }
     }
  }
}
close IN;

## Release the memory of hash table
%link=();


### Genotype file for each group 
my @cluster=$g->connected_components();
my $i=0;
my $groupNum=0;
open OUT, ">$outputDir/secondLevelRun.sh";
foreach my $x (@cluster){
	$i=$i+1;
	$groupNum=$i;
	open O, ">$outputDir/$i.genotype";
	my $markerNum=0;
	 foreach my $y (@$x){ $markerNum++;}
	 print O "population_type DH\n";
	 print O "population_name LG\n";
	 print O "distance_function kosambi\n";
	 print O "cut_off_p_value 1\n";
	 print O "no_map_dist 15.0\n";
	 print O "no_map_size 0\n";
	 print O "missing_threshold 1.00\n";
	 print O "estimation_before_clustering no\n";
	 print O "detect_bad_data yes\n";
	 print O "objective_function COUNT\n";
	 print O "number_of_loci $markerNum\n";
	 print O "number_of_individual $sampleSize\n";
	 foreach my $y (@$x){ $markerNum++;}
	 print O "ID";
	 for(my $j=1;$j<=$sampleSize;$j++){print O " P$j";}
	 print O "\n";
   my $markerCNT=0;
	 foreach my $y (@$x){ 
  	 	my @line=split(/\s+|\t/,$ID{$y});
  	 	print O "$line[0]";
      $markerCNT++;
  	 	for(my $j=1;$j<@line;$j++){
  	 		print O " $numToStr{$ line[$j]}";
  	 	}
  	 	print O "\n";
  }
  close O;
  if($markerCNT>5000){
    print OUT "The linkageMap $i has more than $markerCNT markers, further split will be executed!\n"; 
    system("mkdir $outputDir/$i.out");
    if($markerCNT>20000){
       system("perl  $Bin/GroupFast.pl   -i $outputDir/$i.genotype  -o $outputDir/$i.out  -s 3 -T 0.9  -c 100");
    }
    else{
       system("perl   $Bin/GroupFast.pl   -i $outputDir/$i.genotype  -o $outputDir/$i.out  -s 3 -T 0.9  -c 100");
    }
  }
}
close OUT;



sub setParametersFromCommandLine{
    my $result = $tf->TIGR_GetOptions(
           "i=s"                     => \$inputFile,
           "o=s"                     => \$outputDir,
           "s=s"                     => \$samplingNum,
           "T=f"                     => \$groupThreshold,
           "c=n"                     => \$clusterNum
   );
   $tf  -> printUsageInfoAndExit() if (!$result);  
   $tf  -> printUsageInfoAndExit() if ($inputFile eq "NA" || $outputDir eq "NA");
}


sub setDefaults{
    my $HELPTEXT = qq~
Sort high-density markers in the ultra-fast mode.
Author : Jinzhuang Dou
Usage:   $PRG file -i  -o  -l -L [options]


Main parameters:

  -i    <s>     The input genotype file . 
  -o    <s>     The output file storing the generated results.
  -s    <s>     The maximal sampling times allowed for markers grouping and markers sorting.
  -T    <s>     Threshold over which pairs of contigs will be linked together under the specified sampling times. [0.6] 
  -c    <s>     A scalar denoting the number of clusters. [10]
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
