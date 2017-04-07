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

my $dir_name="NA";
my $outputDir="NA";

setDefaults()                          ;
setParametersFromCommandLine()         ;
printPrnLog("Loading the parameters succesfully\n");

### A|B to 0|1
opendir(DIR, $dir_name) || die "Can't open directory $dir_name"; 
my @dots = readdir(DIR); 
open O, ">$dir_name/mst.map";
foreach (@dots){ 
        if($_=~/genotype/){
           my @id=split(/\./,$_);
           my $cmd="$Bin/MSTMap.exe $dir_name/$id[0].genotype  $dir_name/$id[0].map > out.log";
           system($cmd);
           open F, "<$dir_name/$id[0].map ";
           print O "lg$id[0]\n";
           while(<F>){
             chomp;
             if($_=~/-/){print O "$_\n";}
           }
        }
} 
closedir DIR; 
close O;


sub setParametersFromCommandLine{
    my $result = $tf->TIGR_GetOptions(
           "i=s"                     => \$dir_name,
           "o=s"                     => \$outputDir
   );
   $tf  -> printUsageInfoAndExit() if (!$result);  
   $tf  -> printUsageInfoAndExit() if ($dir_name eq "NA" || $outputDir eq "NA");
}


sub setDefaults{
    my $HELPTEXT = qq~
Sort high-density markers in the ultra-fast mode.
Author : Jinzhuang Dou
Usage:   $PRG file -i  -o  -l -L [options]


Main parameters:

  -i    <s>     The input genotype file . 
  -o    <s>     The output file storing the generated results.
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
