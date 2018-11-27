#!/usr/bin/perl -w

###########################################################################
# Present the result of genome assembly using 2bRAD map.
# The input is the <primaryGroup>  <mapFile> <ctgFile> <genotypeFile>
# The format of *.map is as followings:
# MarkerID        groupID   ctgID  geneticDis 
# consensus1        lg1      ctg1      0.000  
# consensus2        lg1      ctg2      0.000  
# consensus3        lg2      ctg2      17.000
# consensus4        lg2      ctg3      0.000         
# consensus5        lg2      ctg3      1.000
# consensus6        lg2      ctg3      7.000
#    ....          ....      ....      ....
#
# The format of *.link is as followings:
#
# Pairs of markers  geneticDis  linkCnt  (ctgID, ctgID)  groupID 
#  (cons1 conse3)      23.000     14       (ctg3, ctg7)    lg5
#    ....              ....      ....          ....       .....
#
# The format of *.path is as followings:
# groupID  ctgID 
# lg1      (ctg1,ctg4,ctg2)
# ...      ......
#
# The format of *.fa is as followings:
# >lg1      (ctg1,ctg4,ctg2)
# ATGCCCCC{NNNNNNNNN}AAAACCCTTTGGTTAA
###########################################################################






use strict;
use warnings;
use Cwd 'abs_path';
use Time::localtime;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use ParseFasta;
use Foundation;
use MSTMAP    ;



# Declare default values for configure

my $tf = new Foundation                ;
my $PRG = $tf->getProgramInfo('name')  ;
my $REV="1.0"                          ;
my @DEPENDS=("Foundation","ParseFasta");
my $VERSION = "Version 1.00"           ;
my %tagPosTable                        ;
my %global                             ;
my $TEXT    =  ""                      ;
my %disTable =()                       ;
my $MSTpath= "$Bin/MSTMap.exe"         ;




my $genotypeFile   = "2bRADtyping.ctg.type";
my $primaryGroup   = "markerGroup.map"     ;
my $ctgFile        = "contigs"             ;
my $mapFile        = "mst.map"             ;
my $outFile        = "mst.map.out"         ;
my $proNum         = 164                   ;
my $prefix         = "djzTest"             ;  
my %genotype       =()                     ;
my %ctgSeq         =()                     ;
my %connect        =()                     ;
my $modelA;
my $modelB;

setDefaults()                          ;

setParametersFromCommandLine()         ;

my @tempVar = (1..$proNum); 
$TEXT=$TEXT."number_of_individual $proNum\nID @tempVar\n";
$global{MSTMapheader} = $TEXT;

my $result = fileCheck($genotypeFile, $primaryGroup, $ctgFile, $mapFile, $MSTpath);

$genotypeFile =abs_path ($genotypeFile) ;
$primaryGroup =abs_path ($primaryGroup) ;
$ctgFile      =abs_path ($ctgFile)      ;
$mapFile      =abs_path ($mapFile)      ;
$MSTpath      =abs_path ($MSTpath)      ;

printPrnLog("Loading the parameters succesfully\n") if $result;

# Read the ctg file with the fasta format

# Read the primary group information generated using ctg sequence 

my $primaryMap = MSTMAP::new($primaryGroup);

# Read the final group information generated using SARAD software
my $mapMap     = MSTMAP::new($mapFile);
# Initilazte the genotype file 
initiGenotype ($genotypeFile) ;

# Input the pre-assemblies 

initiFasta($ctgFile) ;

# Input the link information generated using all the iterative process
getConnect();

printPrnLog("Loading the datasets including $primaryMap $mapMap 
  $genotypeFile $ctgFile succesfully\n");

# Calculate the link information using the based on the class of mstmap
# defined in the package MSTMAP.pm, however, I found some complex in 
# defined the class, as I am not sure how to clarify the classes of node,
# group, and map. I will optimize that programs one day. 

open (O5, ">$prefix".".gapDis" )   || die $!;

my $mapLgNum = $mapMap->getLgNum();
for(my $i=1;$i<=$mapLgNum;$i++){
    my $lgName       = $mapMap ->getLgName($i);
    my $lgNum        = $mapMap ->getMarkerNum($lgName);
    if ($lgNum ==0) {
        printPrnLog("Warining: $lgName has no markers, pls check this\n") ;
    }
    else{
    
        my @lgMarkerList = $mapMap ->getLgOfMarker($lgName);
        # If only one marker in the group, nothing to deal with but just set the distance being 0 cM.
        my $var = $primaryMap->getLgName ( $lgMarkerList[0] );
        $mapMap -> setCtgOfLg( $lgMarkerList[0], $var );
        
        my $dis  =0;

        $mapMap -> setDisOfLg( $lgMarkerList[0], $dis); 
        
        for (my $j=1; $j<$lgNum; $j++){
           my $var1 = $lgMarkerList[$j-1] ;
           my $var2 = $lgMarkerList[$j]   ;
           $dis     = pairWiseDisCal ($var1, $var2) ;
           
           my $physicalDis1 = $mapMap -> getCtgPos( $var1 );  
           my $physicalDis2 = $mapMap -> getCtgPos( $var2 ); 
           $disTable{$dis}=abs($physicalDis1-$physicalDis2) if $dis>0;
           #print "$var1\t$var2\t$dis\t$physicalDis1\t$physicalDis2 \n";
	         print O5 "$dis\t$disTable{$dis}\t0\n" if $dis >0 && $disTable{$dis}<40000;      
           $dis     =0 if $dis <0 ;
           $var     = $primaryMap->getLgName ( $lgMarkerList[$j] );
           $mapMap -> setDisOfLg( $lgMarkerList[$j], $dis);
           $mapMap -> setCtgOfLg( $lgMarkerList[$j], $var ); 

        }

    }
    
}
printPrnLog("All link information has been calulated succesfully\n");

LSTestimation ();


# Output. 
open (O1, ">$prefix".".map" )  || die $!;
open (O2, ">$prefix".".link" ) || die $!;
open (O3, ">$prefix".".path" ) || die $!;
open (O4, ">$prefix".".fa" )   || die $!;


for(my $i=1;$i<=$mapLgNum;$i++){
    my $lgName       = $mapMap ->getLgName($i);
    my @lgMarkerList = $mapMap ->getLgOfMarker($lgName);
    my $lgNum        = $mapMap ->getMarkerNum($lgName);
    my $last         = "";
    my $flag         = "";
    # Only one marker in the group, so I have do nothing to deal with that but just print 
    # print out with the distance being 0.000 cm.
 
    print O3 "$lgName\t";
    print O4 ">scaffold$i\n";
    if($lgNum >0) {
        
        for(my $j=0;$j< @lgMarkerList; $j++){
                 $flag = "";
                 my $index = $lgMarkerList[$j];
                 my $ctgID = $mapMap -> getCtgOfLg( $index );
                
                 my $physicalDis = $mapMap -> getCtgPos( $index );
                 my $var1  = $lgMarkerList[$j-1];
                 my $var2  = $lgMarkerList[$j];
                 my $dis1  = $mapMap -> getDisOfLg( $var1 );
                 my $dis2  = $mapMap -> getDisOfLg( $var2 );

                 if ($ctgID ne $last){
                      $flag = "BreakPoint";
                      my $tempStr = $ctgSeq{$ctgID} ;

                      print O4 "$tempStr".$global{N};
                      print O3 "->($ctgID)" ;
                      my $tempVar = abs($dis1-$dis2);
                      my $gapDis = LSTprediction ($tempVar);
                      printf O2 ("($var1  $var2)\t(%.2f  %.2f)\t($ctgID $last)\t$lgName\n",$dis1, $dis2, $gapDis ) if $last ne "";
                      print O5 "$tempVar\t$gapDis\t1\n" if (int($gapDis)<40000 && int($tempVar)>0);      
     
                 }
                 $last = $ctgID ;
                 printf O1 ("$index\t$lgName\t$ctgID\t%.2f\t$physicalDis\t$flag\n",$dis2);
        }


    }
    print O3 "\n";
    print O4 "\n";

}

close(O1);
close(O2);
close(O3);
close(O4);
close(O5);


printPrnLog("Output succesfully\n");



sub LSTestimation
{
    my $meanX =0;
    my $meanY =0;
    my $meanXX =0;
    my $meanXY =0;
    my $num    =0;
    foreach my $w (keys %disTable){
           if($w > 0){
              $num ++;
              $meanX = $meanX + $w;
              $meanY = $meanY + $disTable{$w};
              $meanXX = $meanXX + $w*$w;
              $meanXY = $meanXY + $w*$disTable{$w};
           }
    }
    $meanX  = $meanX/$num;
    $meanY  = $meanY/$num;
    $meanXX = $meanXX/$num;
    $meanXY = $meanXY/$num;

    $modelB = ($meanXY - $meanX*$meanY) / ($meanXX -  $meanX*$meanX);
    $modelA = $meanY - $modelB* $meanX;
}

sub LSTprediction{

    my $x = shift ;
    return int($modelA+$modelB*$x);

}

sub getConnect
{
   my $cmd="cat ./step*/marker.connect > demo";
   system($cmd);
   open (IN, "<./demo" ) || die $!;
   while(<IN>) {
         chomp;
         my @line =split (/\s+|\t/, $_);
         $line[5] =$line[5]/$line[4];
         $connect{$line[2]}{$line[3]} = "$line[4]\t$line[5]";
         $connect{$line[3]}{$line[2]} = "$line[4]\t$line[5]";

   }
   close(IN);
}



sub printPrnLog
{
    my $msg = shift ;
    $tf -> logLocal($msg,1);
    print STDERR $msg; 

}


sub setParametersFromCommandLine{

    my $result = $tf->TIGR_GetOptions(
           "g=s"                   => \$genotypeFile,
           "p=s"                   => \$primaryGroup,
           "c=s"                   => \$ctgFile,
           "m=s"                   => \$mapFile,
           "o=s"                   => \$outFile,
   );

   $tf  -> printUsageInfoAndExit() if (!$result);
   
}



sub fileCheck{

    for my $var (@_){
      my $tempVar = "$Bin/$var" ;
      if ( (!defined $var || ! -e "$var" ) && (! -e "$tempVar" )) {
          printPrnLog("\n\nError: $var not defined or not found, pls check it!\n\n");
          $tf  -> printUsageInfoAndExit()
      }
    }

    return 1;
}



sub initiGenotype {

    my $var = shift ;
    open (IN, "<$var" ) || die $!;
    while(<IN>) {
         chomp;
         #sleep(10);
         my $ID = myFind ($_, "-", 0) ;
         $genotype{$ID} =$_ ;
    }
    close(IN);
}


sub initiFasta {

    my $name;

    open (IN, "<$_[0]" ) || die $!;
    while(<IN>){ 
          if($_ =~m/^(>.*?)\n/){ 
                 my @line =split(/\s+|\t/, $_);
                 $name =substr($line[0],1); 
                 $name ="lg$name"   ;
                 $ctgSeq{$name} = ''; 
          } 
          else{ 
                 $ctgSeq{$name}.=$_; 
          } 
    }
    close(IN);
} 




sub getGenotype {

    my $var = shift ;
    return $genotype{$var} if exists $genotype{$var};
    return "N" ;
}



sub pairWiseDisCal{

    my $var1 =shift ;
    my $var2 =shift ;
    my $id="$var1-$var2";
    $var1 = getGenotype($var1);
    $var2 = getGenotype($var2);
    my  $diss=0;
    if($var1 ne "N" && $var2 ne "N"){
    system("mkdir -p log");
    open (O, "> ./log/demo.$id.genotype" ) || die $! ;
    print O "$global{MSTMapheader}\n$var1\n$var2\n";
    close O;
 
    my $cmd= $MSTpath." ./log/demo.$id.genotype   ./log/demo.$id.map";
    system("$cmd >> ./log/mst.log");
   
    open (IN, "./log/demo.$id.map" ) || die $!;
    while(<IN>) {
         chomp;
         $diss = myFind ($_, "-", 1) if $_=~/-/;
    }
    close(IN);
    }
    else{
        #print "$id is not in genotype file, pls check it!\n";
    }
    return $diss ;
}



sub myFind{

    my $var1  =shift ;
    my $var2  =shift ;
    my $index =shift ;
    if ($var1 =~ /$var2/ ){
       my @line =split (/\s+|\t/, $var1);
       return $line[$index];
    }
    return undef;
}


sub setParametersFromCommandLine{


    my $result = $tf->TIGR_GetOptions(
           "c=s"                   => \$ctgFile,
           "p=s"                   => \$primaryGroup,
           "g=s"                   => \$genotypeFile,
           "m=s"                   => \$mapFile,
           "o=s"                   => \$prefix,
           "n=i"                   => \$proNum,

   );

   $tf  -> printUsageInfoAndExit() if (!$result);

   
}


sub setDefaults {
    





    my $HELPTEXT = qq~
Output the scaffold file together with the link information generated using SARAD.
Author : Jinzhuang Dou
Usage:   $PRG file -c ctgFile -p primaryGroup -g genoFile -m mst.map -o prefix  -n proNum [options]


Main parameters:

  -c <s>     The input of contig file for scaffolding using SARAD. 
  -p <g>     The primary group determined according to the contigs files.
  -g <s>     The genotype file with marker information being presence and absence generated using HappyMap
             experiment.
  -m <s>     The final map information generated using SARAD.
  -o <s>     The prefix will be added to the filenames of the outputting results.
  -n <n>     Sample size used in HappyMap experiment.    
  -h|help    print this help and exit;
  -V|version print the version and exit;
  -depend    print the program and database dependency list;
  -debug     set the debug <level> (0, non-debug by default); 
        
OUTPUT:  
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;


    $tf->setHelpInfo($HELPTEXT.$MOREHELP);
    $tf->setUsageInfo($HELPTEXT);
    $tf->setVersionInfo($REV);
    $tf->addDependInfo(@DEPENDS);


    
    $TEXT = qq~
population_type DH
population_name LG
distance_function kosambi
cut_off_p_value 0.1
no_map_dist 15.0
no_map_size 0
missing_threshold 1.00
estimation_before_clustering no
detect_bad_data yes
objective_function COUNT
number_of_loci 2
~;


    $global{N} = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
    printPrnLog("Set the default values of parameters\n");

}




