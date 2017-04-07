#!/usr/bin/perl -w
use Cwd 'abs_path';
use Time::localtime;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use ParseFasta;
use Foundation;


# Declare default values for configure


my $tf = new Foundation;
my $PRG = $tf->getProgramInfo('name');
my $REV="1.0";
my @DEPENDS=("Foundation","ParseFasta");
my $VERSION = "Version 1.00";

# help info
my $HELPTEXT = qq~
Program that generate the simulated HAPPYMAP datasets 
Program: 2bRADsim.pl
Author : Jinzhuang Dou
Usage: $PRG file -r refence_path -c enzymes -l len [options]
        
        -r <s>       reference sequence used to generated the 2b-RAD tags.
        -c <s>       type of IIB restriction enzymes, input format is as followings. 
                     [AGTCN]{30}AC[AGTCN]{5}CTCC[ATGCN]{30}-[ATGCN]{30}GGAG[AGTCN]{5}GT[ATGCN]{30};
        -l <s>       length of restriction enzymes, input format is as followings: [71-71];
        -cl<n>       clone size used in the simulated HappyMap experiment. [40000]";
                     The size of BAC clone is 100kb, and 40kb for Fosmid clone.;
        -s <n>       number of samples in the simulated HappyMap experiment. [100];
        -d <f>       average coverage of clone along the genome sequence for each sample. [0.5];
        -o <s>       specified output dir;    

        -h|help      print this help and exit;
        -V|version   print the version and exit;
        -depend      print the program and database dependency list;
        -debug       set the debug <level> (0, non-debug by default); 
 
  OUTPUT:  
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;


    
# Declare default values for variables
my $CHR_LEN       = 3949921;
my $ref           ="./testRef";
my $cloneSize     = 40000 ;
my $cloneCoverage = 0.5;
my $sampleSize    = 100;
my $enzyme        = "[AGTCN]{30}AC[AGTCN]{5}CTCC[ATGCN]{30}-[ATGCN]{30}GGAG[AGTCN]{5}GT[ATGCN]{30}";
my $outdir        = "simtest";


# Set parameter according to the commandline specified by the users
# Configure TIGR Foundation
        
$tf->setHelpInfo($HELPTEXT.$MOREHELP);
$tf->setUsageInfo($HELPTEXT);
$tf->setVersionInfo($REV);
$tf->addDependInfo(@DEPENDS);
        
my $result = $tf->TIGR_GetOptions(
           "r=s"                   => \$ref,
           "c=s"                   => \$enzyme,
           "cl=i"                  => \$cloneSize,
           "d=f"                   => \$cloneCoverage,
           "s=i"                   => \$sampleSize,
           "o=s"                   => \$outdir,
);



$ref    = abs_path($ref);
$outdir = abs_path($outdir);
$tf     ->printUsageInfoAndExit() if (!$result);
$tf     ->logLocal("Reading fasta records from $ref",9);


open (IN, "<$ref" ) ||  $tf->bail("Could't open $ref ($!)");
my $fr        = new ParseFasta(\*IN);
my %tagList   = getTags();


$tf     ->logLocal("Return the hash list from the reference sequence",9);

print "$CHR_LEN\n";
my $cloneNum = $cloneCoverage * $CHR_LEN / $cloneSize;


system("mkdir $outdir") if (! -e "$outdir");

for ($i =1;$i <= $sampleSize;$i ++) 
{

     #Return the random postion of tags 
     my @list = generateSample() ;

     open (O, ">$outdir/$i".".fa") ||  $tf->bail("Could't write $outdir($!)");
    
     foreach my $w(@list) 
     {
            print O ">$i"."-".$tagList{$w}->{header}."\n";
            print O $tagList{$w}->{body}."\n";      
     }
     close O;

     $tf->logLocal("Generate the $i th samples.",9);

}

exit 0;



sub generateSample
{

    my @list =();

    for(my $i = 1;$i <= $cloneNum; $i++) {

        my $tagStart =int (rand(1) * $CHR_LEN );
        my $tagEnd   =$tagStart + $cloneSize  ;

        foreach my $w(keys %tagList ) 
        {
              push @list, $w  if ($w >= $tagStart &&  $w <= $tagEnd)
        }
    }
    
    return @list;

}


#Get the tags that has the enzyme cut sites along the genome sequence
sub getTags
{

   my @tags     =split /-/, $enzyme;
   my %collector=();
  
   die ("Bad reader\n") if (!defined $fr);

   while ((my $head, my $body) = $fr->getRecord())
   {
 
       while($body=~/($tags[0]|$tags[1])/g){

             my $pos=pos($body);
             $collector{$pos}->{header}="$head"."_".$pos;
             $collector{$pos}->{body}  =$&;
       }
   }
   return %collector;

}


