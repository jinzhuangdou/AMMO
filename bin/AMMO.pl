#!/usr/bin/perl -w
###########################################################################
#   Copyright (C) 2015, Ocean University of China (OUC);
#   all rights reserved. Authored by: Jinzhuang Dou
#   
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
###########################################################################

# Program      : AMMO 
# Function     : Scaffolding(sca) assemblies(l) based on low-cost(l) HappySeq (p) technology.
# Version      : 1.1
# Created time : April 4th  2015
# Revised time : June  4    2015
########################################################################### 

use strict;
use warnings;
use Cwd 'abs_path';
use Time::localtime;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use ParseFasta;
use Foundation;



# Declare default values for configure

my $tf = new Foundation                 ;
my $PRG = $tf->getProgramInfo('name')   ;
my $REV="1.0"                           ;
my @DEPENDS=("Foundation","ParseFasta") ;
my $VERSION = "Version 1.00"            ;
my %global                              ;
my %default_value                       ;
my $stage   =0                          ;
my @process = ("all", "templateFilter", "markerType", "markerGroup",
               "markerSort","contigLink","clean");


setDefaults()                           ;

setParametersFromCommandLine()          ;

foreach my $w (keys %default_value){
    $global{$w} = $default_value{$w} if(!defined $global{$w});
}

for (my $i=0; $i<@process; $i++){
    $stage = $i if ( $global{"STAGE"} eq $process[$i]);
}
print "Now $stage: $global{STAGE}\n";


# Retain the DNA from the targeted genome sequence using the templatFile 
# from short reads or pre-assemblies.
if ($stage ==1 || $stage==0){
    my $result = fileCheck("INPUT_DIR", "TEMPLATE_FILE", "CTG_FILE");  
    printPrnLog("Loading the parameters successfully .\n") if $result;
    templeteFilter(); 
    printPrnLog("Marker denosing using the templetFile $global{TEMPLATE_FILE} successfully.\n");
}

# Type the markers according to the sequencing coverage 
if ($stage ==2 || $stage==0){
    markerType_ctg();
    printPrnLog("Marker typing successfully . \n");
}

# Let the contigs as the primary linkage group if groupMethod ==1. If groupMethod ==0, single 
# marker will be regarded as a single primary linkage group, the computaitonal complexity for 
# further markers assemblies will be very very huge. Hence, parameter groupMethod will recommended 
# 1 for the most case.  

if ($stage ==3 || $stage==0){
   markerGroup($global{"GROUP_METHOD"},"2bRADtyping");
   printPrnLog("Primary linkage group was constructed succesfully\n");
}

if ($stage ==4 || $stage==0){

   # Scaffolding the primary linkage group using the hierarchial algorithm toghether with the 
   # sampling technology.
   reGetParameter();
   markerSort(1,"2bRADtyping.ctg.type");
   printPrnLog("Secondary linkage group was constructed iteratively\n");

}

if ($stage ==5 || $stage==0){
   reGetParameter();
   output();
}

if ($stage ==6 || $stage==0){
   myclean();
   printPrnLog("Over! Thanks for your interests in using AMMO\n");
}



sub myclean{
    system("mkdir -p temporaryFile");
    system("rm test");
    #system("mv step* marker* lg*  *mst* permu* 2bRAD* consensus.fasta ctgUnmap.fasta  demo  ./temporaryFile");
}
sub reGetParameter{
    # Reget the value of parameters $global{MARKER_NUM} and $global{SAMPLE_SIZE} 
    $global{MARKER_NUM} = -1 ;
    $global{SAMPLE_SIZE}= 0  ;
    if(-e "2bRADtyping.ctg.type"){;}
    else{
       if(-e "2bRADtyping"){
          printPrnLog("No templeteFile available in the current folder, the original genotype file will be used.\n");
          system("cp 2bRADtyping 2bRADtyping.ctg.type");
        }
        else{
           printPrnLog("Error: No genotype file will be used.\n");
           exit(-1);
        }
    }
    open F, "<2bRADtyping.ctg.type";
    while(<F>){
          chomp;
          $global{MARKER_NUM}++;
          my @f=split(/\s+/,$_);
          my $len = @f;
          $global{SAMPLE_SIZE} = $len -1;
          last;
    }
    close F;  
}


sub output{

  my $cmd = "perl $Bin/disOutput.pl -g 2bRADtyping.ctg.type -c $global{CTG_FILE}  -m mst.map -p markerGroup.map -o $global{OUT_PREFIX} -n $global{SAMPLE_SIZE}";
  system($cmd);
}

sub markerSort{
  my $cmd = "perl $Bin/MSTMap.pl -g 2bRADtyping.ctg.type  -p $global{CPU_NUM} -s $global{SAMPLE_TIME} -t $global{SAMPLE_SIZE} -c $global{THRESHOLD} -a $global{ITERATIVE_START} -b $global{ITERATIVE_OVER} -T $global{GROUP_THRESHOLD}";
  print "$cmd\n";
  system($cmd);
}





sub markerGroup{

    my $soapPath       = "$Bin"."/soap" ;
    my $bwtbuilderPath = "$Bin"."/2bwt-builder" ;
    my $dataDir        = ".";
    my $templeteFile   = $global{TEMPLATE_FILE};
    my $inputDir       = $global{INPUT_DIR};
    my $ctgFile        = $global{CTG_FILE} ;


    my($groupMethod,$inputGenotype)=@_;
    my($groupNum)=0;
   
    $global{MARKER_NUM} =0;
    #$groupMethod=1;  #NOTE CHANGE HERE
    if($groupMethod==0){
        open F, "<$dataDir/$inputGenotype";
        open O, ">$dataDir/markerGroup.map";
        while(<F>){
                chomp;
                $groupNum++;
                my @f=split(/\s+/,$_);
                print O "lg$groupNum\n$f[0]\n";
        }
        close O;
        close F;
    }
    elsif($groupMethod==1){

        if(-e $ctgFile){

              system("$bwtbuilderPath $ctgFile") if($global{INDEX});
              system("$soapPath -a $dataDir/consensus.fasta -D  $ctgFile".".index -M 4 -r 0 -v 2 -o $dataDir/soapMap");
              system("rm $dataDir/*.index*") if($global{INDEX});


              my %hashTable    =();
              my %enzymeMapList=();
              my %enzymeList   =Get_ctg_ATGC("$dataDir/consensus.fasta");
              my %ctgList      =Get_ctg_ATGC("$ctgFile");

              open F, "<$dataDir/soapMap";
              while(<F>){
                chomp;
                my @f=split(/\s+/,$_);
                $hashTable{$f[7]}{$f[8]}=$f[0];
                $enzymeMapList{$f[0]}++;
              }         
              close F;
              system("rm $dataDir/soapMap");


			        open O, ">$dataDir/markerGroup.map";
              foreach my $w1(keys %hashTable){
                 $groupNum++;
                 print O "lg$w1\n";
                 foreach my $w2 (sort {$a<=>$b} keys %{$hashTable{$w1}}){
                      print O "$hashTable{$w1}{$w2}  $w2\n";
                 }
              }
              close O;


              open O, ">$dataDir/markerUnmap.fasta";
              foreach my $w1(keys %enzymeList){
                      if(exists $enzymeMapList{$w1}){;}
                      else{
                            print O ">$w1\n$enzymeList{$w1}\n";
                      }
              }
              close O;


              open O, ">$dataDir/ctgUnmap.fasta";
              foreach my $w1(keys %ctgList){
                      if(exists $hashTable{$w1}){;}
                      else{
                            print O ">$w1\n$ctgList{$w1}\n";
                      }
              }
              close O;

              open F, "<$inputGenotype";
              open O, ">$inputGenotype.ctg.type";
              while(<F>){
                chomp;
                my @f=split(/\s+/,$_);
                if(exists $enzymeMapList{$f[0]}){
                                 print O "$_\n";
                                 $global{MARKER_NUM}++;
                }
              }         
              close F;
              close O;

        }
        else{
             print "Warning: if you want to use the contig sequence to build the\n";
             print "group informaiton, you must specify the pre-assemblies!\n";
             exit();
        }
    }

    elsif($groupMethod==2){
          #Finish this part one day!
    }
    system("cp $dataDir/markerGroup.map $dataDir/mst.map");


}

sub markerType_ctg{
  my $soapPath       = "$Bin"."/soap" ; 
	my $inputDir          = $global{INPUT_DIR}; 
	my $minSample         = $global{MIN_TYPING} *$global{SAMPLE_SIZE};
	my $maxSample         = $global{MAX_TYPING} *$global{SAMPLE_SIZE};
  my $templeteFile=$global{TEMPLATE_FILE};
 	system("cat $inputDir/*filter* > $inputDir/temp.merge.fa");
  system("mkdir $inputDir/soap/");

  #... extract [AC.....CTCC] from ctgs ... 
  print "$templeteFile\n";
  open (IN, "<$templeteFile" );
	my $fr        = new ParseFasta(\*IN);
	my $enzyme="[AGTCN]{100}AC[AGTCN]{5}CTCC[ATGCN]{100}-[ATGCN]{100}GGAG[AGTCN]{5}GT[ATGCN]{100}";   # 
  #my $enzyme="[AGTCN]{100}CCGG[ATGCN]{100}-[ATGCN]{100}CC[AGTCN]{1}GG[ATGCN]{100}";   # 
	my @tags     =split /-/, $enzyme;
	my %collector=();
  my $line=0;
        die ("Bad reader\n") if (!defined $fr);
         while ((my $head, my $body) = $fr->getRecord()){
         $body =~ tr/[a-z]/[A-Z]/;
         my @tmp=split(/\s+|-/,$head);
         while($body=~/($tags[0]|$tags[1])/g){
              my $pos=pos($body);
              my $id="$tmp[0]"."-".$pos;
              $line++;
              if($line%1000==0){print "Extract $line tags from reference sequence\n";}
              $collector{$id}->{header}=$id;
              $collector{$id}->{body}  =$&;
	       }
  }

	open O, ">./consensus.fasta";
	foreach my $w(keys %collector){
		print O ">$collector{$w}->{header}\n$collector{$w}->{body}\n";
	}
	close O;
	system("cp ./consensus.fasta   $inputDir/soap/ctg.BsaXI.fa");
	system("$Bin/2bwt-builder   $inputDir/soap/ctg.BsaXI.fa"); 
  system("$soapPath -a  $inputDir/temp.merge.fa   -D  $inputDir/soap/ctg.BsaXI.fa.index  -M 4 -r 0 -p 40 -v 2 -o $inputDir/soap/soapMap");

	my %hashTable=();
	my %sampleID=();
  my $lineCNT=0;
	open F, "<$inputDir/soap/soapMap";
    while(<F>){
		  chomp;
        $lineCNT++;
        if($lineCNT%10000==0){print "$lineCNT reads have been finished!\n";} 
        my @file=split(/\s+|\t/,$_); 
        my @tmp=split(/-/,$file[0]); 
			  $hashTable{$file[7]}{$tmp[0]}++;
			  $sampleID{$tmp[0]}++;
  }
	close F;
    
  open O,  ">./2bRADtyping";
	my $missingCnt=0;
	my @sampleList=keys %sampleID;
	foreach my $w1(keys %hashTable){ 
			$missingCnt=0;
			if($w1=~/HASH/){;}
			else{
		     for(my $i=0;$i<@sampleList;$i++){
			     if(exists $hashTable{$w1}{$sampleList[$i]}){$missingCnt++;}
		     }
			 if($missingCnt>$minSample && $missingCnt<$maxSample){
				print O "$w1 ";
			    for(my $i=0;$i<@sampleList;$i++){
					if(exists $hashTable{$w1}{$sampleList[$i]}){print O " A";}
					else{
							print O " B";
					}
				}
			    print O "\n";
            }
		  }
	}
	close(O);

}


sub markerType{

    my $ustacksPath       = "$Bin"."/ustacks" ;
    my $inputDir          = $global{INPUT_DIR};
    my $minSample         = $global{MIN_TYPING} *$global{SAMPLE_SIZE};
    my $maxSample         = $global{MAX_TYPING} *$global{SAMPLE_SIZE};

    system("cat $inputDir/*filter* > $inputDir/temp.merge.fa");

    if (-e "$inputDir/stacks"){ 
        system("rm -r $inputDir/stacks");
    }
    system("mkdir $inputDir/stacks/");
	  system("$ustacksPath -t fasta -m 20  -M 2 -o  $inputDir/stacks -p 20 -f $inputDir/temp.merge.fa");
    
    


    #Split the merged file according to the fileID into each samples.
    my %hashTable=();
    my $consensusID=0;
    my %consensusRead=();
    my %sampleID=();

    my @file=();
    open F, "<$inputDir/stacks/temp.merge.tags.tsv";
    while(<F>){
         chomp;
         @file=split(/\s+/,$_);
         if($_=~/consensus/){
            
            $consensusID++;
            $consensusRead{$consensusID}=$file[6];
         }
         elsif($_=~/primary/){
              @file=split(/-/,$file[5]);
              $hashTable{$consensusID}{$file[0]}++;
              $sampleID{$file[0]}++;
         }
         elsif($_=~/secondary/){
         #Note I need to change the col ID if there is some wrong with the secondary ID.
              @file=split(/-/,$file[4]);
              $hashTable{$consensusID}{$file[0]}++;
              $sampleID{$file[0]}++;
         }
    }
    close(F);



    open O,  ">./2bRADtyping";
    open O1, ">./consensus.fasta";
    my $missingCnt=0;
    my @sampleList=keys %sampleID;

    foreach my $w1(keys %hashTable){
            
            $missingCnt=0;
            for(my $i=0;$i<@sampleList;$i++){
                if(exists $hashTable{$w1}{$sampleList[$i]}){$missingCnt++;}
            }

            if($missingCnt>$minSample && $missingCnt<$maxSample){
               print O "consensus$w1 ";
               print O1 ">consensus$w1\n$consensusRead{$w1}\n";
               for(my $i=0;$i<@sampleList;$i++){
                   if(exists $hashTable{$w1}{$sampleList[$i]}){print O " A";}
                   else{
                        print O " B";
                   }
               }
               print O "\n";  
            }   
    }
    close(O);
    close(O1);
}



sub templeteFilter{
    
    #print "$templeteFile  $inputDir\n";

    my $soapPath       = "$Bin"."/soap" ;
    my $bwtbuilderPath = "$Bin"."/2bwt-builder" ;
    my $templeteFile   = $global{TEMPLATE_FILE};
    my $inputDir       = $global{INPUT_DIR};
    $global{SAMPLE_SIZE}=0;
    $global{INDEX}=1;                ## Remember to delete this variable
    system("$bwtbuilderPath  $templeteFile") if $global{INDEX};
    opendir (DIR, $inputDir);
    my @dire = readdir DIR;
    my %fastaFile=();
    system("rm $inputDir/*filter*");
    system("rm $inputDir/*soap*");
    for(my $i=0;$i<@dire;$i++){
       if($dire[$i] eq "."||$dire[$i] eq ".." || 
               $dire[$i] =~/filter/ || $dire[$i] =~/stacks/ || $dire[$i] =~/merge/ || $dire[$i]=~/sh/){;}
       else{
           system("$soapPath -a $inputDir/$dire[$i] -D  $templeteFile".".index -M 4 -r 0 -p 40  -v 2 -o $inputDir/soapMap") ;
           open F, "<$inputDir/soapMap";
           open O, ">$inputDir/$dire[$i]".".filter";
           $global{SAMPLE_SIZE}++;
           my $line=0;
           print "$dire[$i]\n";
           while(<F>){
                chomp;
                if(length($_)>10){
                  my @file=split(/\t|\s+/,$_);
                  $line++;
                  print  O ">$i-$line\n$file[1]\n";
                }
           }
           close(F);
           close(O);
           system("rm $inputDir/soapMap");
       }
    }
    

}


sub Get_ctg_ATGC{


    my $file_dir=shift; 
    my %ctg_ATGC=();
    my $str="";
    open (my $fh, $file_dir) or die "Error, cannot read file $file_dir";
    open O, ">temp";
    while(<$fh>){
        chomp;
        if($_=~/>/){
              my @line=split(/\s+/,$_);
              my $ID=substr($line[0],1);
              print O "\n$ID ";
        }
        else{
              print O "$_";
        }                     
    }
    close ($fh);
    close O;

    open ($fh, "temp") or die "Error, cannot read file $file_dir";
    while(<$fh>){
        chomp;
        if($_ ne ""){
           my @line=split(/\s+/,$_);
           $ctg_ATGC{$line[0]}=$line[1];
        }          
    }
    close ($fh);
    # system("rm temp");
    return %ctg_ATGC;   
}


  
sub setParametersFromFile($@) {


    my $specFile  = shift @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "Error: PARAMETER_FILE '$specFile' not found.\n"  if (! -e "$specFile");

    print STDERR "\n";
    print STDERR "###\n";
    print STDERR "###  Reading options from '$specFile'\n";
    print STDERR "###\n";
    print STDERR "\n";

    open(F, "< $specFile") or die("Couldn't open '$specFile'", undef);

    while (<F>) {
        s/^\s+//;
        s/\s+$//;
        next if (m/^\s*\#/);
        next if (m/^\s*$/) ;
        my ($var, $val) = split /\s+/;
        undef $val if ($val eq "undef");
        setGlobal($var, $val);
    }
    close(F);
}


sub printPrnLog
{
    my $msg = shift ;
    $tf -> logLocal($msg,1);
    print STDERR $msg; 

}


sub setParametersFromCommandLine{

    my @specOpts =() ;

    my $result = $tf->TIGR_GetOptions(
           "p=s"                   => \$global{"PARAMETER_FILE"},
           "d=s"                   => \$global{"INPUT_DIR"},
           "t=s"                   => \$global{"TEMPLATE_FILE"},
           "g=s"                   => \$global{"GROUP_METHOD"},
           "c=s"                   => \$global{"CTG_FILE"},
           "o=s"                   => \$global{"OUT_PREFIX"},
           "S=s"                   => \$global{"STAGE"},
           "s=s"                   => \$global{"SAMPLE_TIME"},
           "T=n"                   => \$global{"GROUP_THRESHOLD"},
           "f=f"                   => \$global{"THRESHOLD"},
           "L=n"                   => \$global{"CPU_NUM"},
           "a=n"                   => \$global{"ITERATIVE_START"},
           "b=n"                   => \$global{"ITERATIVE_OVER"},
           "P=n"                   => \$global{"CPU_NUM"},
           "l=n"                   => \$global{"MIN_DEPTH"},
           "u=n"                   => \$global{"MIN_TYPING"},
           "v=n"                   => \$global{"MAX_TYPING"},
           "m=n"                   => \$global{"SOAP_MATCH"},
           "M=n"                   => \$global{"USTACKS_M"},
           "C=n"                   => \$global{"USTACKS_C"},
   );

   $tf  -> printUsageInfoAndExit() if (!$result || @specOpts);

   if(defined $global{"PARAMETER_FILE"}) {
      setParametersFromFile($global{"PARAMETER_FILE"});   
   }
   
}



sub fileCheck{

    for my $var (@_){
       print "$global{$var}\n";
      if (!defined $global{$var} || ! -e "$global{$var}") {
          printPrnLog("\n\nError: $var not defined or not found, pls check it!\n\n");
          $tf  -> printUsageInfoAndExit()
      }
      else{
          setGlobal($var, abs_path($global{$var}));
      }
    }
    return 1;
}








sub setDefaults {
    

    my $HELPTEXT = qq~
Scaffold assemblies based on low-cost HappySeq technology. 
Program : AMMO.pl
Author  : Jinzhuang Dou
Released: June 4th, 2015
Usage:    $PRG file -p configFile [options]
or        $PRG file -d inputDir -t templeteFile -c ctgFile -o outprefix  [options]

Main parameters:

  -p <s>     The configure file of all parameters setting required to running this software
  -d <s>     The input directory storing the raw 2b-RAD datasets generated from Happy experiments
  -t <s>     The templete sequence used to discard the noisy reads derived from other species.
  -g <n>     Whether you want to get a 2b-RAD physical map or scaffold the pre-assemblies generated from
             other softwares. 1 for scaffolding pre-assemblies and 0 for de novo physical map construction
             If later, please specify the contigs as an input.\
  -c <s>     reads files storing the pre-assembly for scaffolding.
  -o <s>     The prefix will be added to the filenames of outputting results.
  -S <s>     The stage which you want to run using AMMO, including all, templateFilter, markerType, 
             markerGroup, markerSort, and GapFill. [all]
  -s <n>     The maximal sampling times allowed for markers grouping and markers sorting.
  -T <f>     threshold over which pairs of contigs will be linked together under the specified sampling times. [0.6] 
  -n <n>     The maximal iterative times allowed for markers grouping and markers sorting.  
  -P <n>     Number of processors to use for parallel computation.
  -l <n>     Minimum depth allowed for markers typing 
  -u <f>     The maximal sampling times allowed for markers grouping and markers sorting [0.1].
  -v <f>     The threshold set to determine the linkage information between pairwise markers 
             for all sampling times [0.9].
  -h|help    print this help and exit;
  -V|version print the version and exit;
  -depend    print the program and database dependency list;
  -debug     set the debug <level> (0, non-debug by default); 


Advanced paramters

  -m <n>     match mode for reads mapping using soap software. 
             0: exact mathc only
             1: 1 mismatch match only
             2: 2 mismatch match only
             4: find the best hits
  -M <n>     Minimum depth of coverage required to creat a stack [3]
  -C <n>     Maximum distance allowed between stacks.
        
OUTPUT:  
~;

my $MOREHELP = qq~
Return Codes:   0 - on success, 1 - on failure.
~;



    $global{"HELP"}               = $HELPTEXT; 
    $global{"INPUT_DIR"}          = undef    ;
    $global{"TEMPLATE_FILE"}      = undef    ;
    $global{"GROUP_METHOD"}       = undef    ;
    $global{"CTG_FILE"}           = undef    ;
    $global{"OUT_PREFIX"}         = undef    ;
    $global{"STAGE"}              = undef    ;
    $global{"PARAMETER_FILE"}     = undef    ;
    $global{"SAMPLE_TIME"}        = undef    ;
    $global{"THRESHOLD"}          = undef    ;
    $global{"ITERATIVE_START"}    = undef    ;
    $global{"ITERATIVE_OVER"}     = undef    ;
    $global{"CPU_NUM"}            = undef    ;
    $global{"MIN_DEPTH"}          = undef    ;
    $global{"MIN_TYPING"}         = undef    ;
    $global{"MAX_TYPING"}         = undef    ;
    $global{"SOAP_MATCH"}         = undef    ;
    $global{"USTACKS_M"}          = undef    ;
    $global{"USTACKS_C"}          = undef    ;
    $global{"GROUP_THRESHOLD"}    = undef    ;


    $global{"INDEX"}              = undef      ; # Remember to remove these variables
    $global{"MAX_TYPING"}         = undef      ; #
    $global{"SAMPLE_SIZE"}        = undef      ; #


    $default_value{"HELP"}               = $HELPTEXT; 
    $default_value{"INPUT_DIR"}          = undef    ;
    $default_value{"TEMPLATE_FILE"}      = undef    ;
    $default_value{"GROUP_METHOD"}       = 1        ;
    $default_value{"CTG_FILE"}           = undef    ;
    $default_value{"OUT_PREFIX"}         = undef    ;
    $default_value{"STAGE"}              = "all"    ;
    $default_value{"PARAMETER_FILE"}     = undef    ;
    $default_value{"SAMPLE_TIME"}        = 20       ;
    $default_value{"THRESHOLD"}          = 12       ;
    $default_value{"ITERATIVE_START"}    = 1        ;
    $default_value{"ITERATIVE_OVER"}     = 10       ;
    $default_value{"CPU_NUM"}            = 1        ;
    $default_value{"MIN_DEPTH"}          = 2        ;
    $default_value{"MIN_TYPING"}         = 0.1      ;
    $default_value{"MAX_TYPING"}         = 0.9      ;
    $default_value{"SOAP_MATCH"}         = 4        ;
    $default_value{"USTACKS_M"}          = 1        ;
    $default_value{"USTACKS_C"}          = 2        ;
    $default_value{"GROUP_THRESHOLD"}           = 0        ;


    $default_value{"INDEX"}              = 1        ; # Remember to remove these variables
    $default_value{"MAX_TYPING"}         = 0.9      ; #
    $default_value{"SAMPLE_SIZE"}        = 100      ; #

    
    $tf->setHelpInfo($HELPTEXT.$MOREHELP);
    $tf->setUsageInfo($HELPTEXT);
    $tf->setVersionInfo($REV);
    $tf->addDependInfo(@DEPENDS);

    #printPrnLog("Set the default values of parameters\n");

}



sub getGlobal ($) {
    my $var = shift @_;
    if (!exists($global{$var})) {
      return undef;
    }
    return($global{$var});
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;

    #  If no value, set the field to undefined, the default for many of the options.
    $val = undef  if ($val eq "");
    $global{$var} = $val if(!defined $global{$var});
}
