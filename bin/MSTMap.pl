#!/usr/bin/perl
###############################################################################
#   Copyright (C) 2014, Ocean University of China (OUC);
#   all rights reserved. Authored by: Jinzhuang Dou
#   thinkhighly@163.com
#
#   This Software is used for scaffolding pre-assemblies based on 
#   Rad-map experiment
#
#   Redistribution and use in source and binary forms, with 
#   or without modification, are permitted provided that the 
#   following conditions are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the 
#     documentation and/or other materials provided with the distribution
#
#   * Neither the name of the OUC nor the names of its contributors may be 
#     used to endorse or promote products derived from this software 
#     without specific prior written permission.
#
#############################################################################
use strict;
use warnings;
use Cwd 'abs_path';
use Time::localtime;
use FindBin qw($Bin); 
use lib "$Bin/../lib";
use ParseFasta;
use Foundation;



# Declare default values for configure
sub getParameters();
my $tf = new Foundation                 ;
my $PRG = $tf->getProgramInfo('name')   ;
my $REV="1.0"                           ;
my @DEPENDS=("Foundation","ParseFasta") ;
my $VERSION = "Version 1.00"            ;
my $binDir         =$Bin;
my $dataDir        ;
my $ConfigInput    ;
my $GenotypeInput  ;
my %configParemeter;

my $NumCpu        ;
my $SimulationTime;
my $SimualtionRate;
my $startIterative =1;
my $endIterative   =9;
my $Groupthreshold =0.1;
my $NumMarker     ;
my $NumSample     ;
my $Confidence    ; # confidence interval

getParameters();




my $mstPath=Filefind($binDir, "MSTMap.exe");
if($mstPath eq "NO"){
       print "Warning: no software MSTMap.exe installed!\n";
       exit(); 
}
my $mstGenotype   ="mst.genotype";
my $mstMap        ="mst.map"     ;
#.......................Iterative  assemble........................
our %dis=();
our %marker=();
our %tag=();
our %d=();
our $cnt=0;

my %mapp=();
my %diss=();
my %numMarker=();
my %marker_lg=();
my %LinkMap=();

#The input of this software is the markers informations that has been grouped already.
if($startIterative==1){
  system("cp  ./markerGroup.map ./mst.map");
  system("mkdir step0") unless (-e "step0");
}

my $iterative;
for( $iterative=$startIterative;$iterative<=$endIterative;$iterative++){

   %dis=();
   %marker=();
   %tag=();
   %d=();
   $cnt=0;

   %mapp=();
   %diss=();
   %numMarker=();
   %marker_lg=();
   %LinkMap=();
   $SimulationTime=20; 
   my $id=$iterative-1;

   if($iterative==0){
      system("perl  $Bin/GroupFast.pl   -i ./tmp.geno  -o step0  -s 10 -T 0.55  -c 100");
      system("perl  $Bin/GroupSort.pl   -i  ./step0   -o step0");  
      system("cp ./step0/mst.map ./mst.map");
   }
   if($iterative==1){
     system("cp ./step0/mst.map ./mst.map"); 
   }
   mapSta("start");
   print "$iterative\n";


   system("cp mst.map  ./step$id/mst.map.order");
   boundaryMarker("mst.map",$GenotypeInput);
   permutationGenotype("permutation","tmp.geno");
   multiSort("permutation",$SimulationTime);  
   blockFind();
   blockConnect("Map");
   mapSta("over");
   mycopy($iterative);
}


sub mycopy{
    my $time=$_[0];
    my $t=$time-1;
    system("rm -rf step$time");
    system("mkdir step$time") unless (-e "step$time");
    system("cp -r permutation  ./step$t");
    system("cp *connect*  ./step$t");
    system("cp mst.map  ./step$time");
    system("mv final.map mst.map");
}

sub mapSta{
  open F, "<mst.map";
  my $flag=$_[0];
  my $lgCNT=0;
  my $markerCNT=0;
  while(<F>){
      chomp;
      if($_=~/lg/){$lgCNT++;}
      if($_=~/cons|-/){$markerCNT++;}
  }
  close(F);

  if($lgCNT<=2){
    print "Now no more than 2 groups are detected, and AMMO will stop assemble the groups\n";
    exit(-1);
  }
  if($flag eq "start"){
    print "Iteration $iterative  ... \n";
    print "Started!  Total marker: [$markerCNT] and group number [$lgCNT]\n";
  }
  else{
    print "Finished! Total marker: [$markerCNT] and group number [$lgCNT]\n";
  }
}

sub Filefind {
    my $line=0;
    my $file_dir;
    system("which $_[1] 1>./temp  2>./temp.log");
    system("find $_[0] -name  $_[1] 1>>./temp 2>./temp.log");
    open F, "<./temp";
    while(<F>){
         chomp;
         $line++;
         $file_dir=$_;
    }
    system("rm ./temp*");
    if($line==0){return "NO";}
    else{return $file_dir;}

}

sub blockConnect{
    my $sign=$_[0];
    %dis=();                     
    open F, "<./lg.connect"; 
    while(<F>){
         chomp;
         my @line=split(/\s+/,$_);
         for(my $i=0;$i<@line;$i++){
                $marker{$line[$i]}++;
         }
    }
    close(F);

    open F, "<./lg.connect"; 
    while(<F>){
         chomp;
         my @line=split(/\s+/,$_);
         $dis{$line[0]}{$line[0]}=0;
         for(my $i=1;$i<@line;$i++){
               
                $dis{$line[$i]}{$line[0]}=1;
                $dis{$line[0]}{$line[$i]}=1;

         }
    }
    close(F);
    my @node =keys %dis;
    my @node1=keys %dis;
     
    print "Order the sub-groups...\n";
    my $node_num=@node;
    my $node_num1=@node1;
    my $start;

    open(OUT, ">tree");
    $cnt=0;
    while($cnt<$node_num){ 
  
       %d=();
       for(my $i=0;$i<$node_num;$i++){
           if(exists $tag{$node[$i]}){;}
           else{  
                 $start=$node[$i];          
                 last;
           }
       }
  
       print OUT "\nStart  from  $start  ";

       my @nearby=keys % {$dis{$start}};
       Find_tree($start,"permutations");       
       
    } 

   close(OUT);

   open F, "<tree";
   open OUT, ">tree.map";

   foreach my $w(keys %numMarker){
       if (exists $dis{$w}){;}
       else{
          print OUT "$w\n";
       }
   }

   my $len;
   my $start_node;
   my %link=();
   my $i;
   my $j;
   while(<F>){
         chomp;
         my @line=split(/\s+/,$_);
         my %num=();
         $len=scalar @line;
         if($len==3){
                print OUT "$line[2]\n";
         }
         elsif($len>3){
                my %link=();
                for($i=3;$i<$len;$i++){
                       my @temp=split(/->/,$line[$i]);
                       $link{$temp[0]}{$temp[1]}++;
                       $link{$temp[1]}{$temp[0]}++;
                       $num{$temp[0]}++;
                       $num{$temp[1]}++;
                      
                }   
        
         
                my @node=keys %num;

                # I Need to fix the bug with the situation [/a-b/a-c/a-d/]
                my %rd_marker=();
                my %num_temp=%num;
                for($i=0;$i<@node;$i++){
                       if($num{$node[$i]}>2){ 
                            $num_temp{$node[$i]}=0;
                            foreach my $w(keys %link){
                                    if(exists $link{$w}{$node[$i]} || exists $link{$node[$i]}{$w}){
                                              $num_temp{$w}--;
                                    }
                            }
                            
                       }
                }                    
                for($i=0;$i<@node;$i++){
                   if($num_temp{$node[$i]}<=0){ 
                            $start_node=$node[$i];
                            $rd_marker{$start_node}++;          
                            print OUT "$start_node\n";
                   }
                }
                for($i=0;$i<@node;$i++){
                   if($num_temp{$node[$i]}==1){
                            $start_node=$node[$i];
                            if(exists $rd_marker{$start_node}){;}
                            else{ 
                                   $rd_marker{$start_node}++;          
                                   print OUT "$start_node";
                                   $len=2;
                                   while($len>1){
                                         foreach  $i (keys %{$link{$start_node}}){
                      
                                                   if(exists $rd_marker{$i}){;}
                                                   else{
                                                            $len=$num_temp{$i};
                                                            print OUT "->$i";
                                                            $start_node=$i;
                                                            $rd_marker{$i}++;
                                                   }
                                         }
                                  }
                                  print OUT "\n";
                           }
                  }
               }
             
        }             
   }
   close(F);
   close(OUT);

   #Get the final map according mst.map,tree-map,marker.connect
   if($sign eq "Scaffold"){;}
   else{
   open F,   "<marker.connect.all";  # Note I have changed that 
   open OUT, ">final.map";
   $cnt=0;
   while(<F>){
         chomp;
         $cnt++;
   }
   close F;
   close OUT;
   print "total link $cnt\n";


   open F,   "<tree.map";
   open OUT, ">final.map";
   $cnt=0;
   my @line;
   while(<F>){
          chomp;
          $cnt++;
          @line=split(/->/,$_);
   }
   #Fix the situation when three are only two linkage group
   if($cnt==1){
               #Sort according to the map distance of adjacent markers
               my %s=();
               my %hh=(); 
               my $n11;
               my $n12;
               my $n21;
               my $n22;
               for(my $i=1;$i<=$SimulationTime;$i++){
                  open F, "<./permutation/$i.map";  
                  while(<F>){
                       chomp;
                       if($_=~/-/){
                                    my @line1=split(/\s+/,$_);
                                    $hh{$line1[0]}=$line1[1];
                       }   
                  }  
                  close F;
                  foreach my $w1(keys %hh){
                           foreach my $w2(keys %hh){
                                    if($i==1){$s{$w1}{$w2}=0;}
                                    else{
                                           $s{$w1}{$w2}=abs($hh{$w1}-$hh{$w2})+$s{$w1}{$w2};
                                    }
                           }
                  }       
              }
              my $m11=$mapp{$line[0]}{1};
              my $m12=$mapp{$line[0]}{$numMarker{$line[0]}};
              my $m21=$mapp{$line[1]}{1};
              my $m22=$mapp{$line[1]}{$numMarker{$line[1]}};
              $n11=$s{$m11}{$m21};
              $n12=$s{$m11}{$m22};
              $n21=$s{$m12}{$m21};
              $n22=$s{$m12}{$m22};
              print OUT "group lg$cnt\n;BEGINOFGROUP\n";
              if($n11<$n12 && $n11<$n21 && $n11<$n12 ){

                              for($i=$numMarker{$line[0]};$i>=1;$i--){
                                  print OUT "$mapp{$line[0]}{$i}  $diss{$line[0]}{$i}\n";  
                              }
                              for($i=1;$i<=$numMarker{$line[1]};$i++){
                                  print OUT "$mapp{$line[1]}{$i}  $diss{$line[1]}{$i}\n";  
                              } 
              }
              elsif($n12<$n11 && $n12<$n21 && $n12<$n22 ){

                              for($i=1;$i<=$numMarker{$line[1]};$i++){
                                  print OUT "$mapp{$line[1]}{$i}  $diss{$line[1]}{$i}\n";  
                              } 
                              for($i=1;$i<=$numMarker{$line[0]};$i++){
                                  print OUT "$mapp{$line[0]}{$i}  $diss{$line[0]}{$i}\n";  
                              } 
              }
              elsif($n21<$n11 && $n21<$n12 && $n21<$n22 ){
                              for($i=1;$i<=$numMarker{$line[0]};$i++){
                                  print OUT "$mapp{$line[0]}{$i}  $diss{$line[0]}{$i}\n";  
                              } 

                              for($i=1;$i<=$numMarker{$line[1]};$i++){
                                  print OUT "$mapp{$line[1]}{$i}  $diss{$line[1]}{$i}\n";  
                              } 
                             
             }
             elsif($n22<$n11 && $n22<$n12 && $n22<$n21 ){
                              for($i=1;$i<=$numMarker{$line[1]};$i++){
                                  print OUT "$mapp{$line[0]}{$i}  $diss{$line[0]}{$i}\n";  
                              } 

                              for($i=$numMarker{$line[0]};$i>=1;$i--){
                                  print OUT "$mapp{$line[0]}{$i}  $diss{$line[0]}{$i}\n";  
                              }
        
             }
             print OUT ";ENDOFGROUP\n";
             close OUT;
  }
  else{

   #Add groups that do not have any linkage information with other groups
   
   open F,   "<tree.map";
   open OUT, ">final.map";
   $cnt=0;
   my @tmp=keys %mapp;
   while(<F>){
     chomp;
     $cnt++;
     if($_=~/>/){
       print OUT "group lg$cnt\n;BEGINOFGROUP\n";
       my @line=split(/->/,$_);
       my @con=split(/\s+/,$LinkMap{$line[0]}{$line[1]});

       if($con[0] eq $mapp{$line[0]}{1}){
              for($i=$numMarker{$line[0]};$i>=1;$i--){
                    print OUT "$mapp{$line[0]}{$i}  1 \n";  
              }
       }
       else{
              for($i=1;$i<=$numMarker{$line[0]};$i++){
                    print OUT "$mapp{$line[0]}{$i}  2\n";  
              } 
       }

       for($j=1;$j<@line;$j++){
             my @con=split(/\s+/,$LinkMap{$line[$j-1]}{$line[$j]});
             if($con[1] ne $mapp{$line[$j]}{1}){
                 for($i=$numMarker{$line[$j]};$i>=1;$i--){                   ### Fix me
                    print OUT "$mapp{$line[$j]}{$i}  3\n";  

#                    if($mapp{$line[$j]}{$i} =~/3700719/){
#                            print "Error \n"; print "$line[$j]\n";
                           
#                            for($i=$numMarker{$line[$j]};$i>=1;$i--){  print "$mapp{$line[$j]}{$i}  3\n"; }
#                             sleep(50);

#                     }
                 }
             }
             else{
                 for($i=1;$i<=$numMarker{$line[$j]};$i++){                   ### Fix me
                    print OUT "$mapp{$line[$j]}{$i}  4\n";  
                 } 
             } 
       } 
       print OUT ";ENDOFGROUP\n";

     }
     else{
       print OUT "group lg$cnt\n;BEGINOFGROUP\n";

       if($_=~/HASH/){;}
       else{
           for($i=1;$i<=$numMarker{$_};$i++){
                $diss{$_}{$i}="NA";                    
         
                print OUT "$mapp{$_}{$i}  $diss{$_}{$i}  5\n";
           }
           print OUT ";ENDOFGROUP\n";
       }
     }
   }
   close OUT;
 }
 }
 system("cp final.map mst.map");
}

sub Find_tree{
    my $ss=shift;
    my $Tag=shift;
    my @nearby=keys % {$dis{$ss}};
    $d{$ss}{$ss}=0;

    if(exists $tag{$ss}){;}
    else{
          $tag{$ss}++;
          $cnt++;
          for(my $i=0;$i<@nearby;$i++){
             if(exists $d{$nearby[$i]}{$ss}||$d{$ss}{$nearby[$i]}){;}
             else{  
                       $d{$ss}{$nearby[$i]}=$dis{$ss}{$nearby[$i]};
                       $d{$nearby[$i]}{$ss}=$dis{$ss}{$nearby[$i]};
                       if($Tag eq "group"){
                          	 print OUT "$ss\t$nearby[$i] ";
                       }
                       else{			   
                                 print OUT "$ss->$nearby[$i] ";
                       }
                       Find_tree($nearby[$i],$Tag);
             }
         }
    }
}
    

sub getBlock{
=header
  my $bin=1000;
  for(my $i=1;$i<=$markerNum;$i=$i+$bin){
   $start=$i;
   $over=$start+$bin-1;
   if($over>=$markerNum){$over=$markerNum;}
   $cmd="Rscript /mnt/projects/douj/ssmp/DouCode/AMMO/bin/markerSort.r tmp 20 0.8 0.99 $start $over ./permutation/tmp.$start-$over.link";
   system($cmd);
  }
=cut
  system("cat ./permutation/*.link > ./permutation/marker.link");
  my %hash=();
  open F, "<./permutation/marker.link"; 
  open O, ">./marker.connect";  
  while(<F>){
       chomp;
       my @line=split(/\s+|:/,$_);
       if(exists $hash{"$line[0]$line[1]"} || exists $hash{"$line[1]$line[0]"}){;}
       else{
        $hash{"$line[0]$line[1]"}++;
        print O "$marker_lg{$line[0]}\t$marker_lg{$line[1]}\t$line[0]\t$line[1]\t$line[2]\n";
       }
  }
  close O;
  close F;
}



sub blockFind{
     my %ConnectNum=();
     my %ConnectDis=();
     for(my $i=1;$i<=$SimulationTime;$i++){
        my $last="";
        my @line;
        my $last_dis="";
        open F, "<./permutation/$i".".map";   
        while(<F>){
             chomp;
	           if($_=~/-|consens/){
                   @line=split(/\s+/,$_);
                   if($last ne "" && $last_dis ne ""){                         
  		                     $ConnectNum{$line[0]}{$last}++;
                           $ConnectNum{$last}{$line[0]}++; 
                           if(exists $ConnectDis{$line[0]}{$last}){;}
                           else{
                                         $ConnectDis{$line[0]}{$last}=0;
                                         $ConnectDis{$last}{$line[0]}=0;
                           }
                           $ConnectDis{$line[0]}{$last}=abs($line[1]-$last_dis)+$ConnectDis{$line[0]}{$last};
                           $ConnectDis{$last}{$line[0]}=$ConnectDis{$line[0]}{$last};
		               }
	                 else{;} 
                   $last=$line[0];
                   $last_dis=$line[1];
	          }
    	}
        close(F);
    }  
 
    open O1, ">marker.connect";
    open OUT1, ">marker.connect.all";
   

    foreach my $i(keys %ConnectNum){

         $LinkMap{$marker_lg{$i}}{$marker_lg{$i}}++;
	       foreach my $j(keys %ConnectNum){
		     if(exists $ConnectNum{$i}{$j}){    
            my $threshold=$Confidence*$SimulationTime;
            my $left= $mapp{$marker_lg{$j}}{1};
            my $right=$mapp{$marker_lg{$j}}{$numMarker{$marker_lg{$j}}};
           
            if($ConnectNum{$i}{$j}>$threshold/2){
                    print OUT1 "$marker_lg{$i}\t$marker_lg{$j}\t$i\t$j\t$ConnectNum{$i}{$j}\n";
            }
            if($left ne $right && $numMarker{$marker_lg{$j}}<10){
              if(exists $ConnectNum{$i}{$left} && exists $ConnectNum{$i}{$right}){
                if($ConnectNum{$i}{$left}>$ConnectNum{$i}{$right}){
                   $ConnectNum{$i}{$left}=$ConnectNum{$i}{$left}+$ConnectNum{$i}{$right};
                   $ConnectNum{$i}{$right}=0;
                   $ConnectNum{$left}{$i}=$ConnectNum{$i}{$left};
                   $ConnectNum{$right}{$i}=0;
                }
                else{
                   $ConnectNum{$i}{$right}=$ConnectNum{$i}{$right}+$ConnectNum{$i}{$left};
                   $ConnectNum{$i}{$left}=0;
                   $ConnectNum{$right}{$i}=$ConnectNum{$i}{$right};
                   $ConnectNum{$left}{$i}=0;
                }
              }
            }     

            if($ConnectNum{$i}{$j}>$threshold || $ConnectNum{$j}{$i}>$threshold){
              if(exists $LinkMap{$marker_lg{$i}}{$marker_lg{$j}}){;}
              else{ 
                  $LinkMap{$marker_lg{$i}}{$marker_lg{$j}}="$i $j $ConnectNum{$i}{$j} $ConnectDis{$i}{$j}";
                  $LinkMap{$marker_lg{$j}}{$marker_lg{$i}}="$j $i $ConnectNum{$i}{$j} $ConnectDis{$i}{$j}";
              }
           }
          }
        }
    }

    
    foreach my $i(keys %LinkMap){
	   foreach my $j(keys %{$LinkMap{$i}}){
		       if($i ne $j){
                print O1 "$i\t$j\t$LinkMap{$i}{$j}\n";
		       }
	   }
    }
    close O1;
    removeMulitiNode();
}

sub removeMulitiNode{
    system("sort -k3,3 -k6,6n marker.connect > test");
    open F, "<test";
    my %lgTable=();
    my %hashTable=();
    my %groupCnt=();
    my @line;
    while(<F>){
             chomp;
             @line=split(/\s+/,$_);
             if(exists $hashTable{$line[0]}{$line[2]}){
                
             }
             else{
                $groupCnt{$line[0]}++;
                $hashTable{$line[0]}{$line[2]}++;
             }         
    }
    close F;
    my %index=();
    %hashTable=();

    open F, "<test";
    while(<F>){
             chomp;
             @line=split(/\s+/,$_);
             if(exists $index{$line[2]} && $groupCnt{$line[0]}==2 ){
                   $hashTable{$line[2]}{$line[3]}++;
                   $hashTable{$line[3]}{$line[2]}++;
             }
             else{
                   $index{$line[2]}++;
             }
    }
    close F;
    open F, "<test";
    open O, ">marker.connect.revise";
    while(<F>){
             chomp;
             @line=split(/\s+/,$_);
             if(exists $hashTable{$line[2]}{$line[3]} || exists $hashTable{$line[3]}{$line[2]}){;}
             else{
                     print O "$_\n";
                     $lgTable{$line[0]}{$line[1]}++;
                     $lgTable{$line[1]}{$line[0]}++;
             }
    }
    close F;
    close O;

    system("mv marker.connect.revise marker.connect");
    system("rm test");

    my($i,$j);
    open O, ">lg.connect";
    foreach $i(keys %lgTable){
       print O "$i";
       foreach  $j(keys %{$lgTable{$i}}){
           if($i ne $j){
                 print O "\t$j";
           }
        }
        print O "\n";
    }
    close O;
}

sub createRandNum
{
    my ($MaxNum,$MaxCount) = @_;
    my $i = 0;
    my %rand = ();
    while (1)
    {
                  my $no = int(rand($MaxNum))+1;
                  if ($no>$MaxNum){$no=$MaxNum;}
                  if (!$rand{$no})
                  {
                          $rand{$no}= 1;
                          $i++;                         
                  }
                  last if ($i >= $MaxCount);                    
    }
    my @randnum = keys % rand;
    return @randnum;
}


sub ReorderMSTmap{

     #my $NumSample=164;
     system("rm -rf group_split");
     open(F, "<$GenotypeInput");                
     my %table=(); 
     while(<F>){
             chomp;
             my @line=split(/\s+/,$_);
             $table{$line[0]}=$_;
     }
     close(F);
     open(F, "<mst.map");        # Change this 
     system("mkdir group_split") unless (-e "group_split");
     my $lg=0;
     my $lg_all=0;
     my %num={};
     my @file;
     my %hash={};
     while(<F>){
         chomp;
         if($_=~/lg/){
                $lg++;
                $num{$lg}=0;
         }
         elsif($_=~/consen|-/){
                $num{$lg}++;
                @file=split(/\s+/,$_);
                $hash{$lg}{$num{$lg}}=$file[0];
         }
     }
     close(F);

     $SimulationTime=$lg;
     for($lg=1;$lg<=$SimulationTime;$lg++){
             
             open O, ">./group_split/$lg.genotype";
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
             print O "number_of_loci $num{$lg}\n";
             print O "number_of_individual $NumSample\n";
             print O "ID";
             for(my $i=1;$i<=$NumSample;$i++){print O " P$i";}
             print O "\n";
             for(my $i=1;$i<=$num{$lg};$i++){print O "$table{$hash{$lg}{$i}}\n";}
             close(O);
     }
     
     system("rm mst.map");
     multiSort("group_split",$SimulationTime);
     system("cat ./group_split/*.map > ./mst.map"); 
     print "Finish group markers................\n";
}

sub mycat{
    #Since there are some strange thing that may make number of markers 
    open O, ">>mst.map";
    open O1, ">LOG";
    my $lg_num=shift @_;
    for(my $i=1;$i<=$lg_num;$i++){
            
            open F, "<./group_split/$i.genotype";
            my $rd=0;
            my $marker_num=0;
            while(<F>){
                   chomp;
                   if($_=~/number_of_loci/){ 
                           my @line=split(/\s+/,$_);
                           $marker_num=$line[1];
                   }
            }
            close F;

            open F, "<./group_split/$i.map";
            while(<F>){
                   chomp;
                   if($_=~/-|consen/){ $rd++;}
            }
            close F;
            print O1 "$rd <--------------------> $marker_num\n";
            if($rd eq $marker_num){
                 open F, "<./group_split/$i.map";
                 while(<F>){
                    chomp;
                    print O "$_\n";
                 }
                 close F;
            }           
     }
     close O;
     close O1;
}

sub buildGroup{

     system("rm -rf group_split");
     open(F, "<$GenotypeInput");                  
     my %table=(); 
     while(<F>){
             chomp;
             my @line=split(/\s+/,$_);
             $table{$line[0]}=$_;
     }
     close(F);

     open(F, "<mst.map");     # Change this 
     system("mkdir group_split") unless (-e "group_split");
     my $lg=0;
     my $lg_all=0;
     while(<F>){
         chomp;
         if($_=~/lg/){
             my %hash=();
             my @line=split(/\s+/,$_);
             my $lg_num=0;
             for(my $i=0;$i<@line;$i++){ 
                  if(exists $hash{$line[$i]}){;}
                  else{
                          $hash{$line[$i]}++;
                          $lg_num++;
                          $lg_all++;
                  }
             }
             $lg++;
             open O, ">./group_split/$lg.genotype";
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
             print O "number_of_loci $lg_num\n";
             print O "number_of_individual $NumSample\n";
             print O "ID";
             for(my $i=1;$i<=$NumSample;$i++){print O " P$i";}
             print O "\n";
             foreach my $w(keys %hash){
                   print O "$table{$w}\n";
             }
             close(O);
        }
     }
     $SimulationTime=$lg;
     multiSort("group_split");
     mycat($lg);
     print "Finish group markers................\n";
}  



sub boundaryMarker{
    
    #Cut the boundary of lg for next iterative map construction

    my $inputFile1       =  shift @_;
    my $inputFile2       =  shift @_;
    open F, "<$inputFile1";
    my $idLG     =0;
    my $markerLG =0;
    $NumMarker=0;
    while(<F>){
             chomp;
             my @line=split(/\s+|\t/,$_);
             if($_=~/lg/){
                       $idLG++;
                       $markerLG=0;
             }
             elsif($_=~/-|consens/){ 
                       #Need to be careful of the tag 
                       $markerLG++;
                       $marker_lg{$line[0]}=$idLG;
                       $mapp{$idLG}{$markerLG}=$line[0];
                       #$mapp{$markerLG}{$idLG}=$line[0];
                       $diss{$idLG}{$markerLG}=$line[1];
                       #$diss{$markerLG}{$idLG}=$line[1];
                       $numMarker{$idLG}=$markerLG    ;
             }
    }
    close F;

    my %tmpList=();
    open F, "<$inputFile2";
    while(<F>){
             chomp;
             my @line=split(/\s+/,$_);
             $tmpList{$line[0]}++;
    }
    close F;

    my %rd=();
    my $boundaryMarkerCnt=0;
    for(my $i=1;$i<=$idLG;$i++){
          if(exists $numMarker{$i} && $numMarker{$i}>0){

              my $start=0;
              my $over=0;
              my $getFlag=0;
              for(my $j=1;$j<=$numMarker{$i};$j++){
                if(exists $tmpList{$mapp{$i}{$j}}){
                   if($getFlag==0){
                        $start=$j;
                        $getFlag=1;
                   }
                   $over=$j;
                }
              }

              if($start>0 && $over>0 && $over>$start){
                $rd{$mapp{$i}{$start}}++;
                $rd{$mapp{$i}{$over}}++;
                $boundaryMarkerCnt=$boundaryMarkerCnt+2;
              }
              if($start>0 && $over>0 && $over==$start){
                $rd{$mapp{$i}{$over}}++;
                 $boundaryMarkerCnt=$boundaryMarkerCnt+1;
              }
              if($over==0 ){;}
          }
    }
    open F, "<$inputFile2";
    open O, ">tmp.geno";
    my %char2num=();
    $char2num{"A"}=0;
    $char2num{"B"}=1;
    $NumMarker=0;
    while(<F>){
             chomp;
             my @line=split(/\s+/,$_);

             if(exists $rd{$line[0]}){
                print O "$line[0]";
                $NumMarker++;
                for(my $i=1;$i<@line;$i++){
                   print O "\t$line[$i]";
                  #print O "\t$char2num{$line[$i]}";
                }
                print O "\n";
             }
    }
    close F;
    close O;
    print "Finish getting the boundary markers $boundaryMarkerCnt  [$idLG]\n";
}





sub permutationGenotype{

    my $Tag =shift @_;
    my $file=shift @_;
    
    my $directory="permutation";
    if($Tag eq "group"){ $directory="group";}  
    system("mkdir $directory") unless (-e $directory);
    my $SimualtionRate=0.8;
    my $RandSample  =int($NumSample*$SimualtionRate);
    my $lg2Threshold=2**(-$Groupthreshold);
    print "Sub-sampling parameters: simTime [$SimulationTime] sampleNum [$NumSample] confidence [$SimualtionRate] groupThreshold [$lg2Threshold]\n";


    for(my $i=1;$i<=$SimulationTime;$i++){
           my @RandSampleID=createRandNum($NumSample,$RandSample);
           open  O,">./$directory/$i.genotype";
           print O "population_type DH\n";
           print O "population_name LG\n";
           print O "distance_function kosambi\n";
           print O "cut_off_p_value $lg2Threshold\n";
           print O "no_map_dist 15.0\n";
           print O "no_map_size 0\n";
           print O "missing_threshold 1.00\n";
           print O "estimation_before_clustering no\n";
           print O "detect_bad_data yes\n";
           print O "objective_function COUNT\n";
           print O "number_of_loci $NumMarker\n";
           print O "number_of_individual $RandSample\n";
           print O "ID";
           for(my $i=1;$i<=$RandSample;$i++){print O " P$i";}
           print O "\n";
           open F, "<$file";
           while(<F>){
                chomp;
                my @line=split(/\s+/,$_);
                print O "$line[0]";
                for(my $j=1;$j<@line;$j++){
                           if(grep {$_==$j} @RandSampleID){
                                  print O " $line[$j]";
                           }
                }
                print O "\n";
          }
          close F;
          close O;
   }
}


sub multiSort{

    my $Tag=shift @_; 
    my $SimulationTime=shift @_;    
    my $directory="permutation";
    if($Tag eq "group"){ $directory="group";}
    elsif($Tag eq "group_split"){ $directory="group_split";}
      
    my $time=int($SimulationTime/$NumCpu)+1;
    my $command="";
    for(my $i=1;$i<=$time;$i++){
          $command="";
          my $j=$NumCpu*($i-1)+1;
          if($j<=$SimulationTime){
              $command="($mstPath ./$directory/$j.genotype  ./$directory/$j.map > out.log)";
              for($j=$NumCpu*($i-1)+2;$j<=$NumCpu*$i;$j++){
                  if($j<=$SimulationTime){
                        $command=$command."&($mstPath ./$directory/$j.genotype ./$directory/$j.map > out.log )";
                  }
              }
              system($command);
          }   
    }  
}

sub getParameters(){

    my @specOpts =() ;

    my $result = $tf->TIGR_GetOptions(
           "p=i"                   => \$NumCpu,
           "s=i"                   => \$SimulationTime,
           "t=i"                   => \$NumSample,
           "c=f"                   => \$Confidence,
           "g=s"                   => \$GenotypeInput,
           "a=i"                   => \$startIterative,
           "b=i"                   => \$endIterative,
           "T=f"                   => \$Groupthreshold

   );

   if (@specOpts) {
     print "$0: Argument required.\n";
     usage();
   }
   usage() if  (!$result);
   $tf  -> printUsageInfoAndExit() if (!$result);
}


sub usage {

    print STDERR <<EOQ;

    perl MSTMap.pl  -p  -s -t -c -g -a -b -l
    -p NumCpu
    -s SimulationTime
    -t NumSample
    -c Confidence
    -g GenotypeInput
    -a First IterativeTime 
    -b Last  IterativeTime
    -l GroupThreshold

EOQ
exit(0);
}

