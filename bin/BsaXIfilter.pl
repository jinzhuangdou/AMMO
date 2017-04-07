  ##!/usr/bin/perl -w
##########################################################
#-this progrome is used for high quality reads filtering-#
##########################################################
my $mode1="";
my $l=0;
my $f=2;
my $Q1=20;
my $Q2=10;
my $N=0.3; 
my $S=0.3;
$input=$ARGV[0];
open (OUT, " | bgzip > $input".".QC.fq.gz");
open (FIG, "gzip -dc $input".".filter.fq.gz | ");
$cnt=0;
$base1="[ATGCN]{9}AC[ATGCN]{5}CTCC[ATGCN]{7}";
$base2="[ATGCN]{7}GGAG[ATGCN]{5}GT[ATGCN]{9}";
print "$base1\n";
while (<FIG>){
         chomp;
         $cnt=$cnt+1;
         if($cnt%4==1){$name=$_;}
         if($cnt%4==2){$seq=$_;}
         if($cnt%4==3){$strand=$_;}
         if($cnt%4==0){
             $qual=$_;
             if($seq=~/$base1/ || $seq=~/$base2/){
                print OUT "$name\n$seq\n$strand\n$qual\n";
             }
         }
}
           

