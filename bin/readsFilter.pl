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
parse_command_line();
open FIG,    "<$input";
open OUT,    ">$output";
$cnt=0;
#$base1="[ATGCN]{9}AC[ATGCN]{5}CTCC[ATGCN]{7}";
#$base2="[ATGCN]{7}GGAG[ATGCN]{5}GT[ATGCN]{9}";
print "$base1\t$base2\n";
$all=0;$high_qual=0;
print "$l\n";
while($word= <FIG>){
         chomp($word=$word);
         @file=split(/\s+/,$word);
         $cnt=$cnt+1;
         if($cnt%2==1){$name=substr($word,1);}
         elsif($cnt%2==0){
              $read=$word;
              $len=length $read;
              $tag=0;
              if($read=~/$base1/){
                             while($read=~/$base1/g ){ 
                                  $end=pos($read);
                                  $ss=substr($read,$end-$l,$l);
                                  $tag++;
                             }
                             print OUT ">$name\n$ss\n";
              }
              $read1=$read;
              if($tag==0 && $read1=~/$base2/){
                             while($read1=~/$base2/g ){ 
                                  $end=pos($read1);
                                  $ss=substr($read1,$end-$l,$l);
                                  $tag++;
                             }
                             print OUT ">$name\n$ss\n";
              }
         }
}
           
sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if       ($_ =~ /^-i$/)  { $input   = shift  @ARGV; }
        elsif    ($_ =~ /^-o$/)  { $output  = shift  @ARGV; }
        elsif    ($_ =~ /^-o1$/)  { $output1  = shift  @ARGV; }
        elsif    ($_ =~ /^-b1$/)  { $base1    = shift  @ARGV; }
        elsif    ($_ =~ /^-b2$/)  { $base2    = shift  @ARGV; }
        elsif    ($_ =~ /^-l$/)  { $l       = shift  @ARGV; }
        elsif    ($_ =~ /^-f$/)  { $f       = shift  @ARGV; }
        elsif    ($_ =~ /^-Q1$/) { $Q1      = shift  @ARGV; }
        elsif    ($_ =~ /^-Q2$/) { $Q2      = shift  @ARGV; }
        elsif    ($_ =~ /^-S$/)  { $S       = shift  @ARGV; }
        else {
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	}
        
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl reads_filter.pl -i file path  -o file path -o1 trash -b -f  -T  -Q1  -Q2 -N -[-h]
    i    :input file 
    o    :outputfile
    o1   :trash file
    b    :target restriction site for BsaXI:[ATGC]{9}AC[ATGC]{5}CTCC[ATGC]{7}
    l    :length of read
    f    :f=1 palindromic structure;f=2 no palindromic structure [2]
    Q1   :threshold for low quality score [20].
    Q2   :maximum no.of low-quality bases[10].  
    S    :discard reads with homopolymers >(Sxread length)[0.3].
    h    :display the help information.
EOQ
exit(0);
}

