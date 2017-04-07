$prefix=$ARGV[0];
open F, "<$prefix.mst.map";
open O, ">$prefix.mst.map.Kendall";
while(<F>){
  chomp;
  if($_=~/consensus/){
     @f=split(/-|\t|_|\s+/,$_);
     $len=@f;
     if(@f==4){
        if($f[1] eq $last){
           print O  "\t$f[2]";
        }
        else{
           sleep(1);
           print O  "\n$f[2]";
        }
        $last=$f[1];
     }
  }
}
close F;
close O;
$total=0;
$sum=0;
open F, "<$prefix.mst.map.Kendall";
while(<F>){
   chomp;
   @f=split(/-|\t|_|\s+/,$_);
   if(@f>2){
   $cnt_true=0;
   $cnt_total=0;
   $dif=0;
   $max=0;
   for($i=0;$i<@f;$i++){
    for($j=$i+1;$j<@f;$j++){
         $dif=abs($f[$j]-$f[$i]);
         if($dif>$max){$max=$dif;}
         if($f[$i]<$f[$j]){$cnt_true++;}
         $cnt_total++;
    }
   }
   $accuracy=$cnt_true/$cnt_total;
   if($accuracy<0.5){$accuracy=1-$accuracy;}
   if( $max<1000000){print "$accuracy\n";}
   $sum+=$accuracy;
   $total+=1;
   }
}
close F;
$accuracyAve=$sum/$total;
print "The kendall static is $accuracyAve\n";
