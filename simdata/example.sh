##################################
#find ./ -type f -not -name 'example*' -delete

# Generate the simulate datasets 

#perl ../bin/2bRADsim.pl -r example.ref -d 0.5 -c [AGTCN]{30}AC[AGTCN]{5}CTCC[ATGCN]{30}-[ATGCN]{30}GGAG[AGTCN]{5}GT[ATGCN]{30}  -l 71-71 -o  simtest

# Run AMMO 
#perl ../bin/AMMO.pl  -p example.conf -S templateFilter  


#perl ../bin/AMMO.pl  -p example.conf -S markerType  

#perl ../bin/AMMO.pl  -p example.conf -S markerGroup

#perl ../bin/AMMO.pl  -p example.conf -S markerSort -g 2bRADtyping.ctg.type  -s 20 -t 100 -c 0.6 -a 1 -b 3 -l 0 


