library("picante")
library("matrixStats")
library("amap")
library("bigmemory")
library("biganalytics")
args <- commandArgs(trailingOnly = TRUE)



print(args)
prefix<-args[1];
simTime<-as.numeric(args[2]);
alpha<-as.numeric(args[3]);
sampleRate<-as.numeric(args[4]);
start<-as.numeric(args[5]);
over<-as.numeric(args[6]);
groupNum<-as.numeric(args[7]);
OUT<-args[8]


IN<-prefix
## The input format 
## markerID 0 1 0 1 
geno=read.table(IN)
ncol<-dim(geno)[2]
nrow<-dim(geno)[1]
n_subr<-over-start+1
n_subc<-ceiling((ncol-1)*sampleRate)
cnt<-matrix(0,n_subr,ncol-1)
dt<-geno[,2:ncol]
ID<-geno[,1]
print(dim(dt))


### Cluster markers into different groups using bigmemory packages 
group=matrix(0,nrow,simTime)
for(i in 1:simTime){
    showTxt<-paste("running the",i,"simultion", sep=" ");
    print(showTxt)
    order<-sample(ncol-1)
    dt_sim<-dt[,order[1:n_subc]]
    a<-bigkmeans(as.matrix(dt_sim),groupNum,iter.max=1000, nstart=1, dist="euclid")
	group[,i]<-as.vector(a$cluster)
}

rownames(group)<-ID
print("output the result now ...")
write.table(group,file=OUT,quote = FALSE, sep = "\t",row.names = TRUE, col.names = FALSE)

### Output the group for each permutations
file.remove(OUT)
markerCNT<-dim(group)[1]
simTime<-dim(group)[2]
for(i in seq(1,simTime,1)){
	a<-group[,i]
	for(j in seq(1,groupNum,1)){
		pos<-which(a==j);
		if(length(pos)>1){
			write.table(as.matrix(t(pos)), file = OUT, sep = ",", append=TRUE, col.names=FALSE, row.names=FALSE)
		}
	}
}
