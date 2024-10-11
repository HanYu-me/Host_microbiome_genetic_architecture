library(data.table)
library(pheatmap)
library(stringr)

# get bin data ------------------------------------------------------------
faidx=data.frame(fread("Data//ZY300.fa.fai"))
faidx=faidx[faidx[,1]%in%paste0("Chr",1:24),]
faidx[,1]=as.numeric(sub("Chr","",faidx[,1]))
row.names(faidx)=faidx[,1]
faidx=faidx[order(faidx[,1]),]
chr_start=c(0)
for(j in 2:24){
  chr_start=c(chr_start,sum(faidx[1:(j-1),2]))
}
chr_start=c(chr_start,sum(faidx[,2]))

pwd="WAS_summary/GWAS_filter_summary/root50_NA/top_sug_snp"
files=list.files(pwd)
bin=1000000
seq_bin=data.frame(matrix(ncol=length(seq(1,(sum(faidx[,2])+10000000),bin)),nrow=length(files)))
names(seq_bin)=seq(1,(sum(faidx[,2])+10000000),bin)
pb=txtProgressBar(style=3)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  if(nrow(data)==0){
    next
  }
  for(j in 2:24){
    data[data[,1]==j,3]=as.numeric(data[data[,1]==j,3])+as.numeric(sum(faidx[1:(j-1),2]))
  }
  seq_bin[i,unique(findInterval(data[,3],as.numeric(names(seq_bin))))]=1
}
row.names(seq_bin)=str_split_fixed(files,"_",2)[,1]
root_gwa_sum=seq_bin

pwd="WAS_summary/GWAS_filter_summary/leaf50_NA/top_sug_snp"
files=list.files(pwd)
seq_bin=data.frame(matrix(ncol=length(seq(1,(sum(faidx[,2])+10000000),bin)),nrow=length(files)))
names(seq_bin)=seq(1,(sum(faidx[,2])+10000000),bin)
pb=txtProgressBar(style=3)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  if(nrow(data)==0){
    next
  }
  for(j in 2:24){
    data[data[,1]==j,3]=as.numeric(data[data[,1]==j,3])+as.numeric(sum(faidx[1:(j-1),2]))
  }
  seq_bin[i,unique(findInterval(data[,3],as.numeric(names(seq_bin))))]=1
}
row.names(seq_bin)=str_split_fixed(files,"_",2)[,1]
leaf_gwa_sum=seq_bin



files=list.files("Data//phe_h2//")
files=files[grep("hsq",files)]
phe_h2=c()
phe_num=c()
for(i in files){
  data=data.frame(fread(paste0("Data//phe_h2//",i),fill=T))
  phe_h2=c(phe_h2,data[4,2])
  phe_num=c(phe_num,data[10,2])
}

save(leaf_gwa_sum,root_gwa_sum,file="Results//24.07.01_gwaresults.RData")

# summary lambad ----------------------------------------------------------
leaf_lambda=read.table("WAS_summary//GWAS_filter_summary//leaf50_NA//info.txt")
root_lambda=read.table("WAS_summary//GWAS_filter_summary//root50_NA//info.txt")

par(mfrow=c(2,1))
hist(leaf_lambda[,10],main="leaf lambda",breaks=seq(0,2,0.01))
abline(v=1.1,col="red")
hist(root_lambda[,10],main="root lambda",breaks=seq(0,2,0.01))
abline(v=1.1,col="red")

table(leaf_lambda[,10]>1.1)
table(root_lambda[,10]>1.1)

boxplot(leaf_lambda[leaf_lambda[,10]>1.1,5],
        leaf_lambda[leaf_lambda[,10]<1.1,5],
        main="leaf qtn lambda>1,1 or <1.1")

leaf_lambda[leaf_lambda[,10]<1.1,c(1,5)][order(leaf_lambda[leaf_lambda[,10]<1.1,5]),]

boxplot(root_lambda[root_lambda[,10]>1.1,5],
        root_lambda[root_lambda[,10]<1.1,5],
        main="root qtn lambda>1,1 or <1.1")
tail(root_lambda[root_lambda[,10]<1.1,c(1,5)][order(root_lambda[root_lambda[,10]<1.1,5]),],20)

# remanhattan -------------------------------------------------------------
load("Results//MWAS//24.06.06_MWAS.RData")
leaf_gwa_sum[is.na(leaf_gwa_sum)]=0
root_gwa_sum[is.na(root_gwa_sum)]=0


leaf_gwa_sum[leaf_lambda[leaf_lambda[,10]>1.1,1],]=0
root_gwa_sum[root_lambda[root_lambda[,10]>1.1,1],]=0

chr_start=c(0)
for(j in 2:24){
  chr_start=c(chr_start,sum(faidx[1:(j-1),2]))
}
chr_start=c(chr_start,sum(faidx[,2]))
chr_end=chr_start[2:25]
chr_start=chr_start[1:24]
chr_mean=(chr_start+chr_end)/2



par(mfrow=c(2,1))
leaf_data=apply(leaf_gwa_sum,2,sum)
plot(x=names(leaf_data),y=leaf_data,pch=20,cex=1.5,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(leaf_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated ASVs in leaf",ylim=c(1.5,max(leaf_data)))
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")


root_data=apply(root_gwa_sum,2,sum)
plot(x=names(root_data),y=root_data,pch=20,cex=1.5,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(root_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated ASVs in root",ylim=c(1.5,max(root_data)))
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")

library(ggvenn)
ggvenn(list(leaf=names(leaf_data[leaf_data>=1]),root=names(root_data[root_data>=1])))
