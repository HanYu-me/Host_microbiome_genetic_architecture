library(data.table)
library(pheatmap)
library(stringr)

# get coordinate and data-----------------------------------------------------------
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
seq_bin=data.frame(matrix(ncol=length(seq(1,(sum(faidx[,2])+10000000),10000000)),nrow=length(files)))
names(seq_bin)=seq(1,(sum(faidx[,2])+10000000),10000000)
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
seq_bin=data.frame(matrix(ncol=length(seq(1,(sum(faidx[,2])+10000000),10000000)),nrow=length(files)))
names(seq_bin)=seq(1,(sum(faidx[,2])+10000000),10000000)
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

pwd="WAS_summary/Phe_summary/top_sug_snp"
files=list.files(pwd)
seq_bin=data.frame(matrix(ncol=length(seq(1,(sum(faidx[,2])+10000000),10000000)),nrow=length(files)))
names(seq_bin)=seq(1,(sum(faidx[,2])+10000000),10000000)
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
phe_gwa_sum=seq_bin


# heatmap for intersection ------------------------------------------------
load("Results//MWAS//24.06.06_MWAS.RData")
load("Results//24.06.11_gwaresults.RData")


pwd="WAS_summary/Phe_summary/sub_mlma"
files=list.files(pwd)
files=setdiff(files,paste0(c("TMV","RBSH","KC"),"_sub.txt"))
seq_bin=data.frame(matrix(ncol=length(seq(1,(sum(faidx[,2])+10000000),10000000)),nrow=length(files)))names(seq_bin)=seq(1,(sum(faidx[,2])+10000000),10000000)
pb=txtProgressBar(style=3)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  print(files[i])
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  for(j in 2:24){
    data[data[,1]==j,3]=as.numeric(data[data[,1]==j,3])+as.numeric(sum(faidx[1:(j-1),2]))
  }
  if(i==1){
    plot(data$bp,-log10(data$p),col=ifelse(data$Chr%%2==1,"tomato","skyblue"),pch=20,ylim=c(3,15))
  }else{
    points(data$bp,-log10(data$p),col=ifelse(data$Chr%%2==1,"tomato","skyblue"),pch=20)
  }
}
