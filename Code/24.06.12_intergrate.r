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

pwd="WAS_summary/Phe_summary/top_snp"
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
phe_gwa_sum=seq_bin
leaf_gwa_sum[is.na(leaf_gwa_sum)]=0
root_gwa_sum[is.na(root_gwa_sum)]=0
phe_gwa_sum[is.na(phe_gwa_sum)]=0
save(phe_gwa_sum,leaf_gwa_sum,root_gwa_sum,file="Results//24.6.11_gwaresults.RData")
# all together ------------------------------------------------------------
load("Results//MWAS//24.06.06_MWAS.RData")
load("Results//24.6.11_gwaresults.RData")
leaf_gwa_sum[is.na(leaf_gwa_sum)]=0
root_gwa_sum[is.na(root_gwa_sum)]=0
phe_gwa_sum[is.na(phe_gwa_sum)]=0


# manhattan plot ----------------------------------------------------------
chr_start=c(0)
for(j in 2:24){
  chr_start=c(chr_start,sum(faidx[1:(j-1),2]))
}
chr_start=c(chr_start,sum(faidx[,2]))
chr_end=chr_start[2:25]
chr_start=chr_start[1:24]
chr_mean=(chr_start+chr_end)/2

leaf_data=apply(leaf_gwa_sum,2,sum)
plot(x=names(leaf_data),y=leaf_data,pch=20,cex=1.5,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(leaf_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated ASVs in leaf")
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")



root_data=apply(root_gwa_sum,2,sum)
plot(x=names(root_data),y=root_data,pch=20,cex=1.5,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(root_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated ASVs in root")
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")

plot(leaf_data,root_data,pch=20,col="skyblue",frame.plot = F,
     xlab="Number of associated ASVs in leaf",ylab="Number of associated ASVs in root")
abline(a=0,b=1,col="grey",lwd=2,lty="dashed")
abline(lm(root_data~leaf_data+1),col="cyan2",lwd=2)
cor.test(leaf_data,root_data)
summary(lm(root_data~leaf_data+1))

plot(leaf_data,root_data,pch=20,col="skyblue",frame.plot = F,
     xlab="Number of associated ASVs in leaf",ylab="Number of associated ASVs in root")
abline(lm(root_data~leaf_data+1),col="red",lwd=2)
confint(lm(root_data~leaf_data+1),level=0.95)
abline(a=2.842061,b=0.910948,col="grey",lwd=2,lty="dashed")
abline(a=3.435992,b=1.213019,col="grey",lwd=2,lty="dashed")
table(leaf_data==0&root_data!=0)
table(leaf_data>1&root_data<=1)
table(leaf_data==0&root_data==0)


plot(leaf_data,root_data,pch=20,col="skyblue",frame.plot = F,
     xlab="Number of associated ASVs in leaf",ylab="Number of associated ASVs in root")
abline(h=quantile(root_data,0.5))
abline(v=quantile(leaf_data,0.5))
plot(zscore(leaf_data),zscore(root_data),pch=20,col="skyblue",frame.plot = F,
     xlab="Number of associated ASVs in leaf",ylab="Number of associated ASVs in root")
abline(h=quantile(root_data,0.5))
abline(v=quantile(leaf_data,0.5))
abline(lm(zscore(root_data)~zscore(leaf_data)),col="red",lwd=2)
summary(lm(zscore(root_data)~zscore(leaf_data)))

library(ggplot2)
#bin <- hexbin(leaf_data, root_data, xbins = 50)
#plot(bin, main = 'hexbin')
ggplot(data.frame(x =root_data , y = leaf_data), aes(x = x, y = y)) +
  geom_hex(binwidth = c(2.4, 1.2))+scale_fill_gradient2(low="lightskyblue1",high="royalblue3",transform="log1p")+ 
  coord_fixed()+ylim(c(0,10))+xlim(c(0,35))+theme_classic()+
  xlab("Number of associated ASVs in root")+ylab("Number of associated ASVs in leaf")+
  geom_smooth(method = "lm", se=T, color="cyan2", formula = y ~ x)


ggplot(data.frame(x =zscore(root_data) , y = zscore(leaf_data)), aes(x = x, y = y)) +
  geom_hex(binwidth = c(0.8, 0.8))+scale_fill_gradient2(low="lightskyblue1",high="royalblue3",transform="log1p")+ 
  coord_fixed()+ylim(c(-3,5))+xlim(c(-3,5))+theme_classic()+
  xlab("Number of associated ASVs in root")+ylab("Number of associated ASVs in leaf")+
  geom_smooth(method = "lm", se=T, color="cyan2", formula = y ~ x)

# get asv annotation and analysis the similarity --------------------------
leaf_feature=data.frame(fread("Data//FeatureTable_leaf.tsv"))
root_feature=data.frame(fread("Data//FeatureTable_root.tsv"))
root_feature=data.frame(root_feature$X.OTU.ID,root_feature$taxonomy,str_split_fixed(root_feature$taxonomy,";",7),stringsAsFactors = F)
leaf_feature=data.frame(leaf_feature$X.OTU.ID,leaf_feature$taxonomy,str_split_fixed(leaf_feature$taxonomy,";",7),stringsAsFactors = F)
names(leaf_feature)=c("ASV","taxonomy","k","p","c","o","f","g","s")
names(root_feature)=c("ASV","taxonomy","k","p","c","o","f","g","s")
leaf_feature$final_tax=""
for(i in 1:nrow(leaf_feature)){
  if(leaf_feature[i,"g"]!=" Unclassified"){
    leaf_feature[i,"final_tax"]=leaf_feature[i,"g"]
  }else{
    tax=leaf_feature[i,3:6][leaf_feature[i,3:6]!=" Unclassified"]
    leaf_feature[i,"final_tax"]=tax[length(tax)]
  }
}

sort(table(leaf_feature$final_tax))
root_feature$final_tax=""
for(i in 1:nrow(root_feature)){
  if(root_feature[i,"g"]!=" Unclassified"){
    root_feature[i,"final_tax"]=root_feature[i,"g"]
  }else{
    tax=root_feature[i,3:6][root_feature[i,3:6]!=" Unclassified"]
    root_feature[i,"final_tax"]=tax[length(tax)]
  }
}
row.names(root_feature)=root_feature$ASV
row.names(leaf_feature)=leaf_feature$ASV
# check chr 14 peak -------------------------------------------------------
library(ggvenn)
library(patchwork)

which(leaf_data==max(leaf_data))
which(root_data==max(root_data))
peak_14_leaf=row.names(leaf_gwa_sum)[leaf_gwa_sum[,"2310000001"]==1]
peak_14_root=row.names(root_gwa_sum)[root_gwa_sum[,"2310000001"]==1]

leaf_feature[peak_14_leaf,"final_tax"]
root_feature[peak_14_root,"final_tax"]
for(i in c("k","p","c","o","f","g","s","final_tax")){
  
}
p1=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"p"],root=root_feature[peak_14_root,"p"]))  
p2=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"c"],root=root_feature[peak_14_root,"c"]))  
p3=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"o"],root=root_feature[peak_14_root,"o"]))  
p4=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"f"],root=root_feature[peak_14_root,"f"]))  
p5=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"g"],root=root_feature[peak_14_root,"g"]))  
p6=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"s"],root=root_feature[peak_14_root,"s"]))  
p7=ggvenn(list(leaf=leaf_feature[peak_14_leaf,"final_tax"],root=root_feature[peak_14_root,"final_tax"]))  
(p1|p2)/(p3|p4)/(p5|p6)


pheatmap(leaf_gwa_sum[peak_14_leaf,],cluster_rows = F,cluster_cols = F)
pheatmap(root_gwa_sum[peak_14_root,],cluster_rows = F,cluster_cols = F)


# subgenome genetic bias --------------------------------------------------
library(readxl)
subgenome_info <- data.frame(read_excel("Data/anno/NG-A62700_STables-rev2.xlsx",sheet = "Table-S8"))[,1:4]
subgenome_info=subgenome_info[c(-1),]
names(subgenome_info)=c("chr","start","end","subgenome")
subgenome_info[,1]=as.numeric(sub("Chr","",subgenome_info[,1]))
subgenome_info[,2]=as.numeric(subgenome_info[,2])
subgenome_info[,3]=as.numeric(subgenome_info[,3])
head(subgenome_info)
for(j in 2:24){
  subgenome_info[subgenome_info[,1]==j,3]=as.numeric(subgenome_info[subgenome_info[,1]==j,3])+as.numeric(sum(faidx[1:(j-1),2]))
  subgenome_info[subgenome_info[,1]==j,2]=as.numeric(subgenome_info[subgenome_info[,1]==j,2])+as.numeric(sum(faidx[1:(j-1),2]))
}
bin_sub_info=c()
for(i in as.numeric(names(leaf_gwa_sum))){
  bin_sub_info=c(bin_sub_info,subgenome_info[subgenome_info$start < i & subgenome_info$end > i,4])
}
S_gwa_leaf=apply(leaf_gwa_sum[,bin_sub_info=="S"],2,sum)
T_gwa_leaf=apply(leaf_gwa_sum[,bin_sub_info=="T"],2,sum)
S_gwa_root=apply(root_gwa_sum[,bin_sub_info=="S"],2,sum)
T_gwa_root=apply(root_gwa_sum[,bin_sub_info=="T"],2,sum)

boxplot(S_gwa_leaf,T_gwa_leaf,S_gwa_root,T_gwa_root,
        names=c("Leaf Ssub","Leaf Tsub","Root Ssub","Root Tsub"),
        col=c("tomato","orange","skyblue","#CCEBC5"),
        frame.plot=F,boxwex=0.5,ylab="Number of associated ASVs",ylim=c(0,10))
t.test(S_gwa_leaf,T_gwa_leaf)
t.test(S_gwa_root,T_gwa_root)
# analysis of genetic similarity PVE------------------------------------------
plot(leaf_data,root_data,pch=20,col="skyblue",frame.plot = F,
     xlab="Number of associated ASVs in leaf",ylab="Number of associated ASVs in root")
abline(h=quantile(root_data,0.5))
abline(v=quantile(leaf_data,0.5))
quantile(root_data[root_data>0],0.5)
quantile(leaf_data[leaf_data>0],0.5)
regions=names(leaf_data)[which(leaf_data>1&root_data>1)]
###get_common_snp_set
snp_pos=data.frame(fread(file.choose(),header=F)) #load .map file
for(j in 2:24){
  snp_pos[snp_pos[,1]==j,4]=as.numeric(snp_pos[snp_pos[,1]==j,4])+as.numeric(sum(faidx[1:(j-1),2]))
}
snp_index=findInterval(snp_pos[,4],as.numeric(names(leaf_data)))
comon_snp=snp_pos[,2][snp_index%in%which(leaf_data>1&root_data>4)]
uncommon_snp=setdiff(snp_pos[,2],comon_snp)
write.table(comon_snp,file="Data//PVE//common_snp.txt",row.names = F,col.names = F,quote=F)
write.table(uncommon_snp,file="Data//PVE//uncommon_snp.txt",row.names = F,col.names = F,quote=F)

plot(x=regions,y=rep(1,length(regions)))
plot(x=as.numeric(names(leaf_data)),y=rep(1,length(names(leaf_data))))
points(x=regions,y=rep(1.2,length(regions)))
###get grm
grmfile=data.frame(fread("Data//PVE//uncommon_snp.grm.gz"),stringsAsFactors = F)
grm_data=matrix(ncol=max(grmfile[,1]),nrow=max(grmfile[,1]))
for(i in 1:nrow(grmfile)){
  a=as.numeric(grmfile[i,1:4])
  grm_data[a[1],a[2]]=a[4] 
}
grm_data[upper.tri(grm_data)] <- t(grm_data)[upper.tri(grm_data)]
id=data.frame(fread("Data//PVE//uncommon_snp.grm.id"),stringsAsFactors = F)
row.names(grm_data)=id[,1]
uncommon_snp_grm=grm_data

grmfile=data.frame(fread("Data//PVE//common_snp.grm.gz"),stringsAsFactors = F)
grm_data=matrix(ncol=max(grmfile[,1]),nrow=max(grmfile[,1]))
for(i in 1:nrow(grmfile)){
  a=as.numeric(grmfile[i,1:4])
  grm_data[a[1],a[2]]=a[4] 
}
grm_data[upper.tri(grm_data)] <- t(grm_data)[upper.tri(grm_data)]
id=data.frame(fread("Data//PVE//common_snp.grm.id"),stringsAsFactors = F)
row.names(grm_data)=id[,1]
common_snp_grm=grm_data


### root PVE 
library(hglm)
root_phe=data.frame(fread("Data//phe_filtered/root_50_nosick_NA.txt"))
row.names(root_phe)=root_phe[,1]
pves_root=data.frame(matrix(ncol=3,nrow=length(3:ncol(root_phe))))
pves_root[,1]=names(root_phe)[3:ncol(root_phe)]
for(i in 3:ncol(root_phe)){
  setTxtProgressBar(pb, i/ncol(root_phe))
  ids=intersect(row.names(common_snp_grm),root_phe[,1])
  index=row.names(common_snp_grm)%in%ids
  phe=root_phe[ids,i]
  Z1=c.z.hglm(common_snp_grm[index,index])
  Z2=c.z.hglm(uncommon_snp_grm[index,index])
  na_index=!is.na(phe)
  model=hglm(X=matrix(rep(1,length(ids))[na_index]),y=phe[na_index],Z=cbind(Z1[na_index,],Z2[na_index,]),RandC=c(ncol(Z1),ncol(Z1)))
  pves_root[(i-2),2:3]=model$varRanef/(sum(model$varRanef)+model$varFix)
}

leaf_phe=data.frame(fread("Data//phe_filtered/leaf_50_nosick_NA.txt"))
row.names(leaf_phe)=leaf_phe[,1]
pves_leaf=data.frame(matrix(ncol=3,nrow=length(3:ncol(leaf_phe))))
pves_leaf[,1]=names(leaf_phe)[3:ncol(leaf_phe)]
for(i in 3:ncol(leaf_phe)){
  setTxtProgressBar(pb, i/ncol(leaf_phe))
  ids=intersect(row.names(common_snp_grm),leaf_phe[,1])
  index=row.names(common_snp_grm)%in%ids
  phe=leaf_phe[ids,i]
  Z1=c.z.hglm(common_snp_grm[index,index])
  Z2=c.z.hglm(uncommon_snp_grm[index,index])
  na_index=!is.na(phe)
  model=hglm(X=matrix(rep(1,length(ids))[na_index]),y=phe[na_index],Z=cbind(Z1[na_index,],Z2[na_index,]),RandC=c(ncol(Z1),ncol(Z1)))
  pves_leaf[(i-2),2:3]=model$varRanef/(sum(model$varRanef)+model$varFix)
}
save(pves_root,pves_leaf,file="Results//PVE//24.06.17_common_snp_PVE.RData")

par(mfrow=c(2,1))
hist(pves_leaf[,2]/(pves_leaf[,2]+pves_leaf[,3]),xlab="Ratio",main="leaf")
hist(pves_root[,2]/(pves_root[,2]+pves_root[,3]),xlab="Ratio",main="root")
pves_root$sum=pves_root[,2]+pves_root[,3]
data=t(pves_root[order(pves_root$sum,pves_root[,2],decreasing=T),2:3])
barplot(data,col=c("tomato","skyblue"),border=NA)

PVE_ratio=data.frame(c(pves_leaf[,2]/(pves_leaf[,2]+pves_leaf[,3]),pves_root[,2]/(pves_root[,2]+pves_root[,3])))
PVE_ratio$class=c(rep("Leaf",nrow(pves_leaf)),rep("Root",nrow(pves_root)))
names(PVE_ratio)[1]="value"
library(ggplot2)
library(gghalves) 
ggplot(PVE_ratio,aes(x=class,y=value,fill=class))+
  geom_half_violin(position = position_nudge(x=0.25),side = "r",width=0.8,color=NA)+
  geom_boxplot(alpha=1,width=0.4,,outlier.color = "black",outlier.size=1)+
  scale_fill_manual(values=c(Leaf="skyblue",Root="tomato"))+
  theme_bw()+ylim(c(0,1))+ ylab("PVE ratio")+xlab("")+
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 


# gene function -----------------------------------------------------------
library(AnnotationHub)
library(clusterProfiler)
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
source("Code//Def_function.r")
package = "AnnotationHub"
oldcache = path.expand(rappdirs::user_cache_dir(appname=package))
setAnnotationHubOption("CACHE", oldcache)
ah <- AnnotationHub() # vpn was needed
query(ah, "Nicotiana")
tob_orgdb <- ah[["AH114234"]]
GO_result <- data.frame(fread("Data//anno//function.stat.xls"))
names(GO_result)
row.names(GO_result)=gsub("(.*)\\.1","\\1",GO_result$GeneID)
source("Code//Def_function.r")

# get gff
gff=data.frame(fread("Data//anno//ZY300.gff3"))
gff=gff[gff$V3=="gene",]
gff$gene=gsub("(.*)\\.1","\\1",sub("ID=","",str_split_fixed(gff$V9,";",2)[,1]))
genes=gff[gff$V1=="Chr14" & gff$V4>(2310000000-2302807202) & gff$V4<(2320000000-2302807202),"gene"]
GO_list=get_go_fun(genes)
res_GO <- enrichGO(GO_list,OrgDb=tob_orgdb,keyType = 'GO',ont='ALL',
                   pvalueCutoff=0.05, qvalueCutoff = 0.05)
dotplot(res_GO,split="ONTOLOGY",showCategory=20,font.size=7)+facet_grid(ONTOLOGY~., scale="free")


# a narrow region for find genes ------------------------------------------
pwd="WAS_summary/GWAS_filter_summary/root50_NA/top_sug_snp"
files=list.files(pwd)
bin=100000
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
root_gwa_sum2=seq_bin

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
leaf_gwa_sum2=seq_bin

root_gwa_sum2=root_gwa_sum2[,as.numeric(names(root_gwa_sum2))<2427589413 & as.numeric(names(root_gwa_sum2))>2302807202]
leaf_gwa_sum2=leaf_gwa_sum2[,as.numeric(names(leaf_gwa_sum2))<2427589413 & as.numeric(names(leaf_gwa_sum2))>2302807202]
leaf_gwa_sum2[is.na(leaf_gwa_sum2)]=0
root_gwa_sum2[is.na(root_gwa_sum2)]=0

par(mfrow=c(2,1))
leaf_data=apply(leaf_gwa_sum2,2,sum)
plot(x=names(leaf_data),y=leaf_data,pch=20,cex=1.5,frame.plot=F,
     col="tomato",xlab="Chromosome 14",ylab="Number of associated ASVs in leaf")
root_data=apply(root_gwa_sum2,2,sum)
plot(x=names(root_data),y=root_data,pch=20,cex=1.5,frame.plot=F,
     col="skyblue",xlab="Chromosome 14",ylab="Number of associated ASVs in root")

which(leaf_data==max(leaf_data))
which(root_data==max(root_data))

root_gwa_sum2=root_gwa_sum2[,as.numeric(names(root_gwa_sum2))<2330589413 & as.numeric(names(root_gwa_sum2))>2302807202]
leaf_gwa_sum2=leaf_gwa_sum2[,as.numeric(names(leaf_gwa_sum2))<2330589413 & as.numeric(names(leaf_gwa_sum2))>2302807202]
pheatmap(leaf_gwa_sum2,cluster_rows = F,cluster_cols = F)
pheatmap(root_gwa_sum2,cluster_rows = F,cluster_cols = F)

genes=gff[gff$V1=="Chr14" & gff$V4>(2315800001-2302807202) & gff$V4<(2316800001-2302807202),"gene"]
GO_list=get_go_fun(genes)
res_GO <- enrichGO(GO_list,OrgDb=tob_orgdb,keyType = 'GO',ont='ALL',
                   pvalueCutoff=0.05, qvalueCutoff = 0.05)
dotplot(res_GO,split="ONTOLOGY",showCategory=15,font.size=7)+
  facet_grid(ONTOLOGY~., scale="free")+
  scale_fill_gradient(low="tomato",high="skyblue")


# use less threshold ------------------------------------------------------
sub_snp=snp_pos[snp_pos[,4]>2315800001 & snp_pos[,4]<2316800001,]
pwd="WAS_summary/GWAS_filter_summary/root50_NA/sub_mlma"
files=list.files(pwd)
pb=txtProgressBar(style=3)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  row.names(data)=data[,2]
  sub_snp[,(i+4)]=-log10(data[sub_snp$V2,"p"])
}
names(sub_snp)[5:ncol(sub_snp)]=str_split_fixed(files,"_",2)[,1]
peak_14_root=sub_snp[5:ncol(peak_14_root)]
row.names(peak_14_root)=sub_snp[,2]

sub_snp=snp_pos[snp_pos[,4]>2315800001 & snp_pos[,4]<2316800001,]
pwd="WAS_summary/GWAS_filter_summary/leaf50_NA/sub_mlma"
files=list.files(pwd)
pb=txtProgressBar(style=3)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  row.names(data)=data[,2]
  sub_snp[,(i+4)]=-log10(data[sub_snp$V2,"p"])
}
names(sub_snp)[5:ncol(sub_snp)]=str_split_fixed(files,"_",2)[,1]
peak_14_leaf=sub_snp[5:ncol(peak_14_leaf)]
row.names(peak_14_leaf)=sub_snp[,2]

peak_14_leaf[is.na(peak_14_leaf)]=0
pheatmap(peak_14_leaf,cluster_rows = F,cluster_cols = F)

peak_14_root[is.na(peak_14_root)]=0
pheatmap(peak_14_root,cluster_rows = F,cluster_cols = F)



# class genetic similarity analyisis --------------------------------------
head(leaf_feature)

table(leaf_feature[row.names(leaf_gwa_sum),"p"])
table(leaf_feature[row.names(leaf_gwa_sum),"o"])

table(root_feature[row.names(root_gwa_sum),"p"])
table(root_feature[row.names(root_gwa_sum),"o"])

library(ggvenn)
ggvenn(list(leaf=unique(leaf_feature[row.names(leaf_gwa_sum),"p"]),root=unique(root_feature[row.names(root_gwa_sum),"p"])))  
ggvenn(list(leaf=unique(leaf_feature[row.names(leaf_gwa_sum),"o"]),root=unique(root_feature[row.names(root_gwa_sum),"o"])))  

inter_o=intersect(unique(root_feature[row.names(root_gwa_sum),"o"]),unique(leaf_feature[row.names(leaf_gwa_sum),"o"]))

big_mat=data.frame(matrix(ncol=length(row.names(root_gwa_sum)),nrow=length(row.names(leaf_gwa_sum))))
row.names(big_mat)=row.names(leaf_gwa_sum)
names(big_mat)=row.names(root_gwa_sum)
for(i in row.names(leaf_gwa_sum)){
  for(j in row.names(root_gwa_sum)){
    ids=as.character(intersect(root_phe[,1],leaf_phe[,1]))
    sub=cor.test(zscore(leaf_phe[ids,i]),zscore(root_phe[ids,j]))
    big_mat[i,j]=sub$estimate
    # if(sub$p.value<0.05){
    #   big_mat[i,j]=sub$estimate  
    # }else{
    #   big_mat[i,j]=NA
    # }
  }
}
pheatmap(big_mat,cluster_cols = F,cluster_rows = F)

asvs_leaf=c()
asvs_root=c()
num1=c()
num2=c()
for(i in 1:12){
  asvs_root=c(asvs_root,row.names(root_gwa_sum)[root_feature[row.names(root_gwa_sum),"o"]==inter_o[i]])
  num1=c(num1,length(row.names(root_gwa_sum)[root_feature[row.names(root_gwa_sum),"o"]==inter_o[i]]))
  asvs_leaf=c(asvs_leaf,row.names(leaf_gwa_sum)[leaf_feature[row.names(leaf_gwa_sum),"o"]==inter_o[i]])
  num2=c(num2,length(row.names(leaf_gwa_sum)[leaf_feature[row.names(leaf_gwa_sum),"o"]==inter_o[i]]))
}
annotation_row=data.frame(Class=rep(inter_o,num2))
row.names(annotation_row)=asvs_leaf
annotation_col=data.frame(Class=rep(inter_o,num1))
row.names(annotation_col)=asvs_root
pheatmap(big_mat[asvs_leaf,asvs_root],cluster_cols = T,cluster_rows = T,
         annotation_row=annotation_row,annotation_names_row = F,
         annotation_col=annotation_col,annotation_names_col = F,
         annotation_legend = T)

big_mat2=data.frame(matrix(ncol=length(asvs_root),nrow=length(asvs_root)))
row.names(big_mat2)=asvs_root
names(big_mat2)=asvs_root
for(i in asvs_root){
  for(j in asvs_root){
    sub=cor.test(zscore(root_phe[ids,i]),zscore(root_phe[ids,j]))
    big_mat2[i,j]=sub$estimate
    # if(sub$p.value<0.05){
    #   big_mat[i,j]=sub$estimate  
    # }else{
    #   big_mat[i,j]=NA
    # }
  }
}

big_mat3=data.frame(matrix(ncol=length(asvs_leaf),nrow=length(asvs_leaf)))
row.names(big_mat3)=asvs_leaf
names(big_mat3)=asvs_leaf
for(i in asvs_leaf){
  for(j in asvs_leaf){
    sub=cor.test(zscore(leaf_phe[ids,i]),zscore(leaf_phe[ids,j]))
    big_mat3[i,j]=sub$estimate
    # if(sub$p.value<0.05){
    #   big_mat[i,j]=sub$estimate  
    # }else{
    #   big_mat[i,j]=NA
    # }
  }
}
pheatmap(big_mat2,cluster_rows = F,cluster_cols = F)
pheatmap(big_mat3,cluster_rows = F,cluster_cols = F)
### genetic similarity
par(mfrow=c(3,4))
for(i in 1:12){
  asvs_root=row.names(root_gwa_sum)[root_feature[row.names(root_gwa_sum),"o"]==inter_o[i]]
  asvs_leaf=row.names(leaf_gwa_sum)[leaf_feature[row.names(leaf_gwa_sum),"o"]==inter_o[i]]
  a=apply(leaf_gwa_sum[asvs_leaf,],2,sum)
  b=apply(root_gwa_sum[asvs_root,],2,sum)
  plot(a,b,pch=20,col="skyblue",frame.plot = F,cex=3,
       xlab="Number of associated ASVs in leaf",ylab="Number of associated ASVs in root",main=inter_o[i])
}
