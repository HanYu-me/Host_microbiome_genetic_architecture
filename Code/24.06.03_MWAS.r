library(data.table)
library(pheatmap)
library(hglm)
install.packages("hglm",type="source")
install.packages("lme4",type="source")
source("/data/public/hanyu/Meta_tobacco/Code/Def_function.r")

# data format ----------------------------------------------------------------
phe_blup=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/phe_filtered/predict_traits.csv"))
phe_type=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/phe_filtered/Trait_in_research.csv"))
phe_blup=phe_blup[,!grepl("Cha",names(phe_blup))]
phe_blup=cbind(phe_blup$id,phe_blup[,names(phe_blup)%in%phe_type$abbrevation.in.pre])
phe_info=data.frame("Traits"=names(phe_blup),"Description"="","type"="","Num"=0)
for(i in 1:nrow(phe_info)){
    if(sum(phe_type$abbrevation.in.pre==phe_info[i,1])!=0){
        phe_info[i,2:4]=phe_type[phe_type$abbrevation.in.pre==phe_info[i,1],c(2,4,8)]
    }
}
phe_info=phe_info[-1,]
write.table(phe_info,file="/data/public/hanyu/Meta_tobacco/Data/lines_info/Phe_info.txt",row.names=F,col.names=T,quote=F)

leaf_all=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/phe_filtered/leaf_all_nosick_0.txt"))
leaf_50=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/phe_filtered/leaf_50_nosick_0.txt"))
root_all=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/phe_filtered/root_all_nosick_0.txt"))
root_50=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/phe_filtered/root_50_nosick_NA.txt"))
table(leaf_all$id%in%phe_blup[,1])
table(root_all$id%in%phe_blup[,1])

row.names(root_50)=root_50[,1]
row.names(leaf_50)=leaf_50[,1]
row.names(phe_blup)=phe_blup[,1]
phe_pre=phe_blup[,-1]
root_50=root_50[,c(-1,-2)]
leaf_50=leaf_50[,c(-1,-2)]
phe_s=cbind(row.names(phe_pre),row.names(phe_pre),phe_pre)
names(phe_s)[1:2]=c("family","id")
write.table(phe_s,file="/data/public/hanyu/Meta_tobacco/Data/phe_filtered/Phe_all_predicted.txt",row.names = F,col.names = T,quote = F)
pca=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/Final_0418_pruned_pca.eigenvec"))
ids=data.frame(fread("/data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/Final_0418_pruned.fam"))
row.names(pca)=as.character(ids[,1])
pca=pca[,3:7]
# MWAS ----------------------------------------------------------------
p_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
est_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
for(i in 1:ncol(phe_pre)){
    print(i)
    ids=intersect(row.names(root_50),row.names(phe_pre))
    for(j in 1:ncol(root_50)){
        a=lm(phe_pre[ids,i]~root_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5])
        a=data.frame(summary(a)[4])
        est_mat[i,j]=a[2,1]
        p_mat[i,j]=a[2,4]
    }
}
root_mwas=list("Beta"=est_mat,"P"=p_mat)

p_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
est_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
for(i in 1:ncol(phe_pre)){
    print(i)
    ids=intersect(row.names(leaf_50),row.names(phe_pre))
    for(j in 1:ncol(leaf_50)){
        a=lm(phe_pre[ids,i]~leaf_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5])
        a=data.frame(summary(a)[4])
        est_mat[i,j]=a[2,1]
        p_mat[i,j]=a[2,4]
    }
}
leaf_mwas=list("Beta"=est_mat,"P"=p_mat)

#data=p_mat
data=root_mwas$P
data=leaf_mwas$P
sub=t(apply(data,1,function(x){return(ifelse(x<0.05/length(x),1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.05,1,0))}))
row.names(sub)=names(phe_pre)
pheatmap(sub,cluster_rows = F,cluster_cols = F)

plot(zscore(phe_pre[ids,i]),zscore(leaf_50[ids,j]))
plot(phe_pre[ids,i],leaf_50[ids,j])
cor.test(zscore(phe_pre[ids,i]),zscore(leaf_50[ids,j]))
cor.test(phe_pre[ids,i],leaf_50[ids,j])


# mixed linear reguresion -------------------------------------------------
library(hglm)
library(lme4)
library(pheatmap)
library(lmtest)
pheatmap(ibs)
p_mat_1=matrix(nrow=ncol(phe_pre),ncol=100)
p_mat_2=matrix(nrow=ncol(phe_pre),ncol=100)
p_mat_3=matrix(nrow=ncol(phe_pre),ncol=100)
p_mat_4=matrix(nrow=ncol(phe_pre),ncol=100)
p_mat_5=matrix(nrow=ncol(phe_pre),ncol=100)
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(ibs),intersect(row.names(root_50),row.names(phe_pre)))
  for(j in 1:100){
    data=cbind.data.frame(phe_pre[ids,i],root_50[ids,j],pca[ids,1],pca[ids,2],pca[ids,3],pca[ids,4],pca[ids,5])
    names(data)=paste0("V",1:7)
    
    index=row.names(ibs)%in%ids
    na_index=!is.na(root_50[ids,j])
    Z=c.z.hglm(ibs[index,index])
    hglm_res=hglm(X=matrix(as.numeric(root_50[ids,j][na_index])),y=as.numeric(phe_pre[ids,i][na_index]),Z=Z[na_index,],calc.like = T,bigRR=T)
    hglm_2=hglm(X=matrix(rep(1,length(ids))[na_index]),y=phe_pre[ids,i][na_index],Z=Z[na_index,],calc.like = T,bigRR=T)
    p_mat_1[i,j]=tryCatch({summary(hglm_res)$FixCoefMat[,4]},error=function(e){return(1)})
    p_mat_2[i,j]=tryCatch({lrt(hglm_2,hglm_res)$p.value},error=function(e){return(1)})
    #p_mat_1[i,j]=1
    #p_mat_2[i,j]=1
    
    a=lm(phe_pre[ids,i][na_index]~root_50[ids,j][na_index]+pca[ids,1][na_index]+pca[ids,2][na_index]+pca[ids,3][na_index]+pca[ids,4][na_index]+pca[ids,5][na_index])
    b=lm(phe_pre[ids,i][na_index]~pca[ids,1][na_index]+pca[ids,2][na_index]+pca[ids,3][na_index]+pca[ids,4][na_index]+pca[ids,5][na_index])
    p_mat_3[i,j]=lrtest(a,b)[2,5]
    p_mat_4[i,j]=data.frame(summary(a)[4])[2,4]
    
    a=lm(phe_pre[ids,i]~root_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3])
    p_mat_5[i,j]=data.frame(summary(a)[4])[2,4]
  }
}
data=p_mat_1
data=p_mat_2
data=p_mat_3
data=p_mat_5
sub=t(apply(data,1,function(x){return(ifelse(x<0.05,1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.05/length(x),1,0))}))
row.names(sub)=names(phe_pre)
pheatmap(sub,cluster_rows = F,cluster_cols = F)


# OSCA ORM matrix hglm ----------------------------------------------------
#osca --efile /data/public/hanyu/Meta_tobacco/Data/phe_filtered/leaf_all_nosick_NA.txt --gene-expression --make-bod --no-fid --out /data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/orm/leaf_bod
#osca --befile /data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/orm/leaf_bod --make-orm-gz --out /data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/orm/leaf_orm
#osca --efile /data/public/hanyu/Meta_tobacco/Data/phe_filtered/root_all_nosick_NA.txt --gene-expression --make-bod --no-fid --out /data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/orm/root_bod
#osca --befile /data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/orm/root_bod --make-orm-gz --out /data/public/hanyu/Meta_tobacco/Data/final_vcf_0418/orm/root_orm
data=data.frame(fread(file.choose()),stringsAsFactors = F)
ids=data.frame(fread(file.choose()))
root_orm=matrix(ncol=max(data[,1]),nrow=max(data[,1]))
for(i in 1:nrow(data)){
  a=as.numeric(data[i,1:4])
  root_orm[a[1],a[2]]=a[4] 
}
root_orm[upper.tri(root_orm)] <- t(root_orm)[upper.tri(root_orm)]
row.names(root_orm)=ids$V1

data=data.frame(fread(file.choose()),stringsAsFactors = F)
ids=data.frame(fread(file.choose()))
leaf_orm=matrix(ncol=max(data[,1]),nrow=max(data[,1]))
for(i in 1:nrow(data)){
  a=as.numeric(data[i,1:4])
  leaf_orm[a[1],a[2]]=a[4] 
}
leaf_orm[upper.tri(leaf_orm)] <- t(leaf_orm)[upper.tri(leaf_orm)]
row.names(leaf_orm)=ids$V1

ids=intersect(row.names(ibs),row.names(leaf_orm))
index1=row.names(ibs)%in%ids
index2=row.names(leaf_orm)%in%ids
data=leaf_orm[index2,index2]
data2=ibs[index1,index1]
index3=c()
for(i in row.names(data)){index3=c(index3,which(row.names(data2)==i))}
data2=data2[index3,index3]
pheatmap(data,cluster_cols = F,cluster_rows = F,labels_row = "")
pheatmap(data2,cluster_cols = F,cluster_rows = F,labels_row = "",labels_col = "")
plot(zscore(data[upper.tri(data)]),zscore(data2[upper.tri(data2)]),ylab="ibs",xlab="leaf asv orm")
cor.test(zscore(data[upper.tri(data)]),zscore(data2[upper.tri(data2)]))


ids=intersect(row.names(ibs),row.names(root_orm))
index1=row.names(ibs)%in%ids
index2=row.names(root_orm)%in%ids
data=root_orm[index2,index2]
data2=ibs[index1,index1]
index3=c()
for(i in row.names(data)){index3=c(index3,which(row.names(data2)==i))}
data2=data2[index3,index3]
pheatmap(data,cluster_cols = F,cluster_rows = F,labels_row = "")
pheatmap(data2,cluster_cols = F,cluster_rows = F,labels_row = "",labels_col = "")
plot(zscore(data[upper.tri(data)]),zscore(data2[upper.tri(data2)]),ylab="ibs",xlab="root asv orm")
cor.test(zscore(data[upper.tri(data)]),zscore(data2[upper.tri(data2)]))

### random effect
p_mat_1=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(root_orm),intersect(row.names(root_50),row.names(phe_pre)))
  for(j in 1:ncol(root_50)){
    index=row.names(root_orm)%in%ids
    na_index=!is.na(root_50[ids,j])
    Z=c.z.hglm(root_orm[index,index])
    hglm_res=hglm(X=matrix(as.numeric(root_50[ids,j][na_index])),y=as.numeric(phe_pre[ids,i][na_index]),Z=Z[na_index,],calc.like = T,bigRR=T)
    p_mat_1[i,j]=tryCatch({summary(hglm_res)$FixCoefMat[,4]},error=function(e){return(1)})
  }
}

p_mat_2=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(leaf_orm),intersect(row.names(leaf_50),row.names(phe_pre)))
  for(j in 1:ncol(leaf_50)){
    index=row.names(leaf_orm)%in%ids
    na_index=!is.na(leaf_50[ids,j])
    Z=c.z.hglm(leaf_orm[index,index])
    hglm_res=hglm(X=matrix(as.numeric(leaf_50[ids,j][na_index])),y=as.numeric(phe_pre[ids,i][na_index]),Z=Z[na_index,],calc.like = T,bigRR=T)
    p_mat_2[i,j]=tryCatch({summary(hglm_res)$FixCoefMat[,4]},error=function(e){return(1)})
  }
}




data=p_mat_1
data=p_mat_2
sub=t(apply(data,1,function(x){return(ifelse(x<0.05,1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.05/length(x),1,0))}))
row.names(sub)=names(phe_pre)
pheatmap(sub,cluster_rows = F,cluster_cols = F)

###pc 


p_mat_1=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
e_mat_1=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(root_orm),intersect(row.names(root_50),row.names(phe_pre)))
  for(j in 1:ncol(root_50)){
    index=row.names(root_orm)%in%ids
    na_index=!is.na(root_50[ids,j])
    Z=c.z.hglm(root_orm[index,index])
    a=lm(phe_pre[ids,i]~root_50[ids,j]+Z[,1]+Z[,2]+Z[,3]+Z[,4]+Z[,5])
    p_mat_1[i,j]=data.frame(summary(a)[4])[2,4]
    e_mat_1[i,j]=data.frame(summary(a)[4])[2,1]
    
  }
}
row.names(p_mat_1)=names(phe_pre)
names(p_mat_1)=names(root_50)
row.names(e_mat_1)=names(phe_pre)
names(e_mat_1)=names(root_50)
root_mwas=list("Beta"=e_mat_1,"P"=p_mat_1)


p_mat_2=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
e_mat_2=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(leaf_orm),intersect(row.names(leaf_50),row.names(phe_pre)))
  for(j in 1:ncol(leaf_50)){
    index=row.names(leaf_orm)%in%ids
    na_index=!is.na(leaf_50[ids,j])
    Z=c.z.hglm(leaf_orm[index,index])
    a=lm(phe_pre[ids,i]~leaf_50[ids,j]+Z[,1]+Z[,2]+Z[,3]+Z[,4]+Z[,5])
    p_mat_2[i,j]=data.frame(summary(a)[4])[2,4]
    e_mat_2[i,j]=data.frame(summary(a)[4])[2,1]
  }
}
row.names(p_mat_2)=names(phe_pre)
names(p_mat_2)=names(leaf_50)
row.names(e_mat_2)=names(phe_pre)
names(e_mat_2)=names(leaf_50)
leaf_mwas=list("Beta"=e_mat_2,"P"=p_mat_2)
save(leaf_mwas,root_mwas,file="Results/MWAS/24.06.06_MWAS.RData")

data=p_mat_1
data=p_mat_2
sub=t(apply(data,1,function(x){return(ifelse(x<0.001,1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.05,1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.05/length(x),1,0))}))
row.names(sub)=names(phe_pre)
pheatmap(sub,cluster_rows = F,cluster_cols = F)


plot(phe_pre[ids,i],leaf_50[ids,j])
plot(phe_pre[ids,i],leaf_50[ids,2])

cor.test(phe_pre[ids,i],leaf_50[ids,j])
cor.test(phe_pre[ids,i],leaf_50[ids,2])

png("test.png",width=1000,height=800)
par(mfrow=c(5,6))
for(i in 1:29){
  hist(phe_pre[row.names(root_50),1],main=names(phe_pre)[i])  
}
dev.off()
