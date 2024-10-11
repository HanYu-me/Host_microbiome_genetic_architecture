library(data.table)


# get manhattanlike plot --------------------------------------------------
for(i in phe_info$name){
  system(paste0("cp ","WAS_summary/phe_summary/plot/",i,"_manhattan.png WAS_summary/phe_plot_inresearch/"))
}

names(table(phe_info$type))
pchs=c(15,16,17,18,19)
names(pchs)=names(table(phe_info$type))

cols=c("#C5E8C7","#CDE1A6","#F8DA94","#EF8376","#C49A84")
names(cols)=names(table(phe_info$type))

faidx=data.frame(fread("Data//ZY300.fa.fai"))
faidx=faidx[faidx[,1]%in%paste0("Chr",1:24),]
faidx[,1]=as.numeric(sub("Chr","",faidx[,1]))
row.names(faidx)=faidx[,1]
faidx=faidx[order(faidx[,1]),]
faidx[,2]=faidx[,2]+60000000
for(i in setdiff(phe_info$name,c("TMV","KC"))){
  data=data.frame(fread(paste0("WAS_summary/phe_summary/sub_mlma/",i,"_sub.txt")))
  for(j in 2:24){
    data[data[,1]==j,3]=as.numeric(data[data[,1]==j,3])+as.numeric(sum(faidx[1:(j-1),2]))
  }
  all_pch=rep(20,nrow(data))
  all_pch[data$p<5.76e-7]=pchs[phe_info[i,"type"]]
  all_col=rep("grey",nrow(data))
  all_col[data$p<5.76e-7]=cols[phe_info[i,"type"]]
  all_cex=rep(1,nrow(data))
  all_cex[data$p<5.76e-7]=1.3
  
  if(i=="BT"){
    plot(x=data[,3],y=-log10(data[,9]),col=all_col,cex=all_cex,
         frame.plot=F,pch=all_pch,ylab="-Log10(p)",xlab="Chrosome",xaxt="n",ylim=c(3,10))  
  }else{
    points(x=data[,3],y=-log10(data[,9]),col=all_col,pch=all_pch,cex=all_cex,) 
  }
}
chr_mean=aggregate(bp~Chr,data=data,mean)[,2]
axis(side=1,at=chr_mean,1:24)
abline(h=-log10(5.76E-7),col="red",lty="dashed")
legend("topright",names(cols),pch=pchs,col=cols)


# get bin info ------------------------------------------------------------
faidx=data.frame(fread("Data//ZY300.fa.fai"))
faidx=faidx[faidx[,1]%in%paste0("Chr",1:24),]
faidx[,1]=as.numeric(sub("Chr","",faidx[,1]))
row.names(faidx)=faidx[,1]
faidx=faidx[order(faidx[,1]),]
bin=1000000
pwd="WAS_summary/Phe_summary/top_sug_snp"
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
phe_gwa_sum[is.na(phe_gwa_sum)]=0
phe_gwa_sum=phe_gwa_sum[ setdiff(phe_info$name,c("TMV","KC")),]
save(phe_gwa_sum,file="Results//24.07.08_phe_gwaresults.RData")


chr_start=c(0)
for(j in 2:24){
  chr_start=c(chr_start,sum(faidx[1:(j-1),2]))
}
chr_start=c(chr_start,sum(faidx[,2]))
chr_end=chr_start[2:25]
chr_start=chr_start[1:24]
chr_mean=(chr_start+chr_end)/2

phe_data=apply(phe_gwa_sum,2,sum)
plot(x=names(phe_data),y=phe_data,pch=20,cex=1.5,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(phe_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated traits")
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")


# three bin analysis ------------------------------------------------------
library(ggvenn)
ggvenn(list(leaf=names(leaf_data[leaf_data>=1]),root=names(root_data[root_data>=1]),phe=names(phe_data[phe_data>=1])))  

set1=intersect(intersect(names(leaf_data[leaf_data>=1]),names(root_data[root_data>=1])),names(phe_data[phe_data>=1]))
set2=setdiff(intersect(names(leaf_data[leaf_data>=1]),names(phe_data[phe_data>=1])),set1)
set3=setdiff(intersect(names(root_data[root_data>=1]),names(phe_data[phe_data>=1])),set1)

bin_inter=data.frame(chr=0,ps_chr=0,ps_genom=0,
                     type="",trait="",asv_inleaf="",leaf_feature_asv="",
                     asv_inroot="",root_feature_asv="",genes="")

gff=data.frame(fread("Data//anno//ZY300.gff3"))
gff=gff[gff$V3=="gene",]
gff$gene=gsub("(.*)\\.1","\\1",sub("ID=","",str_split_fixed(gff$V9,";",2)[,1]))

all_pos=rep(c("Phe_Leaf_Root","Phe_Leaf","Phe_Root"),c(14,5,29))
names(all_pos)=c(set1,set2,set3)
for(i in c(set1,set2,set3)){
  chr=findInterval(i,chr_start)
  if(chr!=1){
    ps_chr=as.numeric(i)-sum(faidx[1:(chr-1),2])  
  }else{
    ps_chr=i
  }
  ps_genom=as.numeric(i)
  type=all_pos[i]
  trait=row.names(phe_gwa_sum)[phe_gwa_sum[,i]>=1]
  if(length(trait)>1){
    trait=paste(trait,collapse = "_")
  }
  asv_inleaf=row.names(leaf_gwa_sum)[leaf_gwa_sum[,i]>=1]
  leaf_feature_asv=leaf_feature[asv_inleaf,"o"]
  if(length(asv_inleaf)>1){
    asv_inleaf=paste(asv_inleaf,collapse = "_")
    leaf_feature_asv=paste(leaf_feature_asv,collapse = "_")
  }else if(length(asv_inleaf)==0){
    asv_inleaf=""
    leaf_feature_asv=""
  }
  
  asv_inroot=row.names(root_gwa_sum)[root_gwa_sum[,i]>=1]
  root_feature_asv=root_feature[asv_inroot,"o"]
  if(length(asv_inroot)>1){
    asv_inroot=paste(asv_inroot,collapse = "_")
    root_feature_asv=paste(root_feature_asv,collapse = "_")
  }else if(length(asv_inroot)==0){
    asv_inroot=""
    root_feature_asv=""
  }
  
  genes=gff[gff$V1==paste0("Chr",chr) & gff$V4>(ps_chr-bin) & gff$V4<(ps_chr+bin),"gene"]
  if(length(genes)>1){
    genes=paste(genes,collapse = "_")  
  }
  bin_inter=rbind(bin_inter,c(chr,ps_chr,ps_genom,type,trait,asv_inleaf,leaf_feature_asv,asv_inroot,root_feature_asv,genes))
}
bin_inter=bin_inter[-1,]
write.table(bin_inter,file="Results/three_bin_interval/bin_info.txt",row.names = F,col.names = T,quote=F)



# get interval go ---------------------------------------------------------
library(AnnotationHub)
library(clusterProfiler)
library(stringr)
library(ggplot2)
GO_result <- data.frame(fread("Data//anno//function.stat.txt"))
names(GO_result)
row.names(GO_result)=gsub("(.*)\\.1","\\1",GO_result$GeneID)
source("Code//Def_function.r")

tob_orgdb <- ah[["AH114234"]]

all_bin_go=lapply(bin_inter$genes,function(x){
  GO_list=get_go_fun(str_split(x,"_")[[1]])
  if(length(GO_list)<=3){
    return("NO go result")
  }else{
    res_GO <- enrichGO(GO_list,OrgDb=tob_orgdb,keyType = 'GO',ont='ALL',
                       pvalueCutoff=0.05, qvalueCutoff = 0.05)
    return(res_GO) 
  }
})

save(all_bin_go,bin_inter,file="Results/three_bin_interval/bins_go_result.RData")
i=0
i=i+1;dotplot(all_bin_go[[i]],split="ONTOLOGY",showCategory=10,font.size=7)+facet_grid(ONTOLOGY~., scale="free")


# correlation of asv and phe in interval bin  -----------------------------
leaf_mat=leaf_mwas$P
colnames(leaf_mat)=names(leaf_50)
leaf_p=c()
for(i in c(1:14,15:19)){
  traits=str_split(bin_inter$trait[i],"_")[[1]]
  asvs=str_split(bin_inter$asv_inleaf[i],"_")[[1]]
  leaf_p=c(leaf_p,c(leaf_mat[traits,asvs]))
}

root_mat=root_mwas$P
colnames(root_mat)=names(root_50)
root_p=c()
for(i in c(1:14,20:48)){
  traits=str_split(bin_inter$trait[i],"_")[[1]]
  asvs=str_split(bin_inter$asv_inroot[i],"_")[[1]]
  root_p=c(root_p,c(root_mat[traits,asvs]))
}

par(mfrow=c(2,1))
hist(-log10(leaf_p),breaks=seq(0,4,0.1),main="Correlation(p) of leaf ASV and trait in same loci")
hist(-log10(root_p),breaks=seq(0,4,0.1),main="Correlation(p) of root ASV and trait in same loci")


# coloc of asv and phe in interval bin ------------------------------------
### get mlm files into one file
asv1=unique(unlist(str_split(bin_inter$asv_inleaf,"_")))
asv2=unique(unlist(str_split(bin_inter$asv_inroot,"_")))
traits=unique(unlist(str_split(bin_inter$trait,"_")))
write.table(paste0("cp /data/public/hanyu/Meta_tobacco/Result/GWAS_filter/leaf50_NA/",asv1,".mlma /data/public/hanyu/Meta_tobacco/temp/coloc_mlms/leaf/ &"),
            file="info.txt",append=F,row.names=F,col.names=F,quote=F)
write.table(paste0("cp /data/public/hanyu/Meta_tobacco/Result/GWAS_filter/root50_NA/",asv2,".mlma /data/public/hanyu/Meta_tobacco/temp/coloc_mlms/root/ &"),
            file="info.txt",append=T,row.names=F,col.names=F,quote=F)
write.table(paste0("cp /data/public/hanyu/Meta_tobacco/Result/pheGWAS/",traits,".mlma /data/public/hanyu/Meta_tobacco/temp/coloc_mlms/phe/ &"),
            file="info.txt",append=T,row.names=F,col.names=F,quote=F)
library(coloc)
bin_inter$ps_chr=as.numeric(bin_inter$ps_chr)
h4_all=list()
for(i in 1:nrow(bin_inter)){
  print(i)
  asv1=unlist(str_split(bin_inter$asv_inleaf[i],"_"))
  asv2=unlist(str_split(bin_inter$asv_inroot[i],"_"))
  traits=unique(unlist(str_split(bin_inter$trait[i],"_")))
  leaf_combine=expand.grid(asv1,traits)
  root_combine=expand.grid(asv2,traits)
  h4_vec=c()
  if(bin_inter$type[i]=="Phe_Leaf"){
    for(j in 1:nrow(leaf_combine)){
      micro=data.frame(fread(paste0("Data/coloc_mlms/leaf/",as.character(leaf_combine[j,1]),".mlma")))
      phe=data.frame(fread(paste0("Data/coloc_mlms/phe/",as.character(leaf_combine[j,2]),".mlma")))
      micro=micro[micro[,1]==bin_inter[i,1]&micro[,3]>bin_inter[i,2]-bin&micro[,3]<bin_inter[i,2]+bin,]
      phe=phe[phe[,1]==bin_inter[i,1]&phe[,3]>bin_inter[i,2]-bin&phe[,3]<bin_inter[i,2]+bin,]
      
      micro=list(beta=micro[,7],varbeta=micro[,8],N=sum(!is.na(leaf_50[,as.character(leaf_combine[j,1])])),
                 sdY=sd(leaf_50[,as.character(leaf_combine[j,1])],na.rm=T),
                 type="quant",MAF=micro[,6],snp=micro[,2],position=micro[,3])
      phe=list(beta=phe[,7],varbeta=phe[,8],N=sum(!is.na(phe_pre[,as.character(leaf_combine[j,2])])),
               sdY=sd(phe_pre[,as.character(leaf_combine[j,2])],na.rm=T),
               type=ifelse(phe_info[as.character(leaf_combine[j,2]),"class"]=="Continues","quant","cc"),
               MAF=phe[,6],snp=phe[,2],position=phe[,3])
      abf_res <- coloc.abf (dataset1=micro,dataset2=phe)
      h4_vec=c(h4_vec,abf_res$summary[6])
    }
  }else if(bin_inter$type[i]=="Phe_Root"){
    for(j in 1:nrow(root_combine)){
      micro=data.frame(fread(paste0("Data/coloc_mlms/root/",as.character(root_combine[j,1]),".mlma")))
      phe=data.frame(fread(paste0("Data/coloc_mlms/phe/",as.character(root_combine[j,2]),".mlma")))
      micro=micro[micro[,1]==bin_inter[i,1]&micro[,3]>bin_inter[i,2]-bin&micro[,3]<bin_inter[i,2]+bin,]
      phe=phe[phe[,1]==bin_inter[i,1]&phe[,3]>bin_inter[i,2]-bin&phe[,3]<bin_inter[i,2]+bin,]
      
      micro=list(beta=micro[,7],varbeta=micro[,8],N=sum(!is.na(root_50[,as.character(root_combine[j,1])])),
                 sdY=sd(root_50[,as.character(root_combine[j,1])],na.rm=T),
                 type="quant",MAF=micro[,6],snp=micro[,2],position=micro[,3])
      phe=list(beta=phe[,7],varbeta=phe[,8],N=sum(!is.na(phe_pre[,as.character(root_combine[j,2])])),
               sdY=sd(phe_pre[,as.character(root_combine[j,2])],na.rm=T),
               type=ifelse(phe_info[as.character(root_combine[j,2]),"class"]=="Continues","quant","cc"),
               MAF=phe[,6],snp=phe[,2],position=phe[,3])
      abf_res <- coloc.abf (dataset1=micro,dataset2=phe)
      h4_vec=c(h4_vec,abf_res$summary[6])
    }
  }else{
    for(j in 1:nrow(leaf_combine)){
      micro=data.frame(fread(paste0("Data/coloc_mlms/leaf/",as.character(leaf_combine[j,1]),".mlma")))
      phe=data.frame(fread(paste0("Data/coloc_mlms/phe/",as.character(leaf_combine[j,2]),".mlma")))
      micro=micro[micro[,1]==bin_inter[i,1]&micro[,3]>bin_inter[i,2]-bin&micro[,3]<bin_inter[i,2]+bin,]
      phe=phe[phe[,1]==bin_inter[i,1]&phe[,3]>bin_inter[i,2]-bin&phe[,3]<bin_inter[i,2]+bin,]
      
      micro=list(beta=micro[,7],varbeta=micro[,8],N=sum(!is.na(leaf_50[,as.character(leaf_combine[j,1])])),
                 sdY=sd(leaf_50[,as.character(leaf_combine[j,1])],na.rm=T),
                 type="quant",MAF=micro[,6],snp=micro[,2],position=micro[,3])
      phe=list(beta=phe[,7],varbeta=phe[,8],N=sum(!is.na(phe_pre[,as.character(leaf_combine[j,2])])),
               sdY=sd(phe_pre[,as.character(leaf_combine[j,2])],na.rm=T),
               type=ifelse(phe_info[as.character(leaf_combine[j,2]),"class"]=="Continues","quant","cc"),
               MAF=phe[,6],snp=phe[,2],position=phe[,3])
      abf_res <- coloc.abf (dataset1=micro,dataset2=phe)
      h4_vec=c(h4_vec,abf_res$summary[6])
    }
    for(j in 1:nrow(root_combine)){
      micro=data.frame(fread(paste0("Data/coloc_mlms/root/",as.character(root_combine[j,1]),".mlma")))
      phe=data.frame(fread(paste0("Data/coloc_mlms/phe/",as.character(root_combine[j,2]),".mlma")))
      micro=micro[micro[,1]==bin_inter[i,1]&micro[,3]>bin_inter[i,2]-bin&micro[,3]<bin_inter[i,2]+bin,]
      phe=phe[phe[,1]==bin_inter[i,1]&phe[,3]>bin_inter[i,2]-bin&phe[,3]<bin_inter[i,2]+bin,]
      
      micro=list(beta=micro[,7],varbeta=micro[,8],N=sum(!is.na(root_50[,as.character(root_combine[j,1])])),
                 sdY=sd(root_50[,as.character(root_combine[j,1])],na.rm=T),
                 type="quant",MAF=micro[,6],snp=micro[,2],position=micro[,3])
      phe=list(beta=phe[,7],varbeta=phe[,8],N=sum(!is.na(phe_pre[,as.character(root_combine[j,2])])),
               sdY=sd(phe_pre[,as.character(root_combine[j,2])],na.rm=T),
               type=ifelse(phe_info[as.character(root_combine[j,2]),"class"]=="Continues","quant","cc"),
               MAF=phe[,6],snp=phe[,2],position=phe[,3])
      abf_res <- coloc.abf (dataset1=micro,dataset2=phe)
      h4_vec=c(h4_vec,abf_res$summary[6])
    }
  }
  h4_all[[i]]=h4_vec
}
save.image()

a=hist(unlist(h4_all),breaks=c(0,0.5,0.7,0.9,1),freq=T,main="H4 of all pairs")
table(unlist(h4_all)>0.9)
table(bin_inter$type)

d1=hist(unlist(h4_all[1:14]),breaks=c(0,0.5,0.7,0.9,1))$counts
d2=hist(unlist(h4_all[15:19]),breaks=c(0,0.5,0.7,0.9,1))$counts
d3=hist(unlist(h4_all[20:48]),breaks=c(0,0.5,0.7,0.9,1))$counts
data_plot=data.frame(class=rep(c("Commons bin of leaf, root, and phenotype",
                                 "Commons bin of leaf and phenotype",
                                 "Commons bin of root and phenotype"),each=4),
                     PP4=rep(c("No [0,0.5]","Low (0.5,0.7]","Mediate (0.7,0.9]","High (0.9,1]"),3),
                     value=c(d1,d2,d3))
data_plot$class=factor(data_plot$class,levels=c("Commons bin of leaf, root, and phenotype",
                                                "Commons bin of leaf and phenotype",
                                                "Commons bin of root and phenotype"))
data_plot$PP4=factor(data_plot$PP4,levels=c("No [0,0.5]","Low (0.5,0.7]","Mediate (0.7,0.9]","High (0.9,1]"))
ggplot(data_plot, aes( x = class,weight=value,fill = PP4))+
  geom_bar( width = 0.8)+
  scale_fill_manual( values = c("#C5E8C7","#CDE1A6","#F8DA94","#EF8376"))+
  ylab("Count")+xlab("")+
  theme(
    axis.title=element_text(size=15,face="plain",color="black"),
    axis.text = element_text(size=15,face="plain",color="black"),
    legend.title=element_text(size=15,face="plain",color="black"),
    legend.position = "right",
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black", size = 0.4))+theme_bw()+
  theme(text=element_text(family="A"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust=0.5))

for(i in 1:48){
  if(sum(h4_all[[i]]>0.9)==1){
    print(i)
  }
}
h4_all[[14]]
h4_all[[41]]
bin_inter[14,]
bin_inter[41,]

root_mat["RS",c("ASV40","ASV141")]


GO_list=get_go_fun(str_split(bin_inter[41,"genes"],"_")[[1]])
res_GO <- enrichGO(GO_list,OrgDb=tob_orgdb,keyType = 'GO',ont='ALL',
                     pvalueCutoff=0.05, qvalueCutoff = 0.05)
dotplot(res_GO,split="ONTOLOGY",showCategory=15,font.size=7)+facet_grid(ONTOLOGY~., scale="free")

# GBLUP of meta -----------------------------------------------------------
data=data.frame(fread("Data//genetic//Final_0607_pruned_grm_gz.grm.gz"),stringsAsFactors = F)
ids=data.frame(fread("Data//genetic//Final_0607_pruned_grm_gz.grm.id"))
snp_grm=matrix(ncol=max(data[,1]),nrow=max(data[,1]))
for(i in 1:nrow(data)){
  a=as.numeric(data[i,1:4])
  snp_grm[a[1],a[2]]=a[4] 
}
snp_grm[upper.tri(snp_grm)] <- t(snp_grm)[upper.tri(snp_grm)]
row.names(snp_grm)=ids$V1
snp_grm=as.data.frame(snp_grm)
names(snp_grm)=ids$V1
snp_grm[1:10,1:10]
pheatmap::pheatmap(snp_grm,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(leaf_orm,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(root_orm,cluster_rows = F,cluster_cols = F)

library(parallel)

phe_five_fold_GS_fun=function(i){
  ids=intersect(row.names(phe_man),intersect(intersect(names(root_orm),names(leaf_orm)),names(snp_grm)))
  phe_all=phe_pre[ids,]
  G1=snp_grm[ids,ids]
  G2=leaf_orm[ids,ids]
  G3=root_orm[ids,ids]
  idnum=length(ids)
  sample_index=split(c(1:(idnum)),sample(rep(1:5,c(idnum/5,idnum/5,idnum/5,idnum/5,idnum/5))))
  phe=phe_all[,i]
  phe_save=phe
  five_list=list()
  five_list[["phe_raw"]]=phe_save
  for(j in 1:5){
    phe=phe_save
    phe[sample_index[[j]]]=NA
    index=!is.na(phe)
    Z1=c.z.hglm(G1)
    Z2=c.z.hglm(G2)
    Z3=c.z.hglm(G3)
    if(phe_info[names(phe_man)[i],"class"]=="Categorical"){
      model=hglm(X=matrix(rep(1,length(phe[index]))),
                 y=phe[index],
                 Z=cbind(Z1[index,],Z2[index,],Z3[index,]),
                 RandC=c(ncol(Z1[index,]),ncol(Z2[index,]),ncol(Z3[index,])),
                 family=binomial(link = "logit"))
      model2=hglm(X=matrix(rep(1,length(phe[index]))),
                  y=phe[index],
                  Z=Z1[index,],
                  family=binomial(link = "logit"))
    }else{
      model=hglm(X=matrix(rep(1,length(phe[index]))),
                 y=phe[index],
                 Z=cbind(Z1[index,],Z2[index,],Z3[index,]),
                 RandC=c(ncol(Z1[index,]),ncol(Z2[index,]),ncol(Z3[index,])))
      model2=hglm(X=matrix(rep(1,length(phe[index]))),
                  y=phe[index],
                  Z=Z1[index,])
    }
    
    pre_all=cbind(Z1,Z2,Z3)%*%model$ranef
    pre_all2=Z1%*%model2$ranef
  
    five_list[[as.character(j)]]=list(three_pre=pre_all,one_pre=pre_all2,index=sample_index[[j]])
  }
  write.table(c(i,names(phe_man)[i]),"info.txt",append=T,col.names=F,row.names=F)
  return(five_list)
}
clnum<-15
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(phe_five_fold_GS_fun)))
clusterEvalQ(cl,library(hglm))
clusterExport(cl,"c.z.hglm")
#phe_pre=phe_pre[,names(phe_man)]
clusterExport(cl,c("phe_man","phe_pre","snp_grm","leaf_orm","root_orm","phe_info"))
result_1_5fold=parLapply(cl,1:27, phe_five_fold_GS_fun)
stopCluster(cl)

names(result_1_5fold)=names(phe_man)
save(result_1_5fold,file="Results/GS/result_5_fold.RData")
save(result_1_5fold,file="Results/GS/result_5_fold_raw_phe.RData")

data_plot=data.frame(trait=rep(names(result_1_5fold),each=10),class=rep(rep(c("SNP only","SNP and micro"),each=5),27),value=0)
for(i in names(result_1_5fold)){
  sub=result_1_5fold[[i]]
  c1=c()
  c2=c()
  for(j in 2:6){
    #c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
    #c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
    c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
    c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
  }
  data_plot[data_plot$trait==i&data_plot$class=="SNP only","value"]=c2
  data_plot[data_plot$trait==i&data_plot$class=="SNP and micro","value"]=c1
  print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
}

library(ggplot2)
library(ggpubr)
data_plot$trait=factor(data_plot$trait,levels=names(phe_man))
data_plot$class=factor(data_plot$class,levels=c("SNP only","SNP and micro"))
ggplot(data_plot, aes(x = trait, y = value, fill = class)) +
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ ylab("Gene methylation")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
  stat_compare_means(method = "t.test")

# PVE for meta ------------------------------------------------------------
pves_rand_data=data.frame(trait="",SNP_joint=0,leaf_joint=0,root_joint=0,SNP_only=0,all_joint=0,var_y=0,var_rand1_joint=0,var_rand2_joint=0,var_rand3_joint=0,var_e_joint=0,var_rand_only=0,var_e_only=0)
for(i in 1:27){
  print(i)
  ids=intersect(row.names(phe_man),intersect(intersect(names(root_orm),names(leaf_orm)),names(snp_grm)))
  phe_all=phe_man[ids,]
  G1=snp_grm[ids,ids]
  G2=leaf_orm[ids,ids]
  G3=root_orm[ids,ids]
  idnum=length(ids)
  sample_index=split(c(1:(idnum)),sample(rep(1:5,c(idnum/5,idnum/5,idnum/5,idnum/5,idnum/5))))
  phe=phe_all[,i]
  

  index=!is.na(phe)
  Z1=c.z.hglm(G1)
  Z2=c.z.hglm(G2)
  Z3=c.z.hglm(G3)
  if(phe_info[names(phe_man)[i],"class"]=="Categorical"){
    model=hglm(X=matrix(rep(1,length(phe[index]))),
               y=phe[index],
               Z=cbind(Z1[index,],Z2[index,],Z3[index,]),
               RandC=c(ncol(Z1[index,]),ncol(Z2[index,]),ncol(Z3[index,])),
               family=binomial(link = "logit"))
    model2=hglm(X=matrix(rep(1,length(phe[index]))),
                y=phe[index],
                Z=Z1[index,],
                family=binomial(link = "logit"))
  }else{
    model=hglm(X=matrix(rep(1,length(phe[index]))),
               y=phe[index],
               Z=cbind(Z1[index,],Z2[index,],Z3[index,]),
               RandC=c(ncol(Z1[index,]),ncol(Z2[index,]),ncol(Z3[index,])))
    model2=hglm(X=matrix(rep(1,length(phe[index]))),
                y=phe[index],
                Z=Z1[index,])
  }
  pves_rand_data[i,]=c(names(phe_man)[i],
                       var(Z1%*%model$ranef[1:ncol(Z1)]),
                       var(Z2%*%model$ranef[(1+ncol(Z1)):(ncol(Z1)+ncol(Z2))]),
                       var(Z3%*%model$ranef[(1+ncol(Z1)+ncol(Z2)):(ncol(Z1)+ncol(Z2)+ncol(Z3))]),
                       var(Z1%*%model2$ranef),
                       var(cbind(Z1,Z2,Z3)%*%model$ranef),
                       var(phe,na.rm=T),
                       model$varRanef[1],
                       model$varRanef[2],
                       model$varRanef[3],
                       model$varFix,
                       model2$varRanef,
                       model2$varFix)
}
names(pves_rand_data)
data=pves_rand_data[,c(8,9,10)]
sum=apply(pves_rand_data[,c(8,9,10,11)],1,function(x)sum(as.numeric(x)))
for(i in 1:3){
  data[,i]=as.numeric(data[,i])/sum
}
par(mfrow=c(2,1))
barplot(t(data),names.arg = pves_rand_data$trait,col=c("tomato","orange","#CCEBC5"),
        xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("root joint","leaf joint","SNP joint"),fill=c("#CCEBC5","orange","tomato"),bty="n")

barplot(as.numeric(pves_rand_data$var_rand_only)/(as.numeric(pves_rand_data$var_e_only)+as.numeric(pves_rand_data$var_rand_only)),
        names.arg = pves_rand_data$trait,xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("SNP only"),fill=c("grey"),bty="n")


a=Z1%*%Z2
pheatmap::pheatmap(Z1,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(Z2,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(a,cluster_rows = F,cluster_cols = F)
# meta net ----------------------------------------------------------------
