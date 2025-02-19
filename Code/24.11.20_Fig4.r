# 4A plot GWAS ------------------------------------------------------------
library(data.table)
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

pdf("Writing//V1_24.09.29//Fig4//4A.pdf",width=12,height=8)
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
  
  if(i=="PH"){
    plot(x=data[,3],y=-log10(data[,9]),col=all_col,cex=all_cex,
         frame.plot=F,pch=all_pch,ylab="-Log10(p)",xlab="Chrosome",xaxt="n",ylim=c(3,9))  
  }else{
    points(x=data[,3],y=-log10(data[,9]),col=all_col,pch=all_pch,cex=all_cex,) 
  }
}
chr_mean=aggregate(bp~Chr,data=data,mean)[,2]
axis(side=1,at=chr_mean,1:24)
abline(h=-log10(5.76E-7),col="red",lty="dashed")
legend("topright",names(cols),pch=pchs,col=cols,bty="n",cex=1.5)
dev.off()

# 4B venn -----------------------------------------------------------------
faidx=data.frame(fread("Data//ZY300.fa.fai"))
faidx=faidx[faidx[,1]%in%paste0("Chr",1:24),]
faidx[,1]=as.numeric(sub("Chr","",faidx[,1]))
row.names(faidx)=faidx[,1]
faidx=faidx[order(faidx[,1]),]
bin=1000000
pwd="WAS_summary/Phe_summary/top_sug_snp"
files=paste0(phe_info$name[c(-22,-16,-18)],"_sig.txt")
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
phe_data=apply(phe_gwa_sum,2,sum)
library(ggvenn)
pdf("Writing//V1_24.09.29//Fig4//4B.pdf",width=6,height=6)
ggvenn(list(Phyllosphere=names(leaf_data[leaf_data>=1]),Rhizosphere=names(root_data[root_data>=1]),"Agronomic traits"=names(phe_data[phe_data>=1])))  
dev.off()


# 4C coloc pp bar ---------------------------------------------------------
set1=intersect(intersect(names(leaf_data[leaf_data>=1]),names(root_data[root_data>=1])),names(phe_data[phe_data>=1]))
set2=setdiff(intersect(names(leaf_data[leaf_data>=1]),names(phe_data[phe_data>=1])),set1)
set3=setdiff(intersect(names(root_data[root_data>=1]),names(phe_data[phe_data>=1])),set1)

bin_inter=data.frame(chr=0,ps_chr=0,ps_genom=0,
                     type="",trait="",asv_inleaf="",leaf_feature_asv="",
                     asv_inroot="",root_feature_asv="",genes="")

gff=data.frame(fread("Data//anno//ZY300.gff3"))
gff=gff[gff$V3=="gene",]
gff$gene=gsub("(.*)\\.1","\\1",sub("ID=","",str_split_fixed(gff$V9,";",2)[,1]))

all_pos=rep(c("Phe_Leaf_Root","Phe_Leaf","Phe_Root"),c(9,5,23))
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
bin_inter[,2]=as.numeric(bin_inter[,2])
library(coloc)
asv1=unique(unlist(str_split(bin_inter$asv_inleaf,"_")))
asv2=unique(unlist(str_split(bin_inter$asv_inroot,"_")))
traits=unique(unlist(str_split(bin_inter$trait,"_")))
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


a=hist(unlist(h4_all),breaks=c(0,0.5,0.7,0.9,1),freq=T,main="H4 of all pairs")
table(unlist(h4_all)>0.8)
table(bin_inter$type)
h4_all[[9]]

# d1=hist(unlist(h4_all[1:9]),breaks=c(0,0.5,0.7,0.9,1))$counts
# d2=hist(unlist(h4_all[10:14]),breaks=c(0,0.5,0.7,0.9,1))$counts
# d3=hist(unlist(h4_all[15:37]),breaks=c(0,0.5,0.7,0.9,1))$counts

d1=unlist(lapply(h4_all[1:9],function(x){return(x[1])}))
d1=c(d1,unlist(h4_all[10:14]))
d2=unlist(lapply(h4_all[1:9],function(x){return(x[-1])}))
d2=c(d2,unlist(h4_all[15:37]))

d1=hist(d1,breaks=c(0,0.5,0.7,0.9,1))$counts
d2=hist(d2,breaks=c(0,0.5,0.7,0.9,1))$counts

data=rbind(d2,d1)
pdf("Writing//V1_24.09.29//Fig4//4C.pdf",width=5,height=6)
barplot(t(data),
        names.arg = c("Trait-Rhizosphere","Trait-Phyllosphere"),
        col=c("#C5E8C7","#CDE1A6","#F8DA94","#EF8376"),ylim=c(0,70),main="",ylab="Pair number")
text(2.3,50,"Coloc PP4")
legend(2,48,c("0-0.5","0.5-0.7","0.7-0.9","0.9-1"),fill=c("#C5E8C7","#CDE1A6","#F8DA94","#EF8376"),bty="n")
dev.off()

# DEF locus zoom ----------------------------------------------------------
i=9
j=3
asv2=unlist(str_split(bin_inter$asv_inroot[i],"_"))
traits=unique(unlist(str_split(bin_inter$trait[i],"_")))
root_combine=expand.grid(asv2,traits)
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

abf_res$results$snp[abf_res$results$SNP.PP.H4>0.2]
snps=abf_res$results$snp
"19:84250541"
micro1=data.frame(fread(paste0("Data/coloc_mlms/root/","ASV141",".mlma")))
row.names(micro1)=micro1$SNP
micro1=micro1[snps,]
# micro2=data.frame(fread(paste0("Data/coloc_mlms/root/","ASV40",".mlma")))
# row.names(micro2)=micro2$SNP
# micro2=micro2[snps,]
phe1=data.frame(fread(paste0("Data/coloc_mlms/phe/","RS",".mlma")))
row.names(phe1)=phe1$SNP
phe1=phe1[snps,]

ld=data.frame(fread("Data/genetic/84250541ld.txt"))
row.names(ld)=ld$SNP_B
ldr2=ld[snps,"R2"]
cols=rep("white",length(ldr2))
cols[ldr2<0.2]="skyblue"
cols[ldr2>=0.2 & ldr2<0.4]="darkblue"
cols[ldr2>=0.4 & ldr2<0.6]="orange"
cols[ldr2>=0.6 & ldr2<0.8]="tomato"
cols[ldr2>=0.8]="red"
# library(RColorBrewer)
# color_assign <- colorRamp2(breaks = c(0,0.4,0.7), 
#                            col =c("skyblue","orange","red"))
# cols=color_assign(ldr2)
pdf("Writing//V1_24.09.29//Fig4//4D.pdf",width=7,height=7)
plot(phe1$bp,-log10(phe1$p),col=cols,pch=19,frame.plot = F,ylab="-log10(p)",
     xlab="Chrosome 19 position (Mb)",ylim=c(0,7),xaxt="n")
axis(side=1,at=seq(83127241,85131686,400000),labels=seq(83.1,85.1,0.4))
points(phe1[c("19:84250541"),"bp"],
       -log10(phe1[c("19:84250541"),"p"]),
       col=c("purple"),pch=18,cex=2)
text(x=84450541,y=6.56,"19:84250541")
legend("topright",c("R2<0.2","0.2<=R2<0.4","0.4<=R2<0.6","0.6<=R2<0.8","R2>=0.8"),
       fill=c("skyblue","darkblue","orange","tomato","red"),bty="n")
dev.off()
pdf("Writing//V1_24.09.29//Fig4//4E.pdf",width=7,height=7)
plot(micro1$bp,-log10(micro1$p),col=cols,pch=19,frame.plot = F,ylab="-log10(p)",
     xlab="Chrosome 19 position (Mb)",ylim=c(0,7),xaxt="n")
axis(side=1,at=seq(83127241,85131686,400000),labels=seq(83.1,85.1,0.4))
points(micro1[c("19:84250541"),"bp"],
       -log10(micro1[c("19:84250541"),"p"]),
       col=c("purple"),pch=18,cex=2)
text(x=84450541,y=6.4,"19:84250541")
dev.off()