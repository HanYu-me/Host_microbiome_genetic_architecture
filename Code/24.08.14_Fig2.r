# h2 distribution with hglm ---------------------------------------------------------
library(hglm)
pb=txtProgressBar(style=3)
root_h2=c()
for(i in 1:ncol(root_50)){
  setTxtProgressBar(pb, i/ncol(root_50))
  ids=intersect(names(snp_grm),row.names(root_50))
  phe=root_50[ids,i]
  G=snp_grm[ids,ids]
  Z=c.z.hglm(G)
  index=!is.na(phe)
  model=hglm(y=phe[index],X=matrix(rep(1,sum(index))),Z=Z[index,index])
  root_h2=c(root_h2,model$varRanef/(model$varRanef+model$varFix))
}

leaf_h2=c()
for(i in 1:ncol(leaf_50)){
  setTxtProgressBar(pb, i/ncol(leaf_50))
  ids=intersect(names(snp_grm),row.names(leaf_50))
  phe=leaf_50[ids,i]
  G=snp_grm[ids,ids]
  Z=c.z.hglm(G)
  index=!is.na(phe)
  model=hglm(y=phe[index],X=matrix(rep(1,sum(index))),Z=Z[index,index])
  leaf_h2=c(leaf_h2,model$varRanef/(model$varRanef+model$varFix))
}

summary(root_h2)
summary(leaf_h2)

names(root_h2)=names(root_50)
names(leaf_h2)=names(leaf_50)
pdf("Figs/SFig/S1.pdf",width=8,height=4)
par(mfrow=c(1,2))
hist(root_h2,breaks=seq(0,1,0.05),main="",xlab="Kinship heritability of rhizosphere microbiome",ylab="Frequency (counts)")
hist(leaf_h2,breaks=seq(0,1,0.05),main="",xlab="Kinship heritability of leaf microbiome",ylab="Frequency (counts)")
dev.off()

### correlation of h2 and finded loci
a=apply(root_gwa_sum,1,sum)
inter=intersect(names(a),names(root_h2))
plot(root_h2[inter],a[inter])

b=apply(leaf_gwa_sum,1,sum)
inter=intersect(names(b),names(leaf_h2))
plot(leaf_h2[inter],b[inter])


# mental ------------------------------------------------------------------

ids=intersect(row.names(root_50),row.names(snp_grm))
a=root_50[ids,]
a[is.na(a)]=0
root_dist=vegan::vegdist(a)
genetic_dist=1-snp_grm[ids,ids]
mantel_Root<-vegan::mantel(genetic_dist, root_dist, method="spearman",permutations=9999,strata=NULL)

ids=intersect(row.names(leaf_50),row.names(snp_grm))
a=leaf_50[ids,]
a[is.na(a)]=0
leaf_dist=vegan::vegdist(a)
genetic_dist=1-snp_grm[ids,ids]
mantel_Leaf<-vegan::mantel(genetic_dist, leaf_dist, method="spearman",permutations=9999,strata=NULL)

mantel_Root
mantel_Leaf




# 2AB Manhattan like plot---------------------------------------------------------------
table(apply(root_gwa_sum,1,sum)>0)
table(apply(leaf_gwa_sum,1,sum)>0)

pdf("Writing//V1_24.09.29//Fig2//2A.pdf",width=12,height=4.5)
leaf_data=apply(leaf_gwa_sum,2,sum)
plot(x=names(leaf_data),y=leaf_data,pch=20,cex=2,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(leaf_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated ASVs in phyllosphere",ylim=c(2,max(leaf_data)))
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")
dev.off()

pdf("Writing//V1_24.09.29//Fig2//2B.pdf",width=12,height=4.5)
root_data=apply(root_gwa_sum,2,sum)
plot(x=names(root_data),y=root_data,pch=20,cex=2,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(root_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated ASVs in rhizosphere",ylim=c(2,max(root_data)))
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")
dev.off()


# 2C venn plot ------------------------------------------------------------
library(ggvenn)
g=ggvenn(list(leaf=names(leaf_data[leaf_data>=1]),root=names(root_data[root_data>=1])))
ggsave(g,filename="Writing//V1_24.09.29//Fig2//2C.pdf",width=6,height=6)
ggvenn(list(leaf=names(leaf_data[leaf_data>=2]),root=names(root_data[root_data>=2])))
# 2D genetic time cor -----------------------------------------------------
plot(root_data,leaf_data)
cor.test(root_data,leaf_data)$p.value
a=summary(lm(leaf_data~root_data))
a$coefficients
library(ggplot2)
g=ggplot(data.frame(x =root_data , y = leaf_data), aes(x = x, y = y)) +
  geom_hex(binwidth = c(1.7, 1.1))+scale_fill_gradient2(low="lightskyblue1",high="royalblue3",transform="log1p")+ 
  coord_fixed()+ylim(c(0,8))+xlim(c(0,15))+theme_classic()+
  xlab("Number of associated ASVs in root")+ylab("Number of associated ASVs in leaf")+
  geom_smooth(method = "lm", se=T, color="cyan2", formula = y ~ x)
ggsave(g,filename="Writing//V1_24.09.29//Fig2//2D.pdf",width=6,height=6)

a=cbind(root_data,leaf_data)

over_degree=matrix(ncol=7,nrow=16)
for(i in 0:15){
  for(j in 0:6){
    over_degree[i,j]=sum(a[,1]>i&a[,2]>j)
  }
}
library(pheatmap)
# pheatmap and display the value of each cell
pheatmap(over_degree,cluster_rows = F,cluster_cols = F,display_numbers = T)

sum(a[,1]>1&a[,2]>1)
sum(a[,1]>10)


# 2E PVE analysis ---------------------------------------------------------
load("Results//PVE//24.06.17_common_snp_PVE.RData")

PVE_ratio=data.frame(c(pves_leaf[,2]/(pves_leaf[,2]+pves_leaf[,3]),pves_root[,2]/(pves_root[,2]+pves_root[,3])))
PVE_ratio$class=c(rep("Leaf",nrow(pves_leaf)),rep("Root",nrow(pves_root)))
names(PVE_ratio)[1]="value"
library(gghalves) 
g=ggplot(PVE_ratio,aes(x=class,y=value,fill=class))+
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
ggsave(g,filename="Writing//V1_24.09.29//Fig2//2E.pdf",width=6,height=3)

summary(pves_leaf[,2]/(pves_leaf[,2]+pves_leaf[,3]))
summary(pves_root[,2]/(pves_root[,2]+pves_root[,3]))


# 2F example of chr14 --------------------------------------------------------
library(AnnotationHub)
library(clusterProfiler)
library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
tob_orgdb <- ah[["AH114234"]]
genes=gff[gff$V1=="Chr14" & gff$V4>(2315800001-2302807202) & gff$V4<(2316800001-2302807202),"gene"]
GO_list=get_go_fun(genes)
res_GO <- enrichGO(GO_list,OrgDb=tob_orgdb,keyType = 'GO',ont='BP',
                   pvalueCutoff=0.05, qvalueCutoff = 0.05)
#display GO result
pdf("Writing//V1_24.09.29//Fig2//2F.pdf",width=7,height=8)
barplot(res_GO,showCategory = 30)+scale_fill_gradient(low="tomato",high="skyblue")+theme_classic()
dev.off()
2315800001-2302807202
2316800001-2302807202
# S fig S data -------------------------------------------
### microbial abundance table
ids=intersect(row.names(root_50),intersect(row.names(leaf_50),row.names(snp_grm)))
write.table(root_50[ids,],file="Writing//V1_24.09.29//STable//S1.txt",quote=F,sep="\t")
write.table(leaf_50[ids,],file="Writing//V1_24.09.29//STable//S2.txt",quote=F,sep="\t")

### get snp sumamry data
all_snp=data.frame(fread("data/ASV1.mlma"))
head(all_snp)
all_snp=all_snp[,c(1,2,3,4,5,6)]
row.names(all_snp)=all_snp$SNP

pwd="WAS_summary/GWAS_filter_summary/leaf50_NA/top_sug_snp"
files=list.files(pwd)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  row.names(data)=data$SNP
  all_snp[row.names(data),sub("_sig.txt","",files[i])]=data[row.names(data),"p"]
}
index=apply(all_snp,1,function(x) sum(!is.na(x))>6)
all_snp=all_snp[index,]
write.table(all_snp,file="Writing//V1_24.09.29//STable//S3.txt",quote=F,sep="\t")

all_snp=data.frame(fread("data/ASV1.mlma"))
head(all_snp)
all_snp=all_snp[,c(1,2,3,4,5,6)]
row.names(all_snp)=all_snp$SNP
all_snp$judge=0

pwd="WAS_summary/GWAS_filter_summary/root50_NA/top_sug_snp"
files=list.files(pwd)
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0(pwd,"/",files[i])))
  row.names(data)=data$SNP
  all_snp[row.names(data),sub("_sig.txt","",files[i])]=data[row.names(data),"p"]
  all_snp[row.names(data),"judge"]=1
}
all_snp=all_snp[all_snp$judge>0,]
all_snp=all_snp[,-7]
write.table(all_snp,file="Writing//V1_24.09.29//STable//S4.txt",quote=F,sep="\t")
