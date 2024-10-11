#install.packages("vegan")
#install.packages("devtools")
#install.packages("profvis")
#devtools::install_github('Hy4m/linkET',force = TRUE)
#devtools::install_local("C:/Users/MSI/Downloads/ggcor-1-master.zip")
library(ggcor)
library(linkET)
library(vegan)
library(data.table)
# get relationship --------------------------------------------------------
phe_pre=data.frame(fread("Data//phe_filtered//Phe_all_predicted.txt"))
phe_pre=phe_pre[!is.na(phe_pre[,1]),]
row.names(phe_pre)=phe_pre[,1]
all_id=intersect(root_phe[,1],leaf_phe[,1])
all_id=intersect(all_id,phe_pre[,1])
all_id=as.character(all_id)
phe_man=phe_pre[as.character(all_id),3:ncol(phe_pre)]
asv_all=cbind(root_phe[as.character(all_id),3:ncol(root_phe)],leaf_phe[as.character(all_id),3:ncol(leaf_phe)])
length(3:ncol(root_phe))
length(3:ncol(leaf_phe))
pheatmap(asv_all,cluster_rows = F,cluster_cols = F)
asv_all[is.na(asv_all)]=0
mantel = mantel_test(
  spec = asv_all, env = phe_man,
  spec_select = list(Root =1:685,Leaf= 686:781), 
  spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), # 样本距离使用的vegdist()计算，可以选择适合自己数据的距离指数。
  env_dist = dist_func(.FUN = "vegdist", method = "euclidean"),
  mantel_fun = 'mantel')
mantel=mantel[mantel$env!="Nic",]
mantel$P_adj_BH <- p.adjust(mantel$p, method = 'BH')
mantel = mutate(mantel, 
                r = cut(r, right = TRUE,# 表示分割区间形式为(a1,a2],前开后闭。
                        breaks = c(-Inf,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,Inf),
                        labels = c('(-1,-0.2]','(-0.2,-0.15]','(-0.15,-0.10]','(-0.10,-0.05]','(-0.05,0]','(0, 0.05]', '(0.05, 0.10]', '(0.10, 0.15]','(0.15, 0.2]', '(0.2, 1]')),
                P_adj_BH = cut(P_adj_BH, right = FALSE,# 表示分割区间形式为[a1,a2)。
                        breaks = c(-Inf, 0.001, 0.01, 0.05,0.1,0.2, Inf), 
                        labels = c('<0.001', '0.001-0.01', '0.01-0.05','0.05-0.1','0.1-0.2', '>0.2'),),
                p = cut(p, right = FALSE,# 表示分割区间形式为[a1,a2)。
                       breaks = c(-Inf, 0.001, 0.01, 0.05,0.1,0.2, Inf), 
                       labels = c('<0.001', '0.001-0.01', '0.01-0.05','0.05-0.1','0.1-0.2', '>0.2'),))
            
mantel = data.frame(mantel,stringsAsFactors = FALSE)

length(c('(-1,-0.2]','(-0.2,-0.15]','(-0.15,-0.10]','(-0.10,-0.05]','(-0.05,0]','(0, 0.05]', '(0.05, 0.10]', '(0.10, 0.15]','(0.15, 0.2]', '(0.2, 1]'))
length(c('<0.001', '0.001-0.01', '0.01-0.05','0.05-0.1','0.1-0.2', '>0.2'))
library(RColorBrewer)
display.brewer.pal(n = 12, name = "Set3")

p = 
  #先计算与绘制土壤因子之间相关性热图
  correlate(phe_man,method = "pearson") %>% #计算相关性
  qcorrplot(type = "lower",#表示位于下三角
            diag = F, #去除自相关
            grid_col = "black", # 网格线颜色
            grid_size = 0.25)+ # 网格大小
  geom_square() + # 以颜色和矩形面积表示相关性大小
  # 添加r值与显著性标记
  geom_mark(
    sep = '\n', 
    size = 2,
    sig_level = c(0.05, 0.01, 0.001), # 显著性水平设置
    sig_thres = 0.05 # 显著性阈值，p值大于0.05的相关性系数不绘制。
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(colors =c("#3FA9F5", "white", "#e74a32"),space="Lab")(10),
    limits = c(-1, 1),
    breaks = seq(-1,1,0.5))+ # 颜色区间
  #绘制连线
  geom_couple(aes(color = p, size = r), data = mantel, 
              label.size = 4,
              drop = TRUE,
              label.colour = "black",
              label.fontface = 1,
              nudge_x = 1,
              curvature = 0.2)+#线的弯曲程度
  scale_size_manual(values = c(0.2, 0.5,0.8,1.2)) +  #根据设置的mantel相关性r值区间设置线粗细
  scale_color_manual(values = c("tomato","#BC80BD","#BC80BD","pink","grey")) +  #根据设置的mantel相关性p值区间设置线颜色,因为p结果只有两个区间，可以只设置两个颜色。
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(legend.key = element_blank())
p

p_fdr = 
  #先计算与绘制土壤因子之间相关性热图
  correlate(phe_man,method = "pearson") %>% #计算相关性
  qcorrplot(type = "lower",#表示位于下三角
            diag = F, #去除自相关
            grid_col = "black", # 网格线颜色
            grid_size = 0.25)+ # 网格大小
  geom_square() + # 以颜色和矩形面积表示相关性大小
  # 添加r值与显著性标记
  geom_mark(
    sep = '\n', 
    size = 2,
    sig_level = c(0.05, 0.01, 0.001), # 显著性水平设置
    sig_thres = 0.05 # 显著性阈值，p值大于0.05的相关性系数不绘制。
  ) +
  scale_fill_gradientn(
    colours = colorRampPalette(colors =c("#3FA9F5", "white", "#e74a32"),space="Lab")(10),
    limits = c(-1, 1),
    breaks = seq(-1,1,0.5))+ # 颜色区间
  #绘制连线
  geom_couple(aes(color = P_adj_BH, size = r), data = mantel, 
              label.size = 4,
              drop = TRUE,
              label.colour = "black",
              label.fontface = 1,
              nudge_x = 1,
              curvature = 0.2)+#线的弯曲程度
  scale_size_manual(values = c(0.2, 0.5,0.8,1.2)) +  #根据设置的mantel相关性r值区间设置线粗细
  scale_color_manual(values = c("tomato","#BC80BD","#BC80BD","pink","grey")) +  #根据设置的mantel相关性p值区间设置线颜色,因为p结果只有两个区间，可以只设置两个颜色。
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(legend.key = element_blank())
p_fdr


table(mantel$p,mantel$spec)
table(mantel$P_adj_BH,mantel$spec)
table(mantel$p,mantel$spec)

save(mantel,file="Results//mantel//24.06.14_mantel.RData")



row.names(phe_info)=phe_info[,1]
data=cbind(mantel,phe_info[mantel$env,"type"],phe_info[mantel$env,"Description"])
names(data)[6:7]=c("type","des")
a=table(data$p,data$spec,data$type)
a=table(data$P_adj_BH,data$spec,data$type)
a


# MWAS analysis -----------------------------------------------------------
library(pheatmap)
leaf_50=data.frame(fread("Data//phe_filtered//leaf_50_nosick_NA.txt"))
root_50=data.frame(fread("Data//phe_filtered//root_50_nosick_NA.txt"))
row.names(root_50)=root_50[,1]
row.names(leaf_50)=leaf_50[,1]
root_50=root_50[,c(-1,-2)]
leaf_50=leaf_50[,c(-1,-2)]

phe_pre=data.frame(fread("Data//phe_filtered//Phe_all_raw.txt"))
phe_pre=phe_pre[!is.na(phe_pre[,1]),]
row.names(phe_pre)=phe_pre[,1]
phe_pre=phe_pre[,c(-1,-2)]

pca=data.frame(fread("Data//genetic//gcta_pca.eigenvec"))
elv=data.frame(fread("Data//genetic//gcta_pca.eigenval"))
row.names(pca)=pca[,1]
pca=pca[,3:7]

p_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
est_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(root_50))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(root_50),row.names(phe_pre))
  for(j in 1:ncol(root_50)){
    if(sum(is.na(phe_pre[ids,i]))==length(ids)){break}
    a=lm(phe_pre[ids,i]~root_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5])
    a=data.frame(summary(a)[4])
    est_mat[i,j]=a[2,1]
    p_mat[i,j]=a[2,4]
  }
}
row.names(p_mat)=names(phe_pre)
names(p_mat)=names(root_50)
row.names(est_mat)=names(phe_pre)
names(est_mat)=names(root_50)
root_mwas=list("Beta"=est_mat,"P"=p_mat)

p_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
est_mat=matrix(nrow=ncol(phe_pre),ncol=ncol(leaf_50))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(leaf_50),row.names(phe_pre))
  for(j in 1:ncol(leaf_50)){
    if(sum(is.na(phe_pre[ids,i]))==length(ids)){break}
    a=lm(phe_pre[ids,i]~leaf_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5])
    a=data.frame(summary(a)[4])
    est_mat[i,j]=a[2,1]
    p_mat[i,j]=a[2,4]
  }
}
row.names(p_mat)=names(phe_pre)
names(p_mat)=names(leaf_50)
row.names(est_mat)=names(phe_pre)
names(est_mat)=names(leaf_50)
leaf_mwas=list("Beta"=est_mat,"P"=p_mat)

data=root_mwas$P
data=leaf_mwas$P
sub=t(apply(data,1,function(x){return(ifelse(x<0.05/41,1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.05,1,0))}))
sub=t(apply(data,1,function(x){return(ifelse(x<0.01,1,0))}))
row.names(sub)=names(phe_pre)
pheatmap(sub,cluster_rows = F,cluster_cols = F)
pheatmap(root_mwas$Beta,cluster_rows = F,cluster_cols = F)


# micro net analysis ------------------------------------------------------
library(WGCNA)
power1=6
cutheight=0.5
deepsplit=2
type = "unsigned"
corType = "pearson"
robustY=T
maxPOutliers=1

dataExpr <- leaf_50
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.80,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

power=7
#run time 1h
dataExpr[is.na(dataExpr)]=0
net_leaf=blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = type, minModuleSize = 5,
                         reassignThreshold = 0, mergeCutHeight = 0.9,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=T,  saveTOMFileBase="Results//WGCNA//leaf_power7",
                         corType = corType, deepSplit = 3,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         verbose = 3)
table(net_leaf$colors)
moduleLabels = net_leaf$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net_leaf$dendrograms[[1]], moduleColors[net_leaf$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



dataExpr <- root_50
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.80,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

power=5
#run time 1h
dataExpr[is.na(dataExpr)]=0
net_root=blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                          TOMType = type, minModuleSize = 5,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs=T,  saveTOMFileBase="Results//WGCNA//root_power5",
                          corType = corType, deepSplit = 4,
                          maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                          verbose = 3)
table(net_root$colors)
moduleLabels = net_root$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net_root$dendrograms[[1]], moduleColors[net_root$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

data=matrix(ncol=ncol(net_root$MEs),nrow=ncol(phe_pre))
for(i in 1:ncol(net_root$MEs)){
  for(j in 1:ncol(phe_pre)){
    ids=as.character(intersect(row.names(net_root$MEs),row.names(phe_pre)))
    data[j,i]=cor.test(phe_pre[ids,j],net_root$MEs[ids,i])$p.value
    
  }
}
a=data
a=apply(data,1,function(x){ifelse(x<0.05/20,1,0)})
pheatmap(a)

table(root_feature[names(net_root$colors)[net_root$colors==2],"o"])


# micro PVE of PHE --------------------------------------------------------
data=data.frame(fread("Data//genetic//root_50_orm.orm.gz"),stringsAsFactors = F)
ids=data.frame(fread("Data//genetic//root_50_orm.orm.id"))
root_orm=matrix(ncol=max(data[,1]),nrow=max(data[,1]))
for(i in 1:nrow(data)){
  a=as.numeric(data[i,1:4])
  root_orm[a[1],a[2]]=a[4] 
}
root_orm[upper.tri(root_orm)] <- t(root_orm)[upper.tri(root_orm)]
row.names(root_orm)=ids$V1
root_orm=as.data.frame(root_orm)
names(root_orm)=ids$V1
root_orm[1:10,1:10]

data=data.frame(fread("Data//genetic//leaf_50_orm.orm.gz"),stringsAsFactors = F)
ids=data.frame(fread("Data//genetic//leaf_50_orm.orm.id"))
leaf_orm=matrix(ncol=max(data[,1]),nrow=max(data[,1]))
for(i in 1:nrow(data)){
  a=as.numeric(data[i,1:4])
  leaf_orm[a[1],a[2]]=a[4] 
}
leaf_orm[upper.tri(leaf_orm)] <- t(leaf_orm)[upper.tri(leaf_orm)]
row.names(leaf_orm)=ids$V1
leaf_orm=as.data.frame(leaf_orm)
names(leaf_orm)=ids$V1
leaf_orm[1:10,1:10]

library(hglm)
pves_data=data.frame(phe=names(phe_pre),root_pve=0,leaf_pve=0)
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(phe_pre),intersect(row.names(leaf_orm),row.names(root_orm)))
  Z1=c.z.hglm(root_orm[ids,ids])
  Z2=c.z.hglm(leaf_orm[ids,ids])
  model=hglm(X=matrix(rep(1,length(ids))),y=phe_pre[ids,i],Z=cbind(Z1,Z2),RandC=c(length(ids),length(ids)))
  pves_data[i,2:3]=model$varRanef/(sum(model$varRanef)+model$varFix)
}
pves_data$sum=pves_data$root_pve+pves_data$leaf_pve
data=t(pves_data[order(pves_data$sum,decreasing=T),1:3])
barplot(data[2:3,],col=c("tomato","skyblue"),border=NA,names.arg = data[1,])

###load phe h2
files=list.files("Data//phe_h2//")
files=files[grep("hsq",files)]
phe_info$h2=0
phe_info$h2_se=0
for(i in 1:ncol(phe_pre)){
  data=data.frame(fread(paste0("Data//phe_h2//",names(phe_pre)[i],".hsq")))
  phe_info[names(phe_pre)[i],"h2"]=data[4,2]
  phe_info[names(phe_pre)[i],"h2_se"]=data[4,3]
}

data=t(pves_data[order(pves_data$sum,decreasing=T),1:3])
barplot(phe_info$h2,col=rgb(0,1,0,0.3),border=NA,names.arg = data[1,],ylim=c(0,1),ylab="PVE",las=2)
par(new = TRUE)
barplot(data[2:3,],col=c("tomato","skyblue"),border=NA,names.arg = data[1,],ylim=c(0,1),yaxt="n",xaxt="n")
legend("topright",c("kinship h2","leaf PVE","root PVE"),fill=c(rgb(0,1,0,0.3),"skyblue","tomato"))


pheatmap(root_orm,cluster_rows = F,cluster_cols = F,main="root")
pheatmap(leaf_orm,cluster_rows = F,cluster_cols = F,main="leaf")

pves_fix_data=data.frame(trait=phe_info$Traits,leaf_fix=0,root_fix=0,leaf_rand=0,root_rand=0,sum_var=0,var_phe=0)
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(phe_pre),intersect(row.names(leaf_orm),row.names(root_orm)))
  index_leaf=leaf_mwas$P[names(phe_pre)[i],]<0.01
  index_root=root_mwas$P[names(phe_pre)[i],]<0.01
  X=cbind(rep(1,length(ids)),leaf_50[ids,index_leaf],root_50[ids,index_root])
  X[is.na(X)]=0
  for(j in 2:ncol(X)){
    X[,j]=zscore(X[,j])
  }
  X=as.matrix(X)
  Z1=c.z.hglm(root_orm[ids,ids])
  Z2=c.z.hglm(leaf_orm[ids,ids])
  model=hglm(X=X,y=phe_pre[ids,i],Z=cbind(Z1,Z2),RandC=c(length(ids),length(ids)))
  
  if(sum(index_leaf)==1){
    leaf_xb=var(X[,2:(sum(index_leaf)+1)]*model$fixef[2:(sum(index_leaf)+1)]  )
  }else if(sum(index_leaf)==0){
    leaf_xb=0
  }else{
    leaf_xb=var(X[,2:(sum(index_leaf)+1)]%*%model$fixef[2:(sum(index_leaf)+1)]  )
  }
  if(sum(index_root)==1){
    root_xb=var(X[,2:(sum(index_leaf)+2):(sum(index_root)+1)]*model$fixef[2:(sum(index_leaf)+2):(sum(index_root)+1)]  )
  }else if(sum(index_root)==0){
    root_xb=0
  }else{
    root_xb=var(X[,2:(sum(index_leaf)+2):(sum(index_root)+1)]%*%model$fixef[2:(sum(index_leaf)+2):(sum(index_root)+1)])
  }
  
  root_zu=var(Z1%*%model$ranef[1:length(ids)])
  leaf_zu=var(Z2%*%model$ranef[(length(ids)+1):(length(ids)*2)])
  
  pves_fix_data[i,"leaf_fix"]=leaf_xb
  pves_fix_data[i,"leaf_rand"]=leaf_zu
  pves_fix_data[i,"root_fix"]=root_xb
  pves_fix_data[i,"root_rand"]=root_zu
  pves_fix_data[i,"sum_var"]=var(X%*%model$fixef)+var(cbind(Z1,Z2)%*%model$ranef)
  pves_fix_data[i,"var_phe"]=var(phe_pre[ids,i])
}


data=pves_fix_data[,c(2,3,4,5)]
for(i in 1:4){
  data[,i]=data[,i]/pves_fix_data$var_phe
}
data=data[,c(3,1,2,4)]
barplot(t(data),names.arg = pves_fix_data$trait,col=c("tomato","orange","#CCEBC5","skyblue"),
        xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("root rand","root fix","leaf fix","leaf rand"),fill=c("skyblue","#CCEBC5","orange","tomato"),bty="n")


pves_fix_data$leaf_fix/pves_fix_data$var_phe
pves_fix_data$leaf_fix/pves_fix_data$sum_var
data=pves_fix_data[,c(2,3,4,5)]
sum=c(apply(data,1,sum))
for(i in 1:4){
  #data[,i]=data[,i]/pves_fix_data$var_phe
  #data[,i]=data[,i]/pves_fix_data$sum_var
  data[,i]=data[,i]/sum

}
data=data[,c(3,1,2,4)]
barplot(t(data),names.arg = pves_fix_data$trait,col=c("tomato","orange","#CCEBC5","skyblue"),
        xlim=c(0,40),las=2,ylab="PVE")
legend("topright",c("root rand","root fix","leaf fix","leaf rand"),fill=c("skyblue","#CCEBC5","orange","tomato"),bty="n")

d1=apply(root_mwas$P,1,function(x){sum(x<0.01)})
d2=apply(leaf_mwas$P,1,function(x){sum(x<0.01)})
barplot(d1,col=rgb(0,0.5,1,0.2),xlim=c(0,40),ylab="Number of ASVs",names.arg = NA,ylim=c(0,40))
par(new = TRUE)
barplot(d2,col=rgb(1,0,0,0.4),xlim=c(0,40),ylab="",names.arg = NA,ylim=c(0,40))
legend("topright",c("Root MWAS ASVs","Leaf MWAS ASVs"),fill=c(rgb(0,0.5,1,0.2),rgb(1,0,0,0.4)),bty="n")
       
plot(d1,data[,3])
plot(d2,data[,2])
