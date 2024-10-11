library(ggcor)
library(linkET)
library(vegan)
library(data.table)


# phe h2 and phe info------------------------------------------------------------------
files=list.files("Data//phe_h2//")
files=files[grep("hsq",files)]
phe_h2=c()
phe_num=c()
for(i in files){
  data=data.frame(fread(paste0("Data//phe_h2//",i),fill=T))
  phe_h2=c(phe_h2,data[4,2])
  phe_num=c(phe_num,data[10,2])
}
par(mfrow=c(2,1))
barplot(phe_h2,ylim=c(0,1),main="h2")
barplot(phe_num,names.arg = sub(".hsq","",files),ylim=c(0,214),las=2,main="id num")
abline(h=100,col="red",lty="dashed")
table(phe_num>106)

phe_info=data.frame(fread("Data/anno/Trait_in_research.csv"))
phe_info[phe_info$class=="Continues",3]

info=data.frame(name=sub(".hsq","",files),h2=phe_h2,num=phe_num)
names(phe_info)[3]="name"
info=merge(info,phe_info,by="name")
info=info[,c(1,5,6,7,2,3)]
phe_info=info
row.names(phe_info)=phe_info$name
table(phe_info$num>100,phe_info$type)
phe_info=phe_info[phe_info$num>100,]
table(phe_info$type,phe_info$class)
phe_info["LL",4]="Continues"
# get mantel relationship --------------------------------------------------------
phe_pre=data.frame(fread("Data//phe_filtered//Phe_all_raw.txt"))
phe_pre=phe_pre[!is.na(phe_pre[,1]),]
row.names(phe_pre)=phe_pre[,1]
root_phe=data.frame(fread("Data//phe_filtered/root_50_nosick_NA.txt"))
row.names(root_phe)=root_phe[,1]
leaf_phe=data.frame(fread("Data//phe_filtered/leaf_50_nosick_NA.txt"))
row.names(leaf_phe)=leaf_phe[,1]
all_id=intersect(root_phe[,1],leaf_phe[,1])
all_id=intersect(all_id,phe_pre[,1])
all_id=as.character(all_id)
phe_man=phe_pre[as.character(all_id),rev(phe_info$name[order(phe_info$type,phe_info$class)])]
table(phe_man[,phe_info$name[phe_info$class=="Categorical"][4]])

plot(phe_man$ID,phe_man$IFD)
cor.test(phe_man$ID,phe_man$IFD,method="pearson")
index=!is.na(phe_man$ID)&!is.na(phe_man$IFD)
cor.test(phe_man[,1],phe_man[,2],method="pearson")
cor(phe_man[,1],phe_man[,2],na.rm=T)
ltm::biserial.cor(phe_man$ID[index],as.factor(phe_man$IFD)[index])
cor(phe_man$ID,phe_man$IFD,use="pairwise.complete.obs")


asv_all=cbind(root_phe[as.character(all_id),3:ncol(root_phe)],leaf_phe[as.character(all_id),3:ncol(leaf_phe)])
length(3:ncol(root_phe))
length(3:ncol(leaf_phe))
asv_all[is.na(asv_all)]=0

### convert continues phe into zscore 
for(i in phe_info$name[phe_info$class=="Continues"]){
  phe_man[,i]=zscore(phe_man[,i])
}


mantel = mantel_test(
  spec = asv_all, env = phe_man,
  spec_select = list(Root =1:685,Leaf= 686:781), 
  spec_dist =  dist_func(.FUN = "vegdist", method = "bray"), # 样本距离使用的vegdist()计算，可以选择适合自己数据的距离指数。
  env_dist = dist_func(.FUN = "vegdist", method = "euclidean"),
  mantel_fun = 'mantel')
View(mantel)
mantel$P_adj_BH <- p.adjust(mantel$p, method = 'BH')
mantel[mantel$p<0.05,]
summary(mantel$r)
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

table(mantel$r)
table(mantel$p)


library(ggplot2)
library(aplot)
p = 
  #先计算与绘制土壤因子之间相关性热图
  correlate(phe_man,method = "pearson",use="pairwise.complete.obs") %>% #计算相关性
  qcorrplot(type = "lower",#表示位于下三角
            diag = F, #去除自相关
            grid_col = "black", # 网格线颜色
            grid_size = 0.25)+ # 网格大小
  geom_square()+# 以颜色和矩形面积表示相关性大小
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
  scale_size_manual(values = c(1,0.7,0.3,0.3,0.7,1)) +  #根据设置的mantel相关性r值区间设置线粗细
  scale_color_manual(values = c("red","tomato","#BC80BD","grey","grey","white")) +  #根据设置的mantel相关性p值区间设置线颜色,因为p结果只有两个区间，可以只设置两个颜色。
  guides(color = guide_legend(title = "Mantel's p", order = 1), #图例标题和排序
         size = guide_legend(title = "Mantel's r", order = 2), 
         fill = guide_colorbar(title = "Pearson's r", order = 3)) +
  theme(legend.key = element_blank(),legend.position = "right")
p

group=names(phe_man) %>% as.data.frame() %>%
  mutate(group=phe_info[names(phe_man),"type"]) %>%
  mutate(p="")%>%
  ggplot(aes(p,.,fill=group))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values = c("#C5E8C7","#CDE1A6","#F8DA94","#EF8376","#C49A84"))+
  theme(axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5))+
  labs(fill = "Group")

p%>%insert_left(group, width = 0.01)




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
phe_pre=phe_pre[as.character(all_id),rev(phe_info$name[order(phe_info$type,phe_info$class)])]

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
    if(phe_info[names(phe_pre)[i],"class"]=="Categorical"){
      a=glm(as.factor(phe_pre[ids,i])~root_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5],family=binomial())
      est_mat[i,j]=summary(a)$coefficients[2,1]
      p_mat[i,j]=summary(a)$coefficients[2,4]
    }else{
      a=lm(phe_pre[ids,i]~root_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5])
      a=data.frame(summary(a)[4])
      est_mat[i,j]=a[2,1]
      p_mat[i,j]=a[2,4] 
    }
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
    if(phe_info[names(phe_pre)[i],"class"]=="Categorical"){
      a=glm(as.factor(phe_pre[ids,i])~leaf_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5],family=binomial())
      est_mat[i,j]=summary(a)$coefficients[2,1]
      p_mat[i,j]=summary(a)$coefficients[2,4]
    }else{
      a=lm(phe_pre[ids,i]~leaf_50[ids,j]+pca[ids,1]+pca[ids,2]+pca[ids,3]+pca[ids,4]+pca[ids,5])
      a=data.frame(summary(a)[4])
      est_mat[i,j]=a[2,1]
      p_mat[i,j]=a[2,4] 
    }
  }
}
row.names(p_mat)=names(phe_pre)
names(p_mat)=names(leaf_50)
row.names(est_mat)=names(phe_pre)
names(est_mat)=names(leaf_50)
leaf_mwas=list("Beta"=est_mat,"P"=p_mat)


data=root_mwas$P
for(i in 1:ncol(data)){
  for(j in 1:nrow(data)){
    data[j,i]=-log10(data[j,i])
  }
}
sub=data.frame(data)
names(sub)=names(root_phe)[3:ncol(root_phe)]
annotation_col=data.frame(Class=root_feature[names(sub),'p'])
row.names(annotation_col)=names(sub)
sub=sub[,order(root_feature[names(sub),'p'])]
pheatmap(sub,cluster_rows = F,cluster_cols = T,color= colorRampPalette(colors = c("white","pink","tomato","red"))(100),
         annotation_col = annotation_col,annotation_names_col = T,labels_col ="")

data=leaf_mwas$P
for(i in 1:ncol(data)){
  for(j in 1:nrow(data)){
    data[j,i]=-log10(data[j,i])
  }
}
sub=data.frame(data)
names(sub)=names(leaf_phe)[3:ncol(leaf_phe)]
annotation_col=data.frame(Class=leaf_feature[names(sub),'p'])
row.names(annotation_col)=names(sub)
sub=sub[,order(leaf_feature[names(sub),'p'])]
pheatmap(sub,cluster_rows = F,cluster_cols = F,,color= colorRampPalette(colors = c("white","pink","tomato","red"))(100),
         annotation_col = annotation_col,annotation_names_col = T,labels_col ="")



# micro PVE of phe --------------------------------------------------------
library(hglm)
### get orm
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

pves_fix_data=data.frame(trait=names(phe_pre),leaf_fix=0,root_fix=0,leaf_rand=0,root_rand=0,sum_var=0,var_phe=0,sum_xb=0,sum_zu=0)
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(phe_pre),intersect(row.names(leaf_orm),row.names(root_orm)))
  index_leaf=leaf_mwas$P[names(phe_pre)[i],]<0.001
  index_root=root_mwas$P[names(phe_pre)[i],]<0.001
  X=cbind(rep(1,length(ids)),leaf_50[ids,index_leaf],root_50[ids,index_root])
  X[is.na(X)]=0
  if(ncol(X)!=1){
    X=as.matrix(X)
    Z1=c.z.hglm(root_orm[ids,ids])
    Z2=c.z.hglm(leaf_orm[ids,ids])
    index=!is.na(phe_pre[,i])
    if(phe_info[names(phe_pre)[i],"class"]=="Categorical"){
      model=hglm(X=X[index,],y=phe_pre[ids,i][index],Z=cbind(Z1[index,index],Z2[index,index]),RandC=c(sum(index),sum(index)),family=binomial(link = "logit"))
    }else{
      model=hglm(X=X[index,],y=phe_pre[ids,i][index],Z=cbind(Z1[index,index],Z2[index,index]),RandC=c(sum(index),sum(index)))  
    }
    
    
    if(sum(index_leaf)==1){
      leaf_xb=var(X[index,2:(sum(index_leaf)+1)]*model$fixef[2:(sum(index_leaf)+1)]  )
    }else if(sum(index_leaf)==0){
      leaf_xb=0
    }else{
      leaf_xb=var(X[index,2:(sum(index_leaf)+1)]%*%model$fixef[2:(sum(index_leaf)+1)]  )
    }
    if(sum(index_root)==1){
      root_xb=var(X[index,(sum(index_leaf)+2):ncol(X)]*model$fixef[(sum(index_leaf)+2):ncol(X)]  )
    }else if(sum(index_root)==0){
      root_xb=0
    }else{
      root_xb=var(X[index,(sum(index_leaf)+2):ncol(X)]%*%model$fixef[(sum(index_leaf)+2):ncol(X)])
    }
    
    root_zu=var(Z1[index,index]%*%model$ranef[1:sum(index)])
    leaf_zu=var(Z2[index,index]%*%model$ranef[(sum(index)+1):(sum(index)*2)])
    
    pves_fix_data[i,"leaf_fix"]=leaf_xb
    pves_fix_data[i,"leaf_rand"]=leaf_zu
    pves_fix_data[i,"root_fix"]=root_xb
    pves_fix_data[i,"root_rand"]=root_zu
    pves_fix_data[i,"sum_var"]=var(X[index,]%*%model$fixef)+var(cbind(Z1[index,index],Z2[index,index])%*%model$ranef)
    pves_fix_data[i,"var_phe"]=var(phe_pre[ids,i][index])
    pves_fix_data[i,"sum_xb"]=var(X[index,]%*%model$fixef)
    pves_fix_data[i,"sum_zu"]=var(cbind(Z1[index,index],Z2[index,index])%*%model$ranef)
    
    next
  }

  Z1=c.z.hglm(root_orm[ids,ids])
  Z2=c.z.hglm(leaf_orm[ids,ids])
  index=!is.na(phe_pre[,i])
  if(phe_info[names(phe_pre)[i],"class"]=="Categorical"){
    model=hglm(X=matrix(rep(1,sum(index))),y=phe_pre[ids,i][index],Z=cbind(Z1[index,index],Z2[index,index]),RandC=c(sum(index),sum(index)),family=binomial(link = "logit"))
  }else{
    model=hglm(X=matrix(rep(1,sum(index))),y=phe_pre[ids,i][index],Z=cbind(Z1[index,index],Z2[index,index]),RandC=c(sum(index),sum(index)))  
  }
  
  leaf_xb=0
  root_xb=0
 
  root_zu=var(Z1[index,index]%*%model$ranef[1:sum(index)])
  leaf_zu=var(Z2[index,index]%*%model$ranef[(sum(index)+1):(sum(index)*2)])
  
  pves_fix_data[i,"leaf_fix"]=leaf_xb
  pves_fix_data[i,"leaf_rand"]=leaf_zu
  pves_fix_data[i,"root_fix"]=root_xb
  pves_fix_data[i,"root_rand"]=root_zu
  pves_fix_data[i,"sum_var"]=var(cbind(Z1[index,index],Z2[index,index])%*%model$ranef)
  pves_fix_data[i,"var_phe"]=var(phe_pre[ids,i][index])
  pves_fix_data[i,"sum_xb"]=0
  pves_fix_data[i,"sum_zu"]=1
}

data=pves_fix_data[,c(2,3,4,5)]
sum=c(apply(data,1,sum))
for(i in 1:4){
  data[,i]=data[,i]/sum
  #data[,i]=data[,i]/pves_fix_data$var_phe
}
data=data[,c(3,1,2,4)]

par(mfrow=c(2,1))
d1=apply(root_mwas$P,1,function(x){sum(x<0.05/27)})
d2=apply(leaf_mwas$P,1,function(x){sum(x<0.05/27)})
barplot(d1,col=rgb(0,0.5,1,0.2),xlim=c(0,40),ylab="Number of ASVs",names.arg = NA,ylim=c(0,10))
par(new = TRUE)
barplot(d2,col=rgb(1,0,0,0.4),xlim=c(0,40),ylab="",names.arg = NA,ylim=c(0,10))
legend("topright",c("Root MWAS ASVs","Leaf MWAS ASVs"),fill=c(rgb(0,0.5,1,0.2),rgb(1,0,0,0.4)),bty="n")

barplot(t(data),names.arg = pves_fix_data$trait,col=c("tomato","orange","#CCEBC5","skyblue"),
        xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("root rand","root fix","leaf fix","leaf rand"),fill=c("skyblue","#CCEBC5","orange","tomato"),bty="n")

names(pves_fix_data)
data=pves_fix_data[,c(11:12)]
sum=c(apply(data,1,sum))
for(i in 1:2){
  #data[,i]=data[,i]/sum
  data[,i]=data[,i]/pves_fix_data$sum_var2
}
d1=apply(root_mwas$P,1,function(x){sum(x<0.05/27)})
d2=apply(leaf_mwas$P,1,function(x){sum(x<0.05/27)})
barplot(d1,col=rgb(0,0.5,1,0.2),xlim=c(0,40),ylab="Number of ASVs",names.arg = NA,ylim=c(0,10))
par(new = TRUE)
barplot(d2,col=rgb(1,0,0,0.4),xlim=c(0,40),ylab="",names.arg = NA,ylim=c(0,10))
legend("topright",c("Root MWAS ASVs","Leaf MWAS ASVs"),fill=c(rgb(0,0.5,1,0.2),rgb(1,0,0,0.4)),bty="n")

barplot(t(data),names.arg = pves_fix_data$trait,col=c("tomato","skyblue"),
        xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("Rand","Fix"),fill=c("skyblue","tomato"),bty="n")


data=pves_fix_data[,c(8:9)]
sum=c(apply(data,1,sum))
for(i in 1:2){
  #data[,i]=data[,i]/sum
  data[,i]=data[,i]/pves_fix_data$var_phe
}

barplot(t(data),names.arg = pves_fix_data$trait,col=c("tomato","skyblue"),
        xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("Rand","Fix"),fill=c("skyblue","tomato"),bty="n")


# micro PVE of phe in two random effect -----------------------------------
pves_leaf_root=data.frame(matrix(ncol=2,nrow=27))
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(phe_pre),intersect(row.names(leaf_orm),row.names(root_orm)))
  Z1=c.z.hglm(root_orm[ids,ids])
  Z2=c.z.hglm(leaf_orm[ids,ids])
  index=!is.na(phe_pre[,i])
  if(phe_info[names(phe_pre)[i],"class"]=="Categorical"){
    model=hglm(X=matrix(rep(1,sum(index))),y=phe_pre[ids,i][index],Z=cbind(Z1[index,index],Z2[index,index]),RandC=c(sum(index),sum(index)),family=binomial(link = "logit"))
  }else{
    model=hglm(X=matrix(rep(1,sum(index))),y=phe_pre[ids,i][index],Z=cbind(Z1[index,index],Z2[index,index]),RandC=c(sum(index),sum(index)))  
  }
  pves_leaf_root[i,]=model$varRanef/(sum(model$varRanef)+model$varFix)
}
barplot(t(pves_leaf_root),names.arg = names(phe_pre),col=c("tomato","skyblue"),
        xlim=c(0,40),las=2,ylab="PVE",ylim=c(0,1))
legend("topright",c("leaf","root"),fill=c("skyblue","tomato"),bty="n")


