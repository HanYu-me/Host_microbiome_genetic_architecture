
# 3A mantel overall relationship -------------------------------------------------
library(ggcor)
library(linkET)
library(vegan)
library(data.table)
phe_info=data.frame(fread("Data/anno/Trait_in_research.csv"))
names(phe_info)[3]="name"
phe_info=phe_info[,c(2:5)]
row.names(phe_info)=phe_info$name

phe_pre=data.frame(fread("Data//phe_filtered//Phe_all_raw.txt"))
phe_pre=phe_pre[!is.na(phe_pre[,1]),]
row.names(phe_pre)=phe_pre[,1]
ids=intersect(row.names(root_50),intersect(row.names(leaf_50),row.names(snp_grm)))
phe_pre=phe_pre[ids,]
phe_pre=phe_pre[,apply(phe_pre,2,function(x)sum(!is.na(x)))>100]
phe_pre=phe_pre[,c(-1,-2)]
phe_pre=phe_pre[,phe_info[names(phe_pre),"class"]=="Continues"]

phe_info=phe_info[names(phe_pre),]

all_id=ids
phe_man=phe_pre[as.character(all_id),rev(phe_info$name[order(phe_info$type)])]

asv_all=cbind(root_phe[as.character(all_id),3:ncol(root_phe)],leaf_phe[as.character(all_id),3:ncol(leaf_phe)])
length(3:ncol(root_phe))
length(3:ncol(leaf_phe))
asv_all[is.na(asv_all)]=0

for(i in 1:ncol(phe_man)){
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
mantel$env=factor(mantel$env,levels = phe_info$name)
library(ggplot2)
library(aplot)
library(linkET)
library(ggcor)
g = 
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
g

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

g=g%>%insert_left(group, width = 0.02)
ggsave(g,filename="Writing//V1_24.09.29//Fig3//3A.pdf",width=12,height=12)

# 3B MWAS -----------------------------------------------------------------
library(RColorBrewer)
leaf_p=data.frame(leaf_mwas$P)
leaf_p=leaf_p[names(phe_man),]
names(leaf_p)=names(leaf_50)
leaf_b=data.frame(leaf_mwas$Beta)
leaf_b=leaf_b[names(phe_man),]
names(leaf_b)=names(leaf_50)

leaf_annotation_col=data.frame(ASV=names(leaf_50),Class=root_feature[names(leaf_p),'p'])
leaf_class=names(rev(sort(table(leaf_annotation_col$Class))))
names(leaf_class)=leaf_class


root_p=data.frame(root_mwas$P)
root_p=root_p[names(phe_man),]
names(root_p)=names(root_50)
root_b=data.frame(root_mwas$Beta)
root_b=root_b[names(phe_man),]
names(root_b)=names(root_50)

root_annotation_col=data.frame(ASV=names(root_50),Class=root_feature[names(root_p),'p'])
root_class=names(rev(sort(table(root_annotation_col$Class))))
names(root_class)=root_class

root_class[1:24]=colorRampPalette(brewer.pal(11, "Spectral"))(length(root_class))
leaf_class=root_class[names(leaf_class)]

### plot leaf MWAS result
order_leaf=c()
for(i in names(leaf_class)){
  order_leaf=c(order_leaf,leaf_annotation_col[which(leaf_annotation_col$Class==i),"ASV"])
}
row.names(leaf_annotation_col)=leaf_annotation_col$ASV
col_leaf=leaf_class[leaf_annotation_col[order_leaf,"Class"]]
xinfo=data.frame(cbind(col_leaf,1:96))
names(xinfo)=c("col","x")
xinfo$x=as.numeric(xinfo$x)
xinfo$col=factor(xinfo$col,levels=leaf_class)
shape=c(15,8,16,17,18)
names(shape)=unique(phe_info$type)
pdf("Writing//V1_24.09.29//Fig3//3B.pdf",width=10,height=4)
plot(x=1:96,y=-log10(leaf_p[1,order_leaf]),col=col_leaf,pch=shape[phe_info[row.names(leaf_p)[1],"type"]],
     xaxt="n",yaxt="n",xlab="",ylab="-log10(P)",ylim=c(-0.5,5),frame.plot=F)
axis(side=1,at=aggregate(x~col,data=xinfo,FUN=mean)[,2],
     labels=rep("|",length(leaf_class)),las=1,tick=F,
     pos=0.6)
axis(side=2,at=0:5,labels=0:5,las=2)
#add rect into plot with color in col_leaf
leaf_str=aggregate(x~col,data=xinfo,FUN=min)[,2]
leaf_end=aggregate(x~col,data=xinfo,FUN=max)[,2]
for(i in 1:length(leaf_class)){
  rect(leaf_str[i]-0.5,-0.6,leaf_end[i]+0.5,-0.25,col=leaf_class[i],border=NA)
}
for(i in 2:nrow(leaf_p)){
  points(x=1:96,y=-log10(leaf_p[i,order_leaf]),col=col_leaf,pch=shape[phe_info[row.names(leaf_p)[i],"type"]])
}
abline(h=-log10(0.05/96),col="red",lty="dashed")
dev.off()

### plot root MWAS result
order_root=c()
for(i in names(root_class)){
  order_root=c(order_root,root_annotation_col[which(root_annotation_col$Class==i),"ASV"])
}
row.names(root_annotation_col)=root_annotation_col$ASV
col_root=root_class[root_annotation_col[order_root,"Class"]]
xinfo=data.frame(cbind(col_root,1:685))
names(xinfo)=c("col","x")
xinfo$col=factor(xinfo$col,levels=root_class)
xinfo$x=as.numeric(xinfo$x)
shape=c(15,8,16,17,18)
names(shape)=unique(phe_info$type)
pdf("Writing//V1_24.09.29//Fig3//3C.pdf",width=10,height=4)
plot(x=1:685,y=-log10(root_p[1,order_root]),col=col_root,pch=shape[phe_info[row.names(root_p)[1],"type"]],
     xaxt="n",yaxt="n",xlab="",ylab="-log10(P)",ylim=c(-0.5,5),frame.plot=F)
axis(side=1,at=aggregate(x~col,data=xinfo,FUN=mean)[,2],
     labels=rep("|",length(root_class)),las=1,tick=F,
     pos=0.6)
axis(side=2,at=0:5,labels=0:5,las=2)
#add rect into plot with color in col_root
root_str=aggregate(x~col,data=xinfo,FUN=min)[,2]
root_end=aggregate(x~col,data=xinfo,FUN=max)[,2]
for(i in 1:length(root_class)){
  rect(root_str[i]-0.5,-0.6,root_end[i]+0.5,-0.25,col=root_class[i],border=NA)
}
for(i in 2:nrow(root_p)){
  points(x=1:685,y=-log10(root_p[i,order_root]),col=col_root,pch=shape[phe_info[row.names(root_p)[i],"type"]])
}
abline(h=-log10(0.05/685),col="red",lty="dashed")
dev.off()

##give legend
#plot a empty plot
pdf("Writing//V1_24.09.29//Fig3//3BC_legend.pdf",width=20,height=5)
plot(100,100,xlim=c(0,1),ylim=c(0,1),col="white",xaxt="n",yaxt="n",xlab="",ylab="")
legend("topleft",sub(" p__","",names(root_class)),fill=root_class,border=NA,ncol=6)
dev.off()
pdf("Writing//V1_24.09.29//Fig3//3BC_legend2.pdf",width=20,height=5)
plot(100,100,xlim=c(0,1),ylim=c(0,1),col="white",xaxt="n",yaxt="n",xlab="",ylab="")
legend("topleft",names(shape),pch=shape,col="black",border=NA,ncol=5)
dev.off()

### give micro and phe information
mwas_info=data.frame(ASV="",Microbiome="",Trait="",P="",Beta="",Taxa="",stringsAsFactors = FALSE)[-1,]
for(i in 1:ncol(leaf_p)){
 for(j in 1:nrow(leaf_p)){
   if(is.na(leaf_p[j,i])){
     next
   }
   if(leaf_p[j,i]<0.05/96){
     
     mwas_info=rbind(mwas_info,data.frame(ASV=names(leaf_p)[i],Microbiome="Phyllosphere",Trait=rownames(leaf_p)[j],P=leaf_p[j,i],Beta=leaf_b[j,i],taxa=leaf_feature[names(leaf_p)[i],2]))
   }
 }
}
for(i in 1:ncol(root_p)){
 for(j in 1:nrow(root_p)){
   if(is.na(root_p[j,i])){
     next
   }
   if(root_p[j,i]<0.05/685){
     mwas_info=rbind(mwas_info,data.frame(ASV=names(root_p)[i],Microbiome="Rhizosphere",Trait=rownames(root_p)[j],P=root_p[j,i],Beta=root_b[j,i],taxa=root_feature[names(root_p)[i],2]))
   }
 }
}
write.table(mwas_info,"Writing//V1_24.09.29//STable//S6.txt",quote = F,row.names = F,col.names = T,sep = "\t")
table(phe_info[mwas_info$Trait,"type"])
# 3E netwrok pve ----------------------------------------------------------
library(microeco)
library(ggridges)
library(ggpubr)
library(gghalves)
data_plot=data.frame(value=0,trait="",class="")[-1,]
for(i in names(net_pve_data)){
  sub=net_pve_data[[i]]
  sub=data.frame(value=unlist(sub),trait=i,class=rep(names(sub),unlist(lapply(sub,length))))
  # plot(density(sub$value[sub$class=="single"]),xlim=c(0,0.25),col="black",main=i)
  # lines(density(sub$value[sub$class=="net"]),col="red")
  # lines(density(sub$value[sub$class=="sample"]),col="blue")
  # lines(density(sub$value[sub$class=="shuffle"]),col="green")
  #boxplot(value~class,data=sub,main=i)
  data_plot=rbind(data_plot,sub)
}
data_plot=data_plot[data_plot$class%in%c("net","single","shuffle"),]
data_plot$class=factor(data_plot$class,levels=rev(c("net","single","shuffle")))
data_plot[data_plot$value<0,"value"]=0
data_plot=data_plot[data_plot$trait%in%names(phe_man)[1:11],]
g=ggplot(data_plot, aes(x = value, y = trait,fill = class)) +
  geom_density_ridges(alpha = .7, color = "white", from = 0, to = 0.15,scale=1)+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  scale_fill_manual(values = c("skyblue","tomato","orange"))+
  theme_ridges()+ theme()
g
ggsave(g,filename="Writing//V1_24.09.29//Fig3//3E.pdf",width=6,height=10)

data_plot=data.frame(value=0,trait="",class="")[-1,]
for(i in names(net_pve_data)){
  sub=net_pve_data[[i]]
  sub=data.frame(value=unlist(sub),trait=i,class=rep(names(sub),unlist(lapply(sub,length))))
  data_plot=rbind(data_plot,sub)
}
data_plot=data_plot[data_plot$class%in%c("net","single","shuffle"),]
data_plot$class=factor(data_plot$class,levels=rev(c("net","single","shuffle")))
data_plot[data_plot$value<0,"value"]=0
data_plot=data_plot[data_plot$trait%in%names(phe_man)[12:22],]
g=ggplot(data_plot, aes(x = value, y = trait,fill = class)) +
  geom_density_ridges(alpha = .7, color = "white", from = 0, to = 0.15,scale=1)+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  scale_fill_manual(values = c("skyblue","tomato","orange"))+
  theme_ridges()+ theme()
g
ggsave(g,filename="Writing//V1_24.09.29//Fig3//3F.pdf",width=6,height=10)


# S fig S table -----------------------------------------------------------
### phe distribution
png("Writing//V1_24.09.29//SFig//S3.png",width=10,height=7.5,unit="in",res=900)
par(mfrow=c(4,6))
par(mar=c(4,4,1,1))
for(i in 1:ncol(phe_pre)){
  hist(phe_pre[,i],main=names(phe_pre)[i],xlab="",col="skyblue",border="white",las=2)
}
dev.off()
### get phe data table
library(hglm)
h2=c()
for(i in 1:ncol(phe_pre)){
  print(i)
  ids=intersect(row.names(phe_man),intersect(intersect(names(root_orm),names(leaf_orm)),names(snp_grm)))
  phe_all=phe_man[ids,]
  G1=snp_grm[ids,ids]
  idnum=length(ids)
  phe=phe_pre[,i]

  index=!is.na(phe)
  Z1=c.z.hglm(G1)
  model2=hglm(X=matrix(rep(1,length(phe[index]))),
              y=phe[index],
              Z=Z1[index,])  
  h2=c(h2,model2$varRanef/(model2$varRanef+model2$varFix))
}
names(h2)=names(phe_pre)
data=cbind(phe_info[names(phe_man),c(2,1,3)],h2[names(phe_man)])
write.table(data,"Writing//V1_24.09.29//STable//S5.txt",quote = F,row.names = F,col.names = T,sep = "\t")

### net work society
t1$cal_sum_links(taxa_level = "Phylum")
pdf("Writing//V1_24.09.29//SFig//S4_1.pdf",width=5,height=5)
t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))
dev.off()

t1$cal_sum_links(taxa_level = "Family")
pdf("Writing//V1_24.09.29//SFig//S4_2.pdf",width=5,height=5)
t1$plot_sum_links(method = "circlize", transparency = 0, annotationTrackHeight = circlize::mm_h(c(5, 5)))
dev.off()

### network module 
t1$res_network_attr
View(all_node)
a=all_node[order(all_node$module),]
head(a)
a$Microbiom=ifelse(grepl("root",rownames(a)),"Rhizosphere","Phyllosphere")
a$ASV=sub("leaf_","",sub("root_","",rownames(a)))
names(a)
a=a[,c(18,17,7,2,11:14)]
length(table(a$module))
sort(table(a$module))
sum(sort(table(a$module))[22:29])
write.table(a,file="Writing//V1_24.09.29//STable//S7.txt",quote=F,sep="\t",row.names=F,col.names=T)
