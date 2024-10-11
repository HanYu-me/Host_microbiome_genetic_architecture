library(data.table)
library(microeco)
library(ape)
library(magrittr)
library(ggplot2)
library(igraph)
library(pheatmap)
# Prepare input data of root ----------------------------------------------
otu_table=data.frame(t(root_50))
names(otu_table)=paste0("id_",row.names(root_50))
otu_table[is.na(otu_table)]=0

taxonomy_table_all=root_feature[sub("id_","",row.names(otu_table)),c(3:6)]
taxonomy_table_all%<>%tidy_taxonomy
names(taxonomy_table_all)=c("Kingdom","Phylum","Class","Order")
table(taxonomy_table_all$Kingdom)

sample_info=data.frame(SampleID=colnames(otu_table),Group="Root",stringsAsFactors = F)
row.names(sample_info)=sample_info$SampleID

phe_data_tax_root=phe_pre[sub("id_","",colnames(otu_table)),]
row.names(phe_data_tax_root)=names(otu_table)

data_micro_root=microtable$new(sample_table=sample_info,otu_table=otu_table,tax_table=taxonomy_table_all)
# Prepare input data of leaf ----------------------------------------------
otu_table=data.frame(t(leaf_50))
names(otu_table)=paste0("id_",row.names(leaf_50))
otu_table[is.na(otu_table)]=0

taxonomy_table_all=leaf_feature[sub("id_","",row.names(otu_table)),c(3:6)]
taxonomy_table_all%<>%tidy_taxonomy
names(taxonomy_table_all)=c("Kingdom","Phylum","Class","Order")
table(taxonomy_table_all$Kingdom)

sample_info=data.frame(SampleID=colnames(otu_table),Group="Leaf",stringsAsFactors = F)
row.names(sample_info)=sample_info$SampleID

phe_data_tax_leaf=phe_pre[sub("id_","",colnames(otu_table)),]
row.names(phe_data_tax_leaf)=names(otu_table)

data_micro_leaf=microtable$new(sample_table=sample_info,otu_table=otu_table,tax_table=taxonomy_table_all)


# Process data ------------------------------------------------------------
### resample of data to make sure the same abundance or reads across samples
range(data_micro_root$sample_sums())
data_micro_root$rarefy_samples(sample.size=8000)
range(data_micro_root$sample_sums())

range(data_micro_leaf$sample_sums())
data_micro_leaf$rarefy_samples(sample.size=9000)
range(data_micro_leaf$sample_sums())
### calculate abandunce in tax levels
#data_micro$cal_abund()
#class(data_micro$taxa_abund)


# Construct network of root-------------------------------------------------------
#t1_root <- trans_network$new(dataset = data_micro, cor_method = "pearson", use_WGCNA_pearson_spearman = TRUE)
#a=unlist(t1_root$res_cor_p[["p"]])
#table(a<0.05/600)
#t1_root$cal_network(COR_p_thres = 0.05,COR_p_adjust="bonferroni")
#t1_root$cal_module(method = "cluster_fast_greedy")
t1_root$cal_network(COR_p_thres = 0.05,COR_p_adjust="fdr")
t1_root <- trans_network$new(dataset = data_micro, cor_method = "bray",  filter_thres = 0)
a=unlist(t1_root$res_cor_p[["cor"]])
quantile(a,0.99)
abline(v=0.67,col="red")
t1_root$cal_network(COR_cut = 0.67)
t1_root$cal_module(method = "cluster_edge_betweenness")
table(as.character(t1_root$res_node_table[,7]))
#t1_root$cal_network_attr()
#t1_root$res_network_attr
#t1_root$save_network(filepath = "network.gexf")
#t1_root$plot_taxa_roles(use_type = 1)
#t1_root$plot_taxa_roles(use_type = 2)
t1_root$cal_eigen()

t2 <- trans_env$new(dataset = data_micro, add_data = phe_data_tax_root)
t2$cal_cor(add_abund_table = t1$res_eigen)
t2$plot_cor()
a=t2$res_cor
table(a[,2])
a=t2$res_cor$Pvalue

# Construct network of leaf-------------------------------------------------------
t1_leaf <- trans_network$new(dataset = data_micro_leaf, cor_method = "spearman", use_WGCNA_pearson_spearman = TRUE)
t1_leaf <- trans_network$new(dataset = data_micro_leaf, cor_method = "bray")
a=unlist(t1_leaf$res_cor_p[["cor"]])
a=a[!a==1]
quantile(a,0.99)
hist(a)
abline(v=0.67,col="red")
t1_leaf$cal_network(COR_cut = 0.67)
t1_leaf$cal_network(COR_p_thres = 0.05,COR_p_adjust="fdr")
t1_leaf$cal_module(method = "cluster_edge_betweenness")
t1_leaf$cal_network_attr()
t1_leaf$res_network_attr
t1_leaf$cal_eigen()
table(as.character(t1_leaf$res_node_table[,7]))
t2 <- trans_env$new(dataset = data_micro_leaf, add_data = phe_data_tax_leaf)
t2$cal_cor(add_abund_table = t1_leaf$res_eigen)
t2$plot_cor()


# make big data for leaf and root -----------------------------------------
all_id=intersect(row.names(root_50),row.names(leaf_50))
otu_table1=data.frame(t(root_50[all_id,]))
names(otu_table1)=paste0("id_",all_id)
otu_table1[is.na(otu_table1)]=0
row.names(otu_table1)=paste0("root_",row.names(otu_table1))

otu_table2=data.frame(t(leaf_50[all_id,]))
names(otu_table2)=paste0("id_",all_id)
otu_table2[is.na(otu_table2)]=0
row.names(otu_table2)=paste0("leaf_",row.names(otu_table2))

all_otu_table=rbind(otu_table1,otu_table2)

a=root_feature[sub("root_","",row.names(otu_table1)),c(3:6)]
b=leaf_feature[sub("leaf_","",row.names(otu_table2)),c(3:6)]
row.names(a)=paste0("root_",row.names(a))
row.names(b)=paste0("leaf_",row.names(b))
taxonomy_table_all=rbind(a,b)
taxonomy_table_all%<>%tidy_taxonomy
names(taxonomy_table_all)=c("Kingdom","Phylum","Class","Order")
taxonomy_table_all$Family=c(rep("Root",nrow(a)),rep("Leaf",nrow(b)))

sample_info=data.frame(SampleID=colnames(all_otu_table),Group="Data",stringsAsFactors = F)
row.names(sample_info)=sample_info$SampleID

phe_data_tax_all=phe_pre[all_id,]
row.names(phe_data_tax_all)=names(all_otu_table)
data_micro_all=microtable$new(sample_table=sample_info,otu_table=all_otu_table,tax_table=taxonomy_table_all)


range(data_micro_all$sample_sums())
data_micro_all$rarefy_samples(sample.size=22000)
range(data_micro_all$sample_sums())

# construct network -------------------------------------------------------
t1<- trans_network$new(dataset = data_micro_all, cor_method = "pearson", use_WGCNA_pearson_spearman = TRUE)
t1$cal_network(COR_p_thres = 0.05,COR_p_adjust="bonferroni")

# t1<- trans_network$new(dataset = data_micro_all, cor_method = "bray")
# a=unlist(t1$res_cor_p[["cor"]])
# a=a[!a==1]
# quantile(a,0.99)
# hist(a)
# abline(v=0.665,col="red")
# table(a>0.665)
# a=t1$res_cor_p[["cor"]]
# a[a<0.665]=0
# a[a>=0.665]=1
#pheatmap(a,cluster_rows = T,cluster_cols = T)
#t1$cal_network(COR_cut = 0.665,COR_p_thres = 1)
#t1$cal_network(COR_cut = 0.665,network_method = "gcoda")

t1$cal_module(method = "cluster_fast_greedy")
t1$cal_network_attr()
t1$res_network_attr
t1$cal_eigen()
t1$get_node_table()

all_node=t1$res_node_table
table(as.character(t1$res_node_table[,7]))
t2 <- trans_env$new(dataset = data_micro_all,add_data = phe_data_tax_all)
t2$cal_cor(add_abund_table = t1$res_eigen)
t2$res_cor[["AdjPvalue"]]=t2$res_cor[["Pvalue"]]
t2$plot_cor()

t1$save_network(filepath = "Data/network/WGCNA_pearson_bonferroni_greedy.gexf")



data=root_mwas$P
data=leaf_mwas$P
sub=t(apply(data,1,function(x){return(ifelse(x<0.05/41,1,0))}))
row.names(sub)=names(phe_pre)
pheatmap(sub,cluster_rows = F,cluster_cols = F)


t1$plot_taxa_roles(use_type = 1)
t1$plot_taxa_roles(use_type = 2)


t1$cal_sum_links(taxa_level = "Phylum")
t1$plot_sum_links(method = "circlize", transparency = 0.2, annotationTrackHeight = circlize::mm_h(c(5, 5)))

t1$cal_sum_links(taxa_level = "Family")
t1$plot_sum_links(method = "circlize", transparency = 0, annotationTrackHeight = circlize::mm_h(c(5, 5)))



t1$get_edge_table()
all_edge=t1$res_edge_table

table(grepl("root",all_edge$node1)&grepl("root",all_edge$node2))
table(grepl("root",all_edge$node1)&grepl("leaf",all_edge$node2))
table(grepl("leaf",all_edge$node1)&grepl("leaf",all_edge$node2))

all_node=t1$res_node_table


# net work pc analysis ----------------------------------------------------
table(all_node$module)
module_data=t2$res_cor
# convert to p value matrix,use dcast
module_data=as.data.table(module_data)
module_data_p=dcast(module_data,Taxa~Env,value.var="Pvalue")
module_data_b=dcast(module_data,Taxa~Env,value.var="Correlation")
row.names(module_data_b)=module_data_b$Taxa
module_data_b=module_data_b[,-1]
row.names(module_data_p)=module_data_p$Taxa
module_data_p=module_data_p[,-1]
module_data_p2=module_data_p
module_data_p2[module_data_p<=0.05]="*"
module_data_p2[module_data_p<=0.01]="**"
module_data_p2[module_data_p<=0.001]="***"
module_data_p2[module_data_p>0.05]=""
pheatmap(module_data_b,
         cluster_rows = F,cluster_cols = F,
         display_numbers = module_data_p2)
#add significance mark


# Contribution of network -------------------------------------------------
View(phe_pre)
View(all_otu_table)

library(parallel)
get_pve_fun=function(x){
  all_phe=phe_pre
  row.names(all_phe)=paste0("id_",row.names(all_phe))
  all_otu=data.frame(t(all_otu_table))
  all_pca=pca
  row.names(all_pca)=paste0("id_",row.names(all_pca))
  all_id=intersect(row.names(all_otu),row.names(all_phe))
  all_id=intersect(all_id,paste0("id_",row.names(pca)))
  ### get pve of all otu
  single_pve=unlist(lapply(names(all_otu),function(asv){
    a=lm(all_phe[all_id,x]~all_pca[all_id,1]+all_pca[all_id,2]+all_pca[all_id,3]+all_pca[all_id,4]+all_pca[all_id,5]+all_otu[all_id,asv])
    return(1-var(all_phe[all_id,x]-a$coefficients[7]*all_otu[all_id,asv],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
  }))
  
  net_pve=unlist(lapply(names(table(all_edge$node1))[table(all_edge$node1)>1],function(module){
    node=all_edge$node2[all_edge$node1==module]
    node=c(node,module)
    lm_data=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,node])
    names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
    a=lm(phe~.,data=lm_data)
    return(1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
  }))
  
  sample_pve=unlist(lapply(names(table(all_edge$node1))[table(all_edge$node1)>1],function(module){
    node=all_edge$node2[all_edge$node1==module]
    sample_num=length(node)+1
    sub_pve=c()
    for(i in 1:100){
      sample_node=sample(names(all_otu),sample_num,replace=F)
      lm_data=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,sample_node])
      names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
      a=lm(phe~.,data=lm_data)
      sub_pve=c(sub_pve,1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,sample_node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
    }
    sub_pve[sub_pve<0]=0
    return(mean(sub_pve))
  }))
  
  shuffle_pve=unlist(lapply(names(table(all_edge$node1))[table(all_edge$node1)>1],function(module){
    node=all_edge$node2[all_edge$node1==module]
    node=c(node,module)
    lm_data_sava=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,node])
    sub_pve=c()
    for(i in 1:100){
      lm_data=lm_data_sava
      lm_data[,node]=apply(lm_data[,node],2,function(xx){x=sample(xx)})
      names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
      a=lm(phe~.,data=lm_data)
      sub_pve=c(sub_pve,1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
    }
    sub_pve[sub_pve<0]=0
    return(mean(sub_pve))
  }))
  ### get net pve according the module in all_node
  net_pve2=unlist(lapply(names(table(all_node[,"module"])),function(module){
    node=all_node$name[all_node[,"module"]==module]
    lm_data=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,node])
    names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
    a=lm(phe~.,data=lm_data)
    return(1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
  }))
  
  sample_pve2=unlist(lapply(names(table(all_node[,"module"])),function(module){
    node=all_node$name[all_node[,"module"]==module]
    sample_num=length(node)
    sub_pve=c()
    for(i in 1:100){
      sample_node=sample(names(all_otu),sample_num,replace=F)
      lm_data=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,sample_node])
      names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
      a=lm(phe~.,data=lm_data)
      sub_pve=c(sub_pve,1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,sample_node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
    }
    sub_pve[sub_pve<0]=0
    return(mean(sub_pve))
  }))
  
  shuffle_pve2=unlist(lapply(names(table(all_node[,"module"])),function(module){
    node=all_node$name[all_node[,"module"]==module]
    lm_data_sava=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,node])
    sub_pve=c()
    for(i in 1:100){
      lm_data=lm_data_sava
      lm_data[,node]=apply(lm_data[,node],2,function(xx){x=sample(xx)})
      names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
      a=lm(phe~.,data=lm_data)
      sub_pve=c(sub_pve,1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
    }
    sub_pve[sub_pve<0]=0
    return(mean(sub_pve))
  }))
  
  return(list("single"=single_pve,"net"=net_pve,"sample"=sample_pve,"shuffle"=shuffle_pve,
              "net2"=net_pve2,"sample2"=sample_pve2,"shuffle2"=shuffle_pve2))
}
core=20
cl <- makeCluster(core)
clusterExport(cl,deparse(substitute(get_pve_fun)))
clusterExport(cl,c("phe_pre","all_edge","all_otu_table","pca","phe_info","all_node"))
net_pve_data=parLapply(cl,phe_info$name[phe_info$class=="Continues"],get_pve_fun)
stopCluster(cl)
names(net_pve_data)=phe_info$name[phe_info$class=="Continues"]
save(net_pve_data,file="Data/network/pve.RData")


# display result ----------------------------------------------------------
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gghalves) 
data_plot=data.frame(value=0,trait="",class="")[-1,]
png("net_pve.png",width=2000,height=4000)
par(mfrow=c(12,4))
for(i in names(net_pve_data)){
  sub=net_pve_data[[i]]
  sub=data.frame(value=unlist(sub),trait=i,class=rep(names(sub),unlist(lapply(sub,length))))
  plot(density(sub$value[sub$class=="single"]),xlim=c(0,0.25),col="black",main=i)
  lines(density(sub$value[sub$class=="net"]),col="red")
  lines(density(sub$value[sub$class=="sample"]),col="blue")
  lines(density(sub$value[sub$class=="shuffle"]),col="green")
  boxplot(value~class,data=sub,main=i)
  #sub=sub[sub$class!="sample",]
  data_plot=rbind(data_plot,sub)
}
dev.off()
#data_plot$class=factor(data_plot$class,levels=rev(c("net","net2","sample","sample2","shuffle","shuffle2","single")))
data_plot=data_plot[data_plot$class%in%c("net","single","shuffle"),]
data_plot$class=factor(data_plot$class,levels=rev(c("net","single","shuffle")))
g=ggplot(data_plot, aes(x = value, y = trait,fill = class)) +
  geom_density_ridges(alpha = .7, color = "white", from = 0, to = 0.2,scale=1)+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  theme_ridges()+ theme()
ggsave(g,filename="test.pdf",width=6,height=10)
table(data_plot$class)
## plot a boxplot about single and net and net2
ggplot(data_plot[data_plot$class%in%c("net","single","shuffle"),], aes(x = trait, y = value,fill = class)) +
  geom_boxplot(alpha=1,width=0.4,,outlier.color = "black",outlier.size=1)+
  scale_fill_manual(values=c(single="skyblue",net="tomato",shuffle="green"))+
  theme_bw()+ylim(c(-0.25,0.25))+ ylab("PVE")+xlab("")+
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

shuffle_pve=list()
for(x in phe_info$name[phe_info$class=="Continues"]){
  shuffle_pve[[x]]=unlist(lapply(names(table(all_edge$node1))[table(all_edge$node1)>1],function(module){
    node=all_edge$node2[all_edge$node1==module]
    node=c(node,module)
    lm_data_sava=data.frame(all_phe[all_id,x],all_pca[all_id,1],all_pca[all_id,2],all_pca[all_id,3],all_pca[all_id,4],all_pca[all_id,5],all_otu[all_id,node])
    sub_pve=c()
    for(i in 1:100){
      lm_data=lm_data_sava
      lm_data[,node]=apply(lm_data[,node],2,function(xx){x=sample(xx)})
      names(lm_data)[1:6]=c("phe","pc1","pc2","pc3","pc4","pc5")
      a=lm(phe~.,data=lm_data)
      sub_pve=c(sub_pve,1-var(all_phe[all_id,x]-as.matrix(all_otu[all_id,node])%*%a$coefficients[7:length(a$coefficients)],na.rm=T)/var(all_phe[all_id,x],na.rm=T))
      sub_pve[sub_pve<0]=0
    }
    return(mean(sub_pve))
  })) 
}
boxplot(shuffle_pve[[x]])

trait_pve=data.frame(row.names=c("Single",paste0("M",1:17)))
for(i in names(net_pve_data)){
  sub=net_pve_data[[i]]
  trait_pve[,i]=c(median(sub$single),net_pve2[[i]][1:17])
}
pheatmap(trait_pve,cluster_rows = F,cluster_cols = F)


View(all_node)
a=all_node[order(all_node$module),]
a=a[a$module%in%paste0("M",1:17),]
write.table(a,file="Data/network/module.txt",quote=F,sep="\t",row.names=F,col.names=T)
