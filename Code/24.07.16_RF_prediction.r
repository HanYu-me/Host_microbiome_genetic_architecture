library(ranger)

pheatmap::pheatmap(snp_grm,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(leaf_orm,cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(root_orm,cluster_rows = F,cluster_cols = F)

ids=intersect(row.names(phe_man),intersect(intersect(names(root_orm),names(leaf_orm)),names(snp_grm)))
phe_all=phe_pre[ids,]
G1=snp_grm[ids,ids]
G2=leaf_orm[ids,ids]
G3=root_orm[ids,ids]
Z1=c.z.hglm(G1)
Z2=c.z.hglm(G2)
Z3=c.z.hglm(G3)
Z1=G1
Z2=G2
Z3=G3
idnum=length(ids)

RF_pre_all_10fold=lapply(1:27,function(i){
  phe=phe_all[,i]
  five_list=list()
  all_RF_joint=data.frame(cbind(phe,Z1,Z2,Z3))
  all_RF_only=data.frame(cbind(phe,Z1))
  names(all_RF_joint)[1]="phe"
  names(all_RF_only)[1]="phe"
  all_RF_joint=all_RF_joint[!is.na(all_RF_joint$phe),]
  all_RF_only=all_RF_only[!is.na(all_RF_only$phe),]
  #sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:5,c(nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5))))
  sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:10,rep(c(nrow(all_RF_joint)/10),10))))
  five_list[["phe_raw"]]=all_RF_joint$phe
  for(j in 1:10){
    test=sample_index[[j]]
    train=setdiff(1:nrow(all_RF_joint),test)
    if(phe_info[names(phe_man)[i],"class"]=="Categorical"){
      RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=1,classification=T)
      RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=1,classification=T)
    }else{
      #RF_joint=ranger(phe~.,data=all_RF_joint,num.threads=10,min.node.size=5,importance="impurity")
      RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=5)
      RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=5)
      
    }
    data=predict(RF_joint,all_RF_joint)
    five_list[[as.character(j)]]=list(three_pre=predict(RF_joint,all_RF_joint)$predictions,
                                      one_pre=predict(RF_only,all_RF_only)$predictions,
                                      index=test)
    #imp=data.frame(imp=RF_joint$variable.importance,type=rep(c("SNP","Leaf","Root"),each=nrow(all_RF_joint)))
    #boxplot(imp~type,imp)
  }
  return(five_list)
})

RF_pre_all_5fold=lapply(1:27,function(i){
  phe=phe_all[,i]
  five_list=list()
  all_RF_joint=data.frame(cbind(phe,Z1,Z2,Z3))
  all_RF_only=data.frame(cbind(phe,Z1))
  names(all_RF_joint)[1]="phe"
  names(all_RF_only)[1]="phe"
  all_RF_joint=all_RF_joint[!is.na(all_RF_joint$phe),]
  all_RF_only=all_RF_only[!is.na(all_RF_only$phe),]
  sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:5,c(nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5))))
  #sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:10,rep(c(nrow(all_RF_joint)/10),10))))
  five_list[["phe_raw"]]=all_RF_joint$phe
  for(j in 1:5){
    test=sample_index[[j]]
    train=setdiff(1:nrow(all_RF_joint),test)
    if(phe_info[names(phe_man)[i],"class"]=="Categorical"){
      RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=1,classification=T)
      RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=1,classification=T)
    }else{
      #RF_joint=ranger(phe~.,data=all_RF_joint,num.threads=10,min.node.size=5,importance="impurity")
      RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=5)
      RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=5)
      
    }
    data=predict(RF_joint,all_RF_joint)
    five_list[[as.character(j)]]=list(three_pre=predict(RF_joint,all_RF_joint)$predictions,
                                      one_pre=predict(RF_only,all_RF_only)$predictions,
                                      index=test)
    #imp=data.frame(imp=RF_joint$variable.importance,type=rep(c("SNP","Leaf","Root"),each=nrow(all_RF_joint)))
    #boxplot(imp~type,imp)
  }
  return(five_list)
})

names(RF_pre_all_10fold)=names(phe_pre)
names(RF_pre_all_5fold)=names(phe_pre)


###plot

data_plot_5fold=data.frame(trait=rep(names(RF_pre_all),each=10),class=rep(rep(c("SNP only","SNP and micro"),each=5),27),value=0)
data_plot_10fold=data.frame(trait=rep(names(RF_pre_all),each=20),class=rep(rep(c("SNP only","SNP and micro"),each=10),27),value=0)
for(i in names(RF_pre_all_10fold)){
  sub=RF_pre_all_10fold[[i]]
  c1=c()
  c2=c()
  if(phe_info[i,"class"]=="Continues"){
    for(j in 2:11){
      #c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
      #c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
      c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
      c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
    } 
  }else{
    for(j in 2:11){
      #c1=c(c1,sum(sub[[j]][[1]]==sub$phe_raw)/(sum(sub[[j]][[1]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
      #c2=c(c2,sum(sub[[j]][[2]]==sub$phe_raw)/(sum(sub[[j]][[2]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
      c1=c(c1,sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
      c2=c(c2,sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
    } 
  }
  data_plot_10fold[data_plot_10fold$trait==i&data_plot_10fold$class=="SNP only","value"]=c2
  data_plot_10fold[data_plot_10fold$trait==i&data_plot_10fold$class=="SNP and micro","value"]=c1
  print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
}
for(i in names(RF_pre_all_5fold)){
  sub=RF_pre_all_5fold[[i]]
  c1=c()
  c2=c()
  if(phe_info[i,"class"]=="Continues"){
    for(j in 2:6){
      #c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
      #c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
      c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
      c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
    } 
  }else{
    for(j in 2:6){
      #c1=c(c1,sum(sub[[j]][[1]]==sub$phe_raw)/(sum(sub[[j]][[1]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
      #c2=c(c2,sum(sub[[j]][[2]]==sub$phe_raw)/(sum(sub[[j]][[2]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
      c1=c(c1,sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
      c2=c(c2,sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
    } 
  }
  data_plot_5fold[data_plot_5fold$trait==i&data_plot_5fold$class=="SNP only","value"]=c2
  data_plot_5fold[data_plot_5fold$trait==i&data_plot_5fold$class=="SNP and micro","value"]=c1
  print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
}


library(ggplot2)
library(ggpubr)

data_plot_5fold$trait=factor(data_plot_5fold$trait,levels=names(phe_man))
data_plot_5fold$class=factor(data_plot_5fold$class,levels=c("SNP only","SNP and micro"))
ggplot(data_plot_5fold, aes(x = trait, y = value, fill = class)) +
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ ylab("R/AC")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
  stat_compare_means(method = "t.test")

data_plot_10fold$trait=factor(data_plot_10fold$trait,levels=names(phe_man))
data_plot_10fold$class=factor(data_plot_10fold$class,levels=c("SNP only","SNP and micro"))
ggplot(data_plot_10fold, aes(x = trait, y = value, fill = class)) +
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ ylab("R/AC")+xlab("")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
  stat_compare_means(method = "t.test")

unlist(lapply(RF_pre_all_5fold,function(x){length(x[[2]][[3]])}))
unlist(lapply(RF_pre_all_10fold,function(x){length(x[[2]][[3]])}))

library(patchwork)

# 5 fold together
{
  RF_pre_all_5fold=lapply(1:27,function(i){
    phe=phe_all[,i]
    five_list=list()
    all_RF_joint=data.frame(cbind(phe,Z1,Z2,Z3))
    all_RF_only=data.frame(cbind(phe,Z1))
    names(all_RF_joint)[1]="phe"
    names(all_RF_only)[1]="phe"
    all_RF_joint=all_RF_joint[!is.na(all_RF_joint$phe),]
    all_RF_only=all_RF_only[!is.na(all_RF_only$phe),]
    sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:5,c(nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5))))
    #sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:10,rep(c(nrow(all_RF_joint)/10),10))))
    five_list[["phe_raw"]]=all_RF_joint$phe
    for(j in 1:5){
      test=sample_index[[j]]
      train=setdiff(1:nrow(all_RF_joint),test)
      if(phe_info[names(phe_man)[i],"class"]=="Categorical"){
        RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=1,classification=T,num.trees=1000)
        RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=1,classification=T,num.trees=1000)
      }else{
        #RF_joint=ranger(phe~.,data=all_RF_joint,num.threads=10,min.node.size=5,importance="impurity")
        RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=5,num.trees=1000)
        RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=5,num.trees=1000)
        
      }
      data=predict(RF_joint,all_RF_joint)
      
      five_list[[as.character(j)]]=list(three_pre=predict(RF_joint,all_RF_joint)$predictions,
                                        one_pre=predict(RF_only,all_RF_only)$predictions,
                                        index=test)
      #imp=data.frame(imp=RF_joint$variable.importance,type=rep(c("SNP","Leaf","Root"),each=nrow(all_RF_joint)))
      #boxplot(imp~type,imp)
    }
    return(five_list)
  })
  names(RF_pre_all_5fold)=names(phe_pre)
  data_plot_5fold=data.frame(trait=rep(names(RF_pre_all),each=10),class=rep(rep(c("SNP only","SNP and micro"),each=5),27),value=0)
  #test and train
  for(i in names(RF_pre_all_5fold)){
    sub=RF_pre_all_5fold[[i]]
    c1=c()
    c2=c()
    if(phe_info[i,"class"]=="Continues"){
      for(j in 2:6){
        #c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
        #c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
        c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
        c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
      } 
    }else{
      for(j in 2:6){
        #c1=c(c1,sum(sub[[j]][[1]]==sub$phe_raw)/(sum(sub[[j]][[1]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        #c2=c(c2,sum(sub[[j]][[2]]==sub$phe_raw)/(sum(sub[[j]][[2]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        c1=c(c1,sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
        c2=c(c2,sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
      } 
    }
    data_plot_5fold[data_plot_5fold$trait==i&data_plot_5fold$class=="SNP only","value"]=c2
    data_plot_5fold[data_plot_5fold$trait==i&data_plot_5fold$class=="SNP and micro","value"]=c1
    print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
  }
  data_plot_5fold$trait=factor(data_plot_5fold$trait,levels=names(phe_man))
  data_plot_5fold$class=factor(data_plot_5fold$class,levels=c("SNP only","SNP and micro"))
  g1=ggplot(data_plot_5fold, aes(x = trait, y = value, fill = class)) +
    geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
                 size=0.75,outlier.color = "white",outlier.size=0.01)+
    theme_bw()+ ylab("test+train R/AC")+xlab("")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
    stat_compare_means(method = "t.test")
  
  data_plot_5fold=data.frame(trait=rep(names(RF_pre_all),each=10),class=rep(rep(c("SNP only","SNP and micro"),each=5),27),value=0)
  #test only
  for(i in names(RF_pre_all_5fold)){
    sub=RF_pre_all_5fold[[i]]
    c1=c()
    c2=c()
    if(phe_info[i,"class"]=="Continues"){
      for(j in 2:6){
        c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
        c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
        #c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
        #c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
      } 
    }else{
      for(j in 2:6){
        c1=c(c1,sum(sub[[j]][[1]]==sub$phe_raw)/(sum(sub[[j]][[1]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        c2=c(c2,sum(sub[[j]][[2]]==sub$phe_raw)/(sum(sub[[j]][[2]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        #c1=c(c1,sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
        #c2=c(c2,sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
      } 
    }
    data_plot_5fold[data_plot_5fold$trait==i&data_plot_5fold$class=="SNP only","value"]=c2
    data_plot_5fold[data_plot_5fold$trait==i&data_plot_5fold$class=="SNP and micro","value"]=c1
    print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
  }
  data_plot_5fold$trait=factor(data_plot_5fold$trait,levels=names(phe_man))
  data_plot_5fold$class=factor(data_plot_5fold$class,levels=c("SNP only","SNP and micro"))
  g2=ggplot(data_plot_5fold, aes(x = trait, y = value, fill = class)) +
    geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
                 size=0.75,outlier.color = "white",outlier.size=0.01)+
    theme_bw()+ ylab("test only R/AC")+xlab("")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
    stat_compare_means(method = "t.test")  
}
g1/g2

#10 fold together
{
  RF_pre_all_10fold=lapply(1:27,function(i){
    phe=phe_all[,i]
    five_list=list()
    all_RF_joint=data.frame(cbind(phe,Z1,Z2,Z3))
    all_RF_only=data.frame(cbind(phe,Z1))
    names(all_RF_joint)[1]="phe"
    names(all_RF_only)[1]="phe"
    all_RF_joint=all_RF_joint[!is.na(all_RF_joint$phe),]
    all_RF_only=all_RF_only[!is.na(all_RF_only$phe),]
    #sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:5,c(nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5,nrow(all_RF_joint)/5))))
    sample_index=split(c(1:nrow(all_RF_joint)),sample(rep(1:10,rep(c(nrow(all_RF_joint)/10),10))))
    five_list[["phe_raw"]]=all_RF_joint$phe
    for(j in 1:10){
      test=sample_index[[j]]
      train=setdiff(1:nrow(all_RF_joint),test)
      if(phe_info[names(phe_man)[i],"class"]=="Categorical"){
        RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=1,classification=T,num.trees=1000)
        RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=1,classification=T,num.trees=1000)
      }else{
        #RF_joint=ranger(phe~.,data=all_RF_joint,num.threads=10,min.node.size=5,importance="impurity")
        RF_joint=ranger(phe~.,data=all_RF_joint[train,],num.threads=10,min.node.size=5,num.trees=1000)
        RF_only=ranger(phe~.,data=all_RF_only[train,],num.threads=10,min.node.size=5,num.trees=1000)
        
      }
      data=predict(RF_joint,all_RF_joint)
      five_list[[as.character(j)]]=list(three_pre=predict(RF_joint,all_RF_joint)$predictions,
                                        one_pre=predict(RF_only,all_RF_only)$predictions,
                                        index=test)
      #imp=data.frame(imp=RF_joint$variable.importance,type=rep(c("SNP","Leaf","Root"),each=nrow(all_RF_joint)))
      #boxplot(imp~type,imp)
    }
    return(five_list)
  })
  names(RF_pre_all_10fold)=names(phe_pre)
  data_plot_10fold=data.frame(trait=rep(names(RF_pre_all),each=20),class=rep(rep(c("SNP only","SNP and micro"),each=10),27),value=0)
  for(i in names(RF_pre_all_10fold)){
    sub=RF_pre_all_10fold[[i]]
    c1=c()
    c2=c()
    if(phe_info[i,"class"]=="Continues"){
      for(j in 2:11){
        c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
        c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
        #c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
        #c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
      } 
    }else{
      for(j in 2:11){
        c1=c(c1,sum(sub[[j]][[1]]==sub$phe_raw)/(sum(sub[[j]][[1]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        c2=c(c2,sum(sub[[j]][[2]]==sub$phe_raw)/(sum(sub[[j]][[2]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        #c1=c(c1,sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
        #c2=c(c2,sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
      } 
    }
    data_plot_10fold[data_plot_10fold$trait==i&data_plot_10fold$class=="SNP only","value"]=c2
    data_plot_10fold[data_plot_10fold$trait==i&data_plot_10fold$class=="SNP and micro","value"]=c1
    print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
  }
  
  data_plot_10fold$trait=factor(data_plot_10fold$trait,levels=names(phe_man))
  data_plot_10fold$class=factor(data_plot_10fold$class,levels=c("SNP only","SNP and micro"))
  g1=ggplot(data_plot_10fold, aes(x = trait, y = value, fill = class)) +
    geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
                 size=0.75,outlier.color = "white",outlier.size=0.01)+
    theme_bw()+ ylab("test+train R/AC")+xlab("")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
    stat_compare_means(method = "t.test")

  data_plot_10fold=data.frame(trait=rep(names(RF_pre_all),each=20),class=rep(rep(c("SNP only","SNP and micro"),each=10),27),value=0)
  for(i in names(RF_pre_all_10fold)){
    sub=RF_pre_all_10fold[[i]]
    c1=c()
    c2=c()
    if(phe_info[i,"class"]=="Continues"){
      for(j in 2:11){
        #c1=c(c1,cor.test(sub$phe_raw,sub[[j]][[1]])$estimate)
        #c2=c(c2,cor.test(sub$phe_raw,sub[[j]][[2]])$estimate)
        c1=c(c1,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[1]][sub[[j]][[3]]])$estimate)
        c2=c(c2,cor.test(sub$phe_raw[sub[[j]][[3]]],sub[[j]][[2]][sub[[j]][[3]]])$estimate)
      } 
    }else{
      for(j in 2:11){
        #c1=c(c1,sum(sub[[j]][[1]]==sub$phe_raw)/(sum(sub[[j]][[1]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        #c2=c(c2,sum(sub[[j]][[2]]==sub$phe_raw)/(sum(sub[[j]][[2]]==sub$phe_raw)+sum(sub[[j]][[1]]!=sub$phe_raw)))
        c1=c(c1,sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[1]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
        c2=c(c2,sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])/(sum(sub[[j]][[2]][sub[[j]][[3]]]==sub$phe_raw[sub[[j]][[3]]])+sum(sub[[j]][[1]][sub[[j]][[3]]]!=sub$phe_raw[sub[[j]][[3]]])))
      } 
    }
    data_plot_10fold[data_plot_10fold$trait==i&data_plot_10fold$class=="SNP only","value"]=c2
    data_plot_10fold[data_plot_10fold$trait==i&data_plot_10fold$class=="SNP and micro","value"]=c1
    print(paste0(i,t.test(c2,c1,alternative = "less")$p.value))
  }
  
  data_plot_10fold$trait=factor(data_plot_10fold$trait,levels=names(phe_man))
  data_plot_10fold$class=factor(data_plot_10fold$class,levels=c("SNP only","SNP and micro"))
  g2=ggplot(data_plot_10fold, aes(x = trait, y = value, fill = class)) +
    geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
                 size=0.75,outlier.color = "white",outlier.size=0.01)+
    theme_bw()+ ylab("test only R/AC")+xlab("")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_fill_manual(values=c("SNP and micro"="#FDE3E1","SNP only"="#E4EECB"))+
    stat_compare_means(method = "t.test")
}
g1/g2

data[names(RF_pre_all)]
ids=intersect(row.names(root_50),row.names(leaf_50))
phe_p=data.frame(fread("Data//phe_filtered//Phe_all_raw.txt"))
dim(phe_p)
phe_p=phe_p[!is.na(phe_p[,1]),]
row.names(phe_p)=phe_p[,1]

phe_p=phe_p[ids,]
data=apply(phe_p,2,function(x){sum(!is.na(x))})
barplot(data[-1:-2],las=2,ylim=c(0,200))
