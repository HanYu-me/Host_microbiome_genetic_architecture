
# GBLUP 2-5fold 100 sample test -------------------------------------------
library(parallel)

phe_five_fold_GS_fun=function(i){
  fold=3
  ids=intersect(row.names(phe_man),intersect(intersect(names(root_orm),names(leaf_orm)),names(snp_grm)))
  phe_all=phe_pre[ids,]
  G1=snp_grm[ids,ids]
  G2=leaf_orm[ids,ids]
  G3=root_orm[ids,ids]
  idnum=length(ids)
  sample_index=split(c(1:(idnum)),sample(rep(1:fold,rep(idnum/fold,fold))))
  phe=phe_all[,i]
  phe_save=phe
  five_list=list()
  five_list[["phe_raw"]]=phe_save
  for(j in 1:fold){
    phe=phe_save
    phe[sample_index[[j]]]=NA
    index=!is.na(phe)
    Z1=c.z.hglm(G1)
    Z2=c.z.hglm(G2)
    Z3=c.z.hglm(G3)
    
    model=try(hglm(X=matrix(rep(1,length(phe[index]))),
               y=phe[index],
               Z=cbind(Z1[index,],Z2[index,],Z3[index,]),
               RandC=c(ncol(Z1[index,]),ncol(Z2[index,]),ncol(Z3[index,]))),silent = T)
    model2=try(hglm(X=matrix(rep(1,length(phe[index]))),
                y=phe[index],
                Z=Z1[index,]),silent = T)
    
    pre_all=try(cbind(Z1,Z2,Z3)%*%model$ranef,silent = T)
    pre_all2=try(Z1%*%model2$ranef,silent = T)
    five_list[[as.character(j)]]=try(list(three_pre=pre_all,one_pre=pre_all2,index=sample_index[[j]]),silent = T)
    
    info=try(paste("Times",times,"Trait",names(phe_all)[i],"Fold",fold,"Round",j,
               "R_G",cor(phe_save[sample_index[[j]]],pre_all2[sample_index[[j]]],use="pairwise.complete.obs"),
               "R_G+M",cor(phe_save[sample_index[[j]]],pre_all[sample_index[[j]]],use="pairwise.complete.obs")),silent = T)
    write.table(info,"Data/GS/3fold_GS_summary.txt",append=T,col.names=F,row.names=F,quote=F)
  }
  write.table(paste(times,c(i,names(phe_man)[i]),"finished"),"info.txt",append=T,col.names=F,row.names=F)
  return(five_list)
}
for(times in 1:100){
  print(times)
  clnum<-5
  cl <- makeCluster(getOption("cl.cores", clnum));
  clusterExport(cl,deparse(substitute(phe_five_fold_GS_fun)))
  clusterEvalQ(cl,library(hglm))
  clusterExport(cl,"c.z.hglm")
  clusterExport(cl,c("phe_man","phe_pre","snp_grm","leaf_orm","root_orm","phe_info","times"))
  result_1_5fold=parLapply(cl,1:22, phe_five_fold_GS_fun)
  save(result_1_5fold,file=paste0("Data/GS/3fold_GS/3fold_GS_",times,".RData"))
  stopCluster(cl)
}

save.image(".RData")
### summary
result=data.frame(fread("Data/GS/5fold_GS_summary.txt",fill=T))
library(stringr)
result=result[,c(2,4,8,10,12)]
names(result)=c("Times","Trait","Round","R_G","R_GM")
result=result[!result$Trait=="",]
result$R_G=abs(result$R_G)
result$R_GM=abs(result$R_GM)

#split result into list object by the Trait
result_list=lapply(split(result,result$Trait),function(x) x[,c("Times","Round","R_G","R_GM")])

library(ggplot2)
library(ggpubr)
data_plot=data.frame("Trait"=c(result$Trait,result$Trait),value=c(result$R_G,result$R_GM),
                     class=c(rep("SNP only",nrow(result)),rep("SNP and microbiome",nrow(result))))
data_plot$Trait=factor(data_plot$Trait,levels=names(phe_man))
data_plot$class=factor(data_plot$class,levels=c("SNP only","SNP and microbiome"))
ggplot(data_plot, aes(x = Trait, y = value, fill = class)) +
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ ylab("Pearson R")+xlab("")+ylim(c(0,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c("SNP and microbiome"="#FDE3E1","SNP only"="#E4EECB"))+
  stat_compare_means(method = "t.test")


# 5 A Prediction -----------------------------------------------------------------
data_plot=data.frame("Trait"=c(result$Trait,result$Trait),value=c(result$R_G,result$R_GM),
                     class=c(rep("SNP only",nrow(result)),rep("SNP and microbiome",nrow(result))))
data_plot$Trait=factor(data_plot$Trait,levels=names(phe_man))
data_plot$class=factor(data_plot$class,levels=c("SNP only","SNP and microbiome"))
data_plot1=data_plot[data_plot$Trait%in%names(phe_man)[1:11],]
g=ggplot(data_plot1, aes(x = Trait, y = value, fill = class)) +
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ ylab("Pearson R")+xlab("")+ylim(c(0,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c("SNP and microbiome"="tomato","SNP only"="skyblue"))+
  stat_compare_means(method = "t.test",label="p.format")
ggsave(g,file="Writing//V1_24.09.29//Fig5//5A_1.pdf",width=15,height=5)

data_plot2=data_plot[data_plot$Trait%in%names(phe_man)[12:22],]
g=ggplot(data_plot2, aes(x = Trait, y = value, fill = class)) +
  geom_boxplot(alpha=1,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ ylab("Pearson R")+xlab("")+ylim(c(0,1))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_fill_manual(values=c("SNP and microbiome"="tomato","SNP only"="skyblue"))+
  stat_compare_means(method = "t.test",label="p.format")
ggsave(g,file="Writing//V1_24.09.29//Fig5//5A_2.pdf",width=15,height=5)
