load("Results//24.6.11_gwaresults.RData")


# plot phe manhattan ------------------------------------------------------
for(i in 1:ncol(phe_pre)){
  if(names(phe_pre)[i]%in%c("KC","TMV")){
    next
  }
  data=data.frame(fread(paste0("WAS_summary//phe_summary//sub_mlma//",names(phe_pre)[i],"_sub.txt")))
  for(j in 2:24){
    data[data[,1]==j,3]=as.numeric(data[data[,1]==j,3])+as.numeric(sum(faidx[1:(j-1),2]))
  }
  if(i==1){
    plot(x=data[,3],y=-log10(data[,9]),col=ifelse(data[,1]%%2==1,"tomato","skyblue"),
         frame.plot=F,pch=20,ylab="-Log10(p)",xlab="Chrosome",xaxt="n",ylim=c(3,20))  
  }else{
    points(x=data[,3],y=-log10(data[,9]),col=ifelse(data[,1]%%2==1,"tomato","skyblue"),pch=20) 
  }
}
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")
abline(h=-log10(5.76E-7),col="red")


phe_data=apply(phe_gwa_sum,2,sum)
plot(x=names(phe_data),y=phe_data,pch=20,cex=1.5,frame.plot=F,xaxt="n",
     col=ifelse(findInterval(names(leaf_data),chr_start)%%2==1,"tomato","skyblue"),
     xlab="Chromosome",ylab="Number of associated traits")
axis(side=1,at=chr_mean,1:24)
abline(v=c(chr_start,chr_end[24]),col="grey",lty="dashed")

library(ggvenn)
ggvenn(list(leaf=names(leaf_data[leaf_data>=1]),root=names(root_data[root_data>=1]),phe=names(phe_data[phe_data>=1])))  
ggvenn(list(leaf=names(leaf_data[leaf_data>1]),root=names(root_data[root_data>1]),phe=names(phe_data[phe_data>1])))  

