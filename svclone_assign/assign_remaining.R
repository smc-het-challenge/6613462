#!/usr/bin/env Rscript

get_likely_clonal<-function(psi,pi1,nv,muv,n0,pi0,ov,bv)
{
  epsilon = (psi*pi1*nv*muv)/(n0*pi0+pi1*nv)
  dbinom(bv,bv+ov,epsilon,log=T)
}

get_likely_subclonal<-function(psi,pi1,nv,muv,n0,pi0,ov,bv,nt,prop)
{
  epsilon = (psi*pi1*nv*muv*prop)/(n0*pi0+pi1*(nv*prop+nt*(1-prop)))
  dbinom(bv,bv+ov,epsilon,log=T)
}

args = commandArgs(trailingOnly=TRUE)

fulldat<-read.table(args[1],sep="\t",header=T,stringsAsFactors = F)
input<-read.table(args[2],sep="\t",header=F,stringsAsFactors = F)
sampdat<-read.table(args[3],sep="\t",header=T,stringsAsFactors = F)
clustassign<-read.table(args[4],sep="\t",header=F)
colnames(clustassign)<-"clust"
cluststruct<-read.table(args[5],sep="\t",header=F)
purity<-as.numeric(read.table(args[6]))
clonal<-which.min(abs(cluststruct[,3]-purity))

#cluststruct<-cluststruct[1:4,]
cluststruct<-cluststruct[cluststruct[,2]/sum(cluststruct[,2])>0.02,]
write.table(nrow(cluststruct),"./smc_1B_number_of_clusters.txt",sep="\t",row.names=F,col.names = F)


dat<-cbind(sampdat[,c(1,2)],clustassign)

#join input with full dat to get correct order
workingdat<-merge(input[,c(1,2,3)],fulldat,by.x = c(1,2),by.y=c(1,2),all.x=T)

#take input SNVs and left join with current assignments
workingdat<-merge(workingdat,dat,by.x = c(1,2),by.y=c(1,2),all.x=T)
colnames(workingdat)<-c("chr","pos","SNP","gtype","ref","var","clust")

workingdat$clust<-NA
pvalthresh<-0.05


for(i in 1:nrow(workingdat))
{
  curr<-workingdat[i,]
  gtypes<-unlist(strsplit(workingdat[i,]$gtype,"|",fixed=T))
  if(is.na(curr$gtype)){
    workingdat[i,]$clust<-0
  }else if(is.na(curr$clust)){
    if (length(gtypes)==2){
      #subclonal copy-number
      llclust<-c()
      for(l in 1:2)
      {
        cns<-unique(1:max(as.numeric(unlist(strsplit(gtypes,","))[-3])))
        totalcn<-sum(as.numeric(unlist(strsplit(gtypes[l],","))[-3]))
        prop<-as.numeric(unlist(strsplit(gtypes[l],","))[3])
        if(l==1)
        {
          nt<-sum(as.numeric(unlist(strsplit(gtypes[2],","))[-3]))
        }else{
          nt<-sum(as.numeric(unlist(strsplit(gtypes[1],","))[-3]))
        }
        for(k in 1:nrow(cluststruct))
        {
          llcn<-c()
          for(j in cns)
          {
            psi<-cluststruct[k,3]/purity
            if(psi>1){psi=1}
            llcn<-c(llcn,get_likely_subclonal(psi,purity,totalcn,j/totalcn,2,1-purity,curr$ref,curr$var,nt,prop))
          }
          llclust<-c(llclust,max(llcn,na.rm=T))
        }
      }
      dim(llclust)<-c(nrow(cluststruct),2)
      llclust<-apply(llclust,1,max)
    }else{
      #clonal copy-number
      cns<-unique(1:max(as.numeric(unlist(strsplit(gtypes,","))[-3])))
      totalcn<-sum(as.numeric(unlist(strsplit(gtypes,","))[-3]))
      llclust<-c()
      for(k in 1:nrow(cluststruct))
      {
        llcn<-c()
        for(j in cns)
        {
          psi<-cluststruct[k,3]/purity
          if(psi>1){psi=1}
          llcn<-c(llcn,get_likely_clonal(psi,purity,totalcn,j/totalcn,2,1-purity,curr$ref,curr$var))
        }
        llclust<-c(llclust,max(llcn,na.rm=T))
      }
    }

    if(all(pchisq(2*(llclust[-clonal]-llclust[clonal]),1,lower.tail=F)<pvalthresh))
    {
      workingdat[i,]$clust<-which.max(llclust)
    }else{
      workingdat[i,]$clust<-clonal
    }
  }
}

if(sum(workingdat$clust==0)>0)
{
cluststruct<-rbind(cluststruct,c(nrow(cluststruct)+1,0,0))
workingdat$clust[workingdat$clust==0]<-nrow(cluststruct)
}
for(i in 1:nrow(cluststruct))
{
  cluststruct[i,2]<-sum(workingdat$clust==i)
}

write.table(cluststruct,"./smc_1C_cluster_structure.txt",sep="\t",row.names=F,col.names = F)
write.table(workingdat$clust,"./smc_2A_mutations_to_clusters.txt",sep="\t",row.names=F,col.names = F)


#truth<-read.table("tumour1/Tumour1.truth.2A.txt")
#truthvcf<-read.table("tumour1/Tumour1.truth.scoring_vcf.vcf",sep="\t",header=F,stringsAsFactors = F)
#compare<-workingdat[truthvcf[,12]=="True",]
#sum(compare$clust==truth[,1])/length(compare$clust)
#h2h<-cbind(compare,truth[,1])
#head(h2h[h2h[,7]!=h2h[,8],])
