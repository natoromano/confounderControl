# To all 2 simulated data sets, compare all 14 methods on all statistics
# 
# 23-Oct-14 fix IV to just x7
# 15-Oct-14 use rfPS in reg mode (set min node size to 5% of sample size, i.e. 100, Malley 2012) 
# 11-Sep-14 moved rdmSamp1x to last row
# 10-Sep-14 Yen Low
##############################################################################
#
#
#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
# RScriptPath="~/scripts"
# ProjectPath="~/projects"
setwd("../models/compareAll")
getwd()


###function for extraction results
getResults<-function(resultsArray,time,fn="mean",bias=TRUE){
  summaryTab=cbind(apply(resultsArray,c(1,2),get(fn),na.rm=T),
                   time=apply(time[,3,],1,get(fn),na.rm=T))
  
  #bias need special handling cos must take abs before meam
  if(bias==TRUE){
    biasCols=apply(abs(resultsArray[,c("bias_matched","inCI_matched","bias_adj","inCI_adj"),]),c(1,2),get(fn),na.rm=T)
    summaryTab[,c("bias_matched","inCI_matched","bias_adj","inCI_adj")]=biasCols     #replace bias
    trueOR=summary(trueORvec,na.rm=T)
    trueBeta=summary(trueCoeffvec,na.rm=T)
    list(summaryTab=summaryTab,trueOR=trueOR,trueBeta=trueBeta)
  }else return(summaryTab)
}

#################
modelOrder=c("expertPS","hdPS","lassoPS","rfPS",
             "euclidean","jaccard","dice","cosine","pearson","spearman",
             "lassoMV","rdmBag","rdmJK","rdmSamp1x")

####load resultsArray and time
load(file="../sim2000x10/sim2000x10_1000MC.RData", verbose=T)
#overwrite resultsArray with trimmed version
load(file="../sim2000x10/resultsArray_trimmed.RData", verbose=T)
mean1=getResults(resultsArray,time,"mean",bias=TRUE)
sd1=getResults(resultsArray,time,"sd",bias=TRUE)
mean1$summaryTab=mean1$summaryTab[modelOrder,]
sd1$summaryTab=sd1$summaryTab[modelOrder,]

load(file="../sim2000x100/sim2000x100_1000MC.RData", verbose=T)
mean2=getResults(resultsArray,time,"mean",bias=TRUE)
sd2=getResults(resultsArray,time,"sd",bias=TRUE)
mean2$summaryTab=mean2$summaryTab[modelOrder,]
sd2$summaryTab=sd2$summaryTab[modelOrder,]

wantedVar=c("bias_matched","ORmatched","se_matched","inCI_matched","n1",
            "bias_adj","ORadj","se_adj","inCI_adj",
            "Cstat","KS","time")

#append rownames with suffix
rownames(mean1$summaryTab)=paste(rownames(mean1$summaryTab),"1",sep=".")
rownames(mean2$summaryTab)=paste(rownames(mean2$summaryTab),"2",sep=".")

#merge means from 2 data sets
temp=rbind(mean1$summaryTab,mean2$summaryTab)
results_mean=as.data.frame(round(temp[,wantedVar],2))
            
#repeat for SD
#append rownames with suffix
rownames(sd1$summaryTab)=paste(rownames(sd1$summaryTab),"1",sep=".")
rownames(sd2$summaryTab)=paste(rownames(sd2$summaryTab),"2",sep=".")

#merge means from 2 data sets
temp=rbind(sd1$summaryTab,sd2$summaryTab)
results_sd=as.data.frame(round(temp[,wantedVar],2))

#put in printable format (mean with sd in brackets)
results=matrix(NA,nrow=nrow(results_mean),ncol=ncol(results_mean))
for(i in 1:nrow(results)) for(j in 1:ncol(results)) results[i,j]=paste(results_mean[i,j]," (",results_sd[i,j],")",sep="")
dimnames(results)=dimnames(results_mean)
results=gsub(".*NA.*","-",results)
results=gsub("\\(0\\)$","\\(0.00\\)",results)
results=gsub("^1 \\(","1.00 \\(",results)
results[,"n1"]=gsub("\\.[0-9]+","",results[,"n1"])
write.table(results,"summaryTab.txt",sep="\t",quote=F,col.names=NA)

