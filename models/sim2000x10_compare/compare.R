###############################################################################
# (Script)
# To sim data 2000 x 10 variables, do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman 
# correl
# 
# 29-Apr-14 Yen Low
##############################################################################

############# IMPORTS #############
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath <- "/home/yenlow/scripts"
source(paste(RScriptPath,"/R/utils.R", sep="")) # installnewpackage
ProjectPath <- "~/projects/confounderControl"
setwd(paste(ProjectPath, "/models/sim2000x10_compare", sep=""))
getwd()

installnewpackage(c('reshape2', 'lattice', 'vioplot', 'ggplot2', 'gridExtra'))
require(reshape2)
require(lattice)
require(vioplot)
require(ggplot2)
require(gridExtra)

load(file="../sim2000x10/sim2000x10_1000MC.RData")
source(file="../../charts.R")
apply(resultsArray, c(1,2), mean, na.rm=T)
results <- apply(abs(resultsArray[,c("bias_matched", "bias_adj"),]), 
                 c(1,2), mean, na.rm=T)

print(results)

############# PLOT RESULTS ARRAY #############
# extract results of interest -> put in list required of vioplot2log
# reverse the order from top to bottom (use 0 for lassoMV which is NA)
ORmatched=ORadj=bias_matched=bias_adj=se_matched=se_adj=list()
inCI_matched=inCI_adj=n1=SMD=Cstat=timelist=list() 
# (use 0 for lassoMV which is NA)

counter <- 0
for (j in 14:1) {  # excluse lassoMV (row11) which is NA
  counter <- counter+1
  ORmatched[[counter]] <- resultsArray[j,"ORmatched",]
  ORadj[[counter]] <- resultsArray[j,"ORadj",]
  bias_matched[[counter]] <- abs(resultsArray[j,"bias_matched",])
  bias_adj[[counter]] <- abs(resultsArray[j,"bias_adj",])
  se_matched[[counter]] <- resultsArray[j,"se_matched",]
  se_adj[[counter]] <- resultsArray[j,"se_adj",]
  n1[[counter]] <- resultsArray[j,"n1",]
  SMD[[counter]] <- resultsArray[j,"SMD",]
  Cstat[[counter]] <- resultsArray[j,"Cstat",]
  timelist[[counter]] <- time[j,3,]
}

# use 0 for lassoMV which is NA
ORmatched[[4]] <- 0   
bias_matched[[4]] <- 0
se_matched[[4]] <- 0
n1[[4]] <- 0

# use 0 for Euc and sim which are NA
for (j in c(1:3,5:10)) {
  ORadj[[j]] <- 0
  bias_adj[[j]] <- 0
  se_adj[[j]] <- 0 
}

for (j in 1:10) Cstat[[j]] <- 0.5  # No Cstat for non-PS methods
names(ORmatched)=rev(dimnames(resultsArray)[[1]])
names(ORadj)=rev(dimnames(resultsArray)[[1]])
names(bias_matched)=rev(dimnames(resultsArray)[[1]])
names(bias_adj)=rev(dimnames(resultsArray)[[1]])
names(se_matched)=rev(dimnames(resultsArray)[[1]])
names(se_adj)=rev(dimnames(resultsArray)[[1]])
names(n1)=rev(dimnames(resultsArray)[[1]])
names(SMD)=rev(dimnames(resultsArray)[[1]])
names(Cstat)=rev(dimnames(resultsArray)[[1]])
names(timelist)=rev(dimnames(resultsArray)[[1]])


############# VIOLIN PLOTS #############
# Function for repeated violin plots
vioplot2list <- function(list1, list2=NULL, xticks, title="", legendtext=NULL,
                         legendpos="bottomleft", log=TRUE, ref=1, outfile=NULL){
  colors=c("#CCCCCC", "#666666")
  if (!is.null(outfile)) png(file=outfile,width=5,height=5,units="in",res=150)
  if (!is.null(list2)) { # if 2 lists (assymmetrical violin)
  vioplot2log(list1,horizontal=TRUE,ylim=range(xticks),col=colors[1],
              xaxt="n",side="right",names=names(list1),mar=c(2,6,2,0.5),
              title=title,cex.main=1,boxsize=0.3,log=log)
  vioplot2log(list2,horizontal=TRUE,col=colors[2],add=T,log=log,
              xaxt="n",side="left",names=names(list2),boxsize=0.3,
              rectCol="white")
  } else {  #if 1 list (symmetrical violin)
    vioplot2log(list1,horizontal=TRUE,ylim=range(xticks),col=colors[1],
                xaxt="n",side="both",names=names(list1),mar=c(2,6,2,0.5),
                title=title,cex.main=1,boxsize=0.3,log=log)
  }
  axis(1, xticks, labels=xticks, las=1)
  abline(v=xticks,col="gray",lty="dotted")
  abline(v=ref,lwd=2)
  if(!is.null(legendtext)) legend(legendpos,legendtext,fill=colors,bty="n",
                                  cex=0.8)
  if(!is.null(outfile)) dev.off()
}

# violin plot for OR
vioplot2list(ORmatched,list2=ORadj,xticks=c(.1,1,5),log=TRUE,ref=1,
             legendpos="bottomright",
             title="OR (simulated data: 2000 rows x 10 variables)",
             legendtext=c("OR_matched","OR_adjusted"),
             outfile="OR_violin_sim2000x10.png")
vioplot2list(bias_matched,list2=bias_adj,xticks=seq(0,2,by=1),log=FALSE,ref=0,
             legendpos="bottomright",
             title="Bias (simulated data: 2000 rows x 10 variables)",
             legendtext=c("Bias_matched","Bias_adjusted"),
             outfile="bias_violin_sim2000x10.png")
vioplot2list(se_matched,list2=se_adj,xticks=seq(0,1,by=0.5),log=FALSE,ref=0,
             legendpos="bottomright",
             title="SE (simulated data: 2000 rows x 10 variables)",
             legendtext=c("SE_matched","SE_adjusted"),
             outfile="SE_violin_sim2000x10.png")
vioplot2list(n1,xticks=seq(1000,2000,by=200),log=FALSE,ref=2000,
             title="Matched cohort size (simulated data: 2000 rows x 10 variables)",
             legendtext=NULL,outfile="nmatched_violin_sim2000x10.png")
vioplot2list(SMD,xticks=seq(-0.3,0.1,by=0.1),log=FALSE,ref=1,
             title="Std Mean Difference (simulated data: 2000 rows x 10 variables)",
             legendtext=NULL,outfile="SMD_violin_sim2000x10.png")
vioplot2list(Cstat,xticks=seq(0.5,1,by=0.1),log=FALSE,ref=0.5,
             title="C-statistic of PS model (simulated data: 2000 rows x 10 variables)",
             legendtext=NULL,outfile="cstat_violin_sim2000x10.png")
vioplot2list(timelist,xticks=c(0.1,1,seq(10,50,by=10)),log=TRUE,ref=100,
             title="Computing time in sec (simulated data: 2000 rows x 10 variables)",
             legendtext=NULL,outfile="time_violin_sim2000x10.png")
# ORmatched
vioplot2list(ORmatched,list2=ORadj,xticks=c(.1,50),log=TRUE,ref=10,
             legendpos="bottomright",
             title="OR (simulated data: 2000 rows x 10 variables)",
             legendtext=c("OR_matched","OR_adjusted"),
             outfile="violinAnatomy.png")


# reverse order for barplot
reportmat <- apply(resultsArray[dim(resultsArray)[1]:1, 
                                c("inCI_adj","inCI_matched"),], c(1,2), sum,
                   na.rm=T)/dim(resultsArray)[3]

png(file="coverage_barplot_sim2000x10.png", width=5, height=5,
    units="in" ,res=150)
par(mar=c(2, 6.5, 2, 0.5))
barplot(t(reportmat), beside=T, horiz=TRUE,
        main="Coverage of true OR within 95% CI (sim data: 2000 x 10 var)",
        cex.main=1, las=1, xlim=c(0,1), axes=TRUE)
dev.off()
# abline(v=1, lwd=5)


############# LEGACY PLOTS #############
# dotplots with base R
metricsOfInterest <- dimnames(resultsArray)[[2]][grep("^OR", 
                                            dimnames(resultsArray)[[2]])]
data <- melt(resultsArray[,metricsOfInterest[1:3],])
colnames(data) <- c('Var1', 'Var2', 'Var3', 'value')
data$Var3 <- rep(11:1, nrow(data)/11)

png("dotplot_sim2000x10.png", width=5, height=5, units="in", res=150)
plot(Var3~value, subset(data, Var2=="ORlow_matched"),
     pch="<", cex=0.5, col="green", yaxt="n", ylab="",
     xlim=range(xticks), log="x") # main = title
axis(2, 11:1,labels=data$Var1[1:11], las=2)
points(Var3~value, subset(data, Var2=="ORupp_matched"),
       pch=">", cex=0.5, col="red")
points(Var3~value, subset(data, Var2=="ORmatched"), pch=16, cex=0.5)
abline(v=xticks, col="gray", lty="dotted")
abline(v=1, lwd=2)
dev.off()

# boxplots
getOption("na.action")
options(na.action="na.pass")
par(mar=c(4, 5, 2, 1))
boxplot(value~Var3,subset(data, Var2=="ORmatched"), horizontal=T,
        outline=T, names=data$Var1[11:1],las=1,log="x",xaxt="n")
xticks=c(0.001, 0.01, 0.1, 1, 5)
axis(1, xticks, labels=xticks, las=1)
abline(v=xticks, col="gray", lty="dotted")
abline(v=1, lwd=2)

# } ??


############# CALIPER TUNING #############
#### extract NsdEucCal results
load(file="../sim2000x10/tuneNsdEuc.RData")

selectedMet <- c("bias_matched", "SMD", "n1", "inCI_matched")
temp <- resultsArrayEuc[,selectedMet,]
means <- apply(temp, c(1,2), mean, na.rm=T)
sds <- apply(temp, c(1,2), sd, na.rm=T)
df <- as.data.frame(cbind(as.numeric(dimnames(resultsArrayEuc)[[1]]),
                          means, sds, means-sds, means+sds))
colnames(df) <- c("nsd",paste(selectedMet, ".mean", sep=""),
               paste(selectedMet,".sd",sep=""),
                     paste(selectedMet, ".low", sep=""),
               paste(selectedMet,".upp",sep=""))

p = ggplot(df[-1,],aes(x=nsd)) + theme_bw() + xlim(0,3) + ylab("Bias") 
                  + geom_hline(yintercept=0)
p = p + geom_ribbon(aes(ymin=bias_matched.low, ymax=bias_matched.upp),
                    fill="gray", alpha=0.3)
p = p + geom_line(aes(y=bias_matched.mean), size=1.5)
plot(p)

q = ggplot(df[-1,],aes(x=nsd)) + theme_bw() + ylim(500,1500) 
          + xlim(0,3) + ylab("Number of patients") + geom_hline(yintercept=0)
q = q + geom_ribbon(aes(ymin=n1.low, ymax=n1.upp), fill="gray", alpha=0.5)
q = q + geom_line(aes(y=n1.mean), size=1.5)
plot(q)

png("NsdEucCal.png", width=3.5, height=5, units="in", res=150)
grid.arrange(p, q, ncol=1)
dev.off()


############# EXTRACT SIMCALIPER RESULTS #############
load(file="../sim2000x10/tuneCalSim.RData")
# apply(resultsArrayCal[,"bias_matched",,],c(1,2),mean,na.rm=T)

temp <- resultsArrayCal[,selectedMet,,]
means <- apply(temp,c(1,2,3),mean,na.rm=T)
ses <- apply(temp,c(1,2,3),sd,na.rm=T)/dim(resultsArrayCal)[4]
