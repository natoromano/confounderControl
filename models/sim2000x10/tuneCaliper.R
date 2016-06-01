###############################################################################
# (Script)
# To simulated data 2000x10, perform similarity search by Euclidean dist, 
# Jaccard, Dice, Cosine similarities, Pearson and Spearman correllation.
# Tune caliper (Euc: mean+3sd, Jac: 0.6, others: 0.7 or 0.8)
# 
# 05-Apr-14 Yen Low
# 22-Feb-16 refactoring, Nathanael Romano
##############################################################################

############# SETUP #############
# Install and load packages
# install.packages(c('Matching', 'glmnet', 'randomForest', 'vegan', 'FNN'))
require(rJava)
require(Matching)
require(Epi) # clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)

setwd("~/Documents/Stanford/Shah Lab/confounderControl/models/sim2000x10")
# required functions
source("../../hdpsFunctions.R")

############# PARAMETERS #############
id <- "id"
exposed <- "exposure" # name of exposure variable (must be string)
outcome <- "outcome"

############# TUNE N. OF STANDARD DEVIATIONS #############
desiredOrder <- c("n1", "ORmatched", "ORlow_matched", "ORupp_matched",
                  "coeff_matched", "se_matched", "bias_matched", "inCI_matched",
                  "n0","ORadj", "ORlow_adj", "ORupp_adj", "coeff_adj", "se_adj",
                  "bias_adj", "inCI_adj", "SMD", "KS", "KL", "Cstat")
nsds <- rev(c(0, 1, 1.5, 2, 2.5, 3, 10, 1000))  # number of standard deviations

Nsim <- 100  # number of simulations
resultsArrayEuc <- array(dim=c(length(nsds), 10, Nsim))
Noutcomesmat <- matrix(nrow=Nsim,ncol=2)
trueORvec <- c() 
trueSMDvec <- c()
trueCoeffvec <- c()

for (i in 1:Nsim) {
  print('Simulation (euclidean) ' + i)
  # generate simulated data
  source("../../data/2000by10.R")
  ds10$id=1:nrow(ds10)
  
  trueORvec[i] <- trueOR
  trueSMDvec[i] <- trueSmd
  trueCoeffvec[i] <- trueCoeff
  Noutcomesmat[i,] <- Noutcomes
  
  xmat_ctrl <- xmat[ds10[, exposed]==0,]
  xmat_trted <- xmat[ds10[, exposed]==1,]
  rownames(xmat_ctrl) <- ds10[ds10[, exposed]==0, id]
  rownames(xmat_trted) <- ds10[ds10[, exposed]==1, id]

  resultsmat <- c()
  # for each caliper (defined by the number of standard deviations)
  for (q in 1:length(nsds)) {
    # match on euclidean distance
    matchedEuclidean <- matchByDist(xmat_ctrl, xmat_trted, method="euclidean",
                                    k_neighbors=5, caliper=0.8, nsd=nsds[q])
    # extract results
    euclideanResults <- extractResults(ps=matchedEuclidean, exposurevec=NULL,
                                    data=ds10, fmod=NULL, id=id,exposed=exposed,
                                    outcome=outcome,logitFlag=F, outfile=NULL,
                                    verbose=FALSE)
    
    # consolidate results
    resultsmat <- rbind(resultsmat, euclideanResults$results)
  }
  
    colnames(resultsmat) <- names(euclideanResults$results)
    resultsmat <- as.data.frame(resultsmat)
    resultsmat$bias_matched <- resultsmat$coeff_matched - trueCoeff
    resultsmat$inCI_matched <- (resultsmat$ORlow_matched<=trueOR &
                                  resultsmat$ORupp_matched>=trueOR)
    resultsArrayEuc[, , i] <- as.matrix(resultsmat[, desiredOrder[c(1:9, 17)]])  
}

dimnames(resultsArrayEuc)[1:2] <- list(nsds, desiredOrder[c(1:9,17)])
save(resultsArrayEuc, trueORvec, trueSMDvec, trueCoeffvec,
     file="tuneNsdEuc.RData")


# get final results matrix
reportmat <- 
rbind(
  # 1-3 sd (1sd is equiv to 1316 patients total)
  c(n, rowMeans(resultsArrayEuc[, "n1",])),
  rowMeans(rbind(trueORvec, resultsArrayEuc[, "ORmatched",])),  # 1-3 sd
  rowMeans(rbind(trueSMDvec, resultsArrayEuc[, "SMD",])),  # any
  c(0, rowMeans(abs(resultsArrayEuc[, "bias_matched",]))),  # 2-3 sd
  c(0,rowMeans(resultsArrayEuc[, "se_matched", ])),  # 1.5-3 sd
  c(1000,  rowSums(resultsArrayEuc[, "inCI_matched", ])) / 1000  # 2-3 sd
)
rownames(reportmat) <- c("n1", "ORmatched", "SMD", "bias_matched",
                         "se_matched", "inCI_matched")
colnames(reportmat)[1] <- "true_value"
write.table(reportmat, file="tuneNsdEuc.txt", sep="\t", col.names=NA)


# plot final results matrix
png(file="tuneNsdEuc.png", width=7, height=11, units="in",
    res=300, bg="transparent")
par(mfrow=c(3,2), mar=c(2,3,3,1))
barplot(reportmat[1,-1], space=0.1, ylim=c(0,2000), 
        main="Number of patients", cex.main=1)
abline(h=reportmat[1,1], lwd=5)
barplot(reportmat[2,-1], space=0.1,
        main="Odds ratio", cex.main=1)
abline(h=reportmat[2,1], lwd=5)
barplot(reportmat[3,-1], space=0.1,
        main="Effect size (standardized mean difference)", cex.main=1)
abline(h=reportmat[3,1], lwd=5)
barplot(reportmat[4,-1], space=0.1,
        main="Bias", cex.main=1)
abline(h=0, lwd=5)
barplot(reportmat[5,-1], space=0.1, ylim=c(0,0.4),
        main="Std error of beta", cex.main=1)
abline(h=0, lwd=5)
barplot(reportmat[6,-1], space=0.1, ylim=c(0,1),
        main="Coverage of true OR within 95% CI", cex.main=1)
abline(h=1, lwd=5)
dev.off()


############# TUNE CALIPER FOR SIMILARITY METHODS #############
calipers <- c(-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

Nsim <- 100  # number of simulations
resultsArrayCal <- array(dim=c(5, 10, length(calipers), Nsim))
Noutcomesmat <- matrix(nrow=Nsim, ncol=2)
trueORvec = trueSMDvec = trueCoeffvec = c()

for (i in 1:Nsim) {
  print('Simulation (similarities) ' + i)
  # generate simulated data
  source("../../data/2000by10.R")
  load("../../data/2000x10.RData", verbose=T)
  rm(list=ls()[!(ls() %in% "ds10")]) # remove all objects except ds10
  ds10$id <- 1:nrow(ds10)
  
  trueORvec[i] <- trueOR
  trueSMDvec[i] <- trueSmd
  trueCoeffvec[i] <- trueCoeff
  Noutcomesmat[i,] <- Noutcomes
  
  xmat_ctrl <- xmat[ds10[,exposed]==0,]
  xmat_trted <- xmat[ds10[,exposed]==1,]
  rownames(xmat_ctrl) <- ds10[ds10[,exposed]==0, id]
  rownames(xmat_trted) <- ds10[ds10[,exposed]==1, id]
  
  # for each caliper
  for(p in 1:length(calipers)){
    # compute matched sets
    matchedJaccard <- matchByDist(xmat_ctrl, xmat_trted, method="jaccard", 
                                  k_neighbors=5, caliper=calipers[p])
    matchedDice <- matchByDist(xmat_ctrl, xmat_trted, method="dice"
                               ,k_neighbors=5, caliper=calipers[p])
    matchedCosine <- matchByDist(xmat_ctrl, xmat_trted,method="cosine",
                                 k_neighbors=5, caliper=calipers[p])
    matchedPearson <- matchByDist(xmat_ctrl, xmat_trted,method="pearson",
                                  k_neighbors=5, caliper=calipers[p])
    matchedSpearman <- matchByDist(xmat_ctrl, xmat_trted,method="spearman",
                                   k_neighbors=5, caliper=calipers[p])
    # extract results
    jaccardResults <- extractResults(ps=matchedJaccard, exposurevec=NULL,
                                  data=ds10, fmod=NULL, id=id, exposed=exposed,
                                  outcome=outcome, logitFlag=F, outfile=NULL,
                                  verbose=FALSE)
    diceResults <- extractResults(ps=matchedDice, exposurevec=NULL, data=ds10,
                                  fmod=NULL, id=id, exposed=exposed, 
                                  outcome=outcome,logitFlag=F, outfile=NULL,
                                  verbose=FALSE)
    cosineResults <- extractResults(ps=matchedCosine, exposurevec=NULL,
                                 data=ds10, fmod=NULL, id=id, exposed=exposed,
                                 outcome=outcome, logitFlag=F, outfile=NULL,
                                 verbose=FALSE)
    pearsonResults <- extractResults(ps=matchedPearson, exposurevec=NULL,
                                  data=ds10, fmod=NULL, id=id,exposed=exposed,
                                  outcome=outcome, logitFlag=F, outfile=NULL,
                                  verbose=FALSE)
    spearmanResults <- extractResults(ps=matchedSpearman, exposurevec=NULL,
                                   data=ds10, fmod=NULL, id=id, exposed=exposed,
                                   outcome=outcome, logitFlag=F, outfile=NULL,
                                   verbose=FALSE)
    
    # consolidate results
    resultsmat <- as.data.frame(rbind(jaccardResults, diceResults, 
                                      cosineResults, pearsonResults,
                                      spearmanResults))
    rownames(resultsmat)=c("jaccard", "dice", "cosine", "pearson", "spearman")
    resultsmat$bias_matched <- sum(resultsmat$coeff_matched, -trueCoeff)
    resultsmat$inCI_matched <- (resultsmat$ORlow_matched<=trueOR &
                                  resultsmat$ORupp_matched>=trueOR)
    resultsArrayCal[,,p,i] <- as.matrix(resultsmat[, desiredOrder[c(1:9,17)]])
  }
  
  for (i in seq(10, 90, by=10)) {
    save(resultsArrayCal, trueORvec, trueSMDvec, trueCoeffvec, 
         file="tuneCalSim.RData")
  }
} # end of Nsim loop for similarities

dimnames(resultsArrayCal)[1:3] <- list(rownames(resultsmat),
                                       desiredOrder[c(1:9,17)], calipers)
save(resultsArrayCal, trueORvec, trueSMDvec, trueCoeffvec,
     file="tuneCalSim.RData")

dim(resultsArrayCal)
apply(resultsArrayCal, c(1,2,3), mean, na.rm=T)

sink("tuneCaliperSim.txt")
print("jaccard")
round(apply(resultsArrayCal["jaccard",,,], c(1,2), mean,na.rm=T), 3)  # 0.6
print("dice")
round(apply(resultsArrayCal["dice",,,], c(1,2), mean, na.rm=T), 3)  # 0.7/0.8
print("cosine")
round(apply(resultsArrayCal["cosine",,,], c(1,2), mean, na.rm=T), 3)  # 0.8/0.7
print("pearson")
round(apply(resultsArrayCal["pearson",,,], c(1,2), mean, na.rm=T), 3)  # 0.7/0.8
print("spearman")
round(apply(resultsArrayCal["spearman",,,], c(1,2), mean, na.rm=T), 3) # 0.7
sink()
