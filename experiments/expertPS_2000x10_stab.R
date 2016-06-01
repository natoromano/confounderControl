###############################################################################
# (Script)
# To simulated data (2000x10), perform expertPS to test various aspects 
# (stability, influence if IV, etc.), and compare with random matching.
#
# Nathanael Romano 
##############################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath="/home/yenlow/scripts"
source(paste(RScriptPath,"/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath,"/R/patientChar.R", sep="")) # patientchar

############# IMPORTS #############
installnewpackage(c("rJava", "Matching", "Epi", "glmnet", "randomForest", 
                    "vegan", "FNN", "Matrix", "doParallel", "foreach"))

require(rJava)
require(Matching)
require(Epi) # clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)
require(Matrix)

# parallel core computing for lassoMV bootstrapping
require(doParallel)
require(foreach)

# specify number of cores req
cl <- makeCluster(10) 
registerDoParallel(cl)

############# SET-UP #############
setwd("/home/naromano/projects/confounderControl/models/sim2000x10")
getwd()

baseDir = getwd();
tempDir = paste(baseDir, "tmp",sep="/") # output directory
patientfile=paste(baseDir, "patients.txt",sep="")
dir.create(tempDir)

# initialize the Java subsystem. include a jdbc object
.jinit(classpath="/home/yenlow/software/pharmacoepi.jar", 
       force.init = TRUE,parameters="-Xmx2g");
.jclassPath()

# required functions
source("../../hdpsFunctions.R")
source("../../rfFunctions.R")

id <- "id"
exposed <- "exposure" # name of exposure variable (must be string)
outcome <- "outcome"
expertvariables <- c("x1", "x2", "x3", "x4", "x5", "x6")  # exclude x7, IV
outcomemodelvariables_exclmatchingvar <- c("x8","x9","x10") 
# empirical variables may be matched depending on HD-PS
empvariables <- paste("x", 1:10, sep="")
empvariables_cat <- c("x1", "x3", "x5", "x6", "x8", "x9")
empvariables_num <- setdiff(empvariables, empvariables_cat)

# formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS <- as.formula(paste(exposed," ~ ",
                        paste(expertvariables,
                              collapse=" + "), sep=""))

# formula for outcome after expertPS
fmod <- paste("outcome ~ exposure +", 
              paste(outcomemodelvariables_exclmatchingvar,
                    collapse=" + ", sep=""))

desiredOrder <- c("n1", "ORmatched", "ORlow_matched", "ORupp_matched",
                  "coeff_matched", "se_matched", "bias_matched", "inCI_matched",
                  "n0", "ORadj", "ORlow_adj", "ORupp_adj", "coeff_adj",
                  "se_adj", "bias_adj", "inCI_adj", "SMD", "KS", "KL", "Cstat")

# to store consolidated results
startAll <- proc.time()[1:3]

Nsim <- 200 # number of simulations
Nmethods <- 2 # 2 methods (incl. 3 rdm sampling),
resultsArray <- array(dim=c(Nmethods, 20, Nsim))  
Noutcomesmat <- matrix(nrow=Nsim, ncol=2)
time <- array(dim=c(Nmethods, 3, Nsim))
dimnames(time)[[1]] <- c("expertPS", "rdmSamp1x")
dimnames(time)[[2]] <- c("user", "system", "elapsed")
dimnames(resultsArray)[[1]] <- dimnames(time)[[1]]

wantedVarList=matchedID=var_corOutcomeList=list()
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=c()


# foreach(i=1:Nsim) %dopar%{
for (i in 1:Nsim) {
  cat("\n\n\n-----------------------", i, "---------------------\n\n")
  
  # generate simulated data
  # source("../../data/2000by10.R")
  # load("../../data/2000x10.RData", verbose=T)
  # rm(list=ls()[!(ls() %in% "ds10")]) # remove all objects except ds10
  ds10$id=1:nrow(ds10)
  nminor <- min(table(ds10[,exposed]))
  
  trueORvec[i] <- trueOR
  trueSMDvec[i] <- trueSmd[2]
  trueCoeffvec[i] <- trueCoeff
  Nexposed[i] <- sum(ds10[,exposed])
  Noutcomesmat[i,] <- Noutcomes
  
  ############# EXPERT PSM #############
  # scale continuous variables (already scaled to standard normal)
  
  # calculate estimated PS
  print("############# expertPS ############")
  start <- proc.time()[1:3]
  
  #xmat <- as.matrix(ds10[,paste("x", 1:10, sep="")])
  #expertPsMod <- cv.glmnet(x=xmat[, expertvariables], y=ds10[, exposed],
  #                      family="binomial", standardize=F, alpha=0)
  #psl <- unlogit(as.numeric(predict(expertPsMod, xmat[, expertvariables], 
  #                                  s=expertPsMod$lambda.1se)))
  #names(psl) <- ds10[, id]
  #expertResults <- extractResults(ps=psl,
  #                                exposurevec=ds10[,exposed], fmod=fmod,
  #                                data=ds10, id=id, exposed=exposed, 
  #                                outcome=outcome, logitFlag=TRUE,
  #                                outfile=NULL, verbose=FALSE)
  
  expertPsMod <- glm(fPS, data=ds10[,c(id, exposed, expertvariables)], 
                     family="binomial")
  expertResults <- extractResults(ps=expertPsMod$fitted.values, 
                                  exposurevec=expertPsMod$y, fmod=fmod,
                                  data=ds10, id=id, exposed=exposed, 
                                  outcome=outcome, logitFlag=TRUE,
                                  outfile=NULL, verbose=FALSE)
  
  end <- proc.time()[1:3]
  time["expertPS",,i] <- end-start
  
  
  ####### random sampling #######
  print("############ random downsampling 1x ############")
  start <- proc.time()[1:3]
  
  matchedID_rdm <- cbind(ds10[ds10[,exposed]==1,id], 
                         sample(ds10[ds10[,exposed]==0,id], 
                                nminor, replace=FALSE))
  rdm1xResults <- extractResultsRdm(matchedID_rdm, ds10, fmod=NULL,
                                    id=id, exposed=exposed, outcome=outcome,
                                    verbose=FALSE)
  
  end <- proc.time()[1:3]
  time["rdmSamp1x",,i] <- end-start
 
  # consolidate matchedIDs
  matchedID[[i]] <- list(expertResults[[2]], rdm1xResults[[2]])
  names(matchedID[[i]]) <- dimnames(time)[[1]]
  
  # check if baseline variables are balanced
  beforeMatching <- matchStats(numvar=empvariables_num, catvar=empvariables_cat,
                               treatment=exposed, data=ds10,
                               outXlsx=NULL,verbose=FALSE)
  afterMatching <- list()
  for(j in 1:2){
      afterMatching[[j]] <- matchStats(numvar=empvariables_num,
                                       catvar=empvariables_cat,
                                       treatment=exposed,
                                       data=ds10[matchedID[[i]][[j]],],
                                       outXlsx=NULL,verbose=FALSE)
  }
  names(afterMatching) <- names(matchedID[[i]])[-(13:14)]
  
  # consolidate results
  resultsmat <- as.data.frame(rbind(expertResults[[1]], rdm1xResults[[1]]))
  rownames(resultsmat) <- dimnames(time)[[1]]
  resultsmat$bias_adj <- resultsmat$coeff_adj-trueCoeff
  resultsmat$bias_matched <- resultsmat$coeff_matched-trueCoeff
  resultsmat$inCI_adj <- (resultsmat$ORlow_adj<=trueOR &
                            resultsmat$ORupp_adj>=trueOR)
  resultsmat$inCI_matched <- (resultsmat$ORlow_matched<=trueOR & 
                                resultsmat$ORupp_matched>=trueOR)
  resultsArray[,,i] <- as.matrix(resultsmat[,desiredOrder])
  
  INSTAB_THR = 10 # instability threshold
  if (abs(resultsmat['expertPS','bias_matched']) > INSTAB_THR) {
    save(ds10, resultsmat, 
         file=paste(paste('instability', i, sep='_'), '.RData', sep=''))
  }
} # end of Nsim Monte Carlo simulations

endAll <- proc.time()[1:3]
endAll-startAll

# remember to terminate the cores when done
stopCluster(cl)

dimnames(resultsArray)[1:2] <- list(rownames(resultsmat), desiredOrder)
save(resultsArray, time, wantedVarList, matchedID, Noutcomesmat, Nexposed,
     trueORvec, trueSMDvec, trueCoeffvec, file="expertPS_2000x10_MS_stab.RData")

round(apply(resultsArray, c(1,2), mean, na.rm=T), 3)
apply(time[,"elapsed",], 1, summary)

rowSums(resultsArray[,"inCI_matched",], na.rm=T)  # excludes lassoMV
rowSums(resultsArray[,"inCI_adj",], na.rm=T)  # for lassoMV
# round(rowMeans(abs(resultsArray["expertPS",c("n1","KS","KL","Cstat","se_adj",
# "bias_adj","se_matched","bias_matched"),]),na.rm=T),3)
rowMeans(rbind(trueORvec, resultsArray[,"ORadj",]), na.rm=T)
rowMeans(rbind(trueORvec, resultsArray[,"ORmatched",]), na.rm=T)
