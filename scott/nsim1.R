################################################################################
# (Script)
# Run Scott Powers' simulation 1, but with extra methods.
# Methods:
# - naive estimation
# - naive lasso
# - direct estimation by machine learning (lassoMV)
# - PS-based: expertPS, lassoPS, hdPS, rfPS
# - similarity-based: euclidean
# - random downsampling
# For the new methods, we require at least 5 samples of each treatment class
#
# /!\ SHOULD BE RAN FROM COMMAND LINE, ELSE CALLING PYTHON SCRIPTS WITH SYSTEM
# DOES NOT WORK
#
################################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath="/home/yenlow/scripts"
source(paste(RScriptPath, "/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath, "/R/patientChar.R", sep="")) # patientchar
name <- 'simulation_1'

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
setwd('~/projects/confounderControl/scott/')

baseDir = getwd();
tempDir = paste(baseDir, "tmp", sep="/") # output directory
patientfile=paste(baseDir,"patients.txt",sep="")
dir.create(tempDir)

# initialize the Java subsystem. include a jdbc object
.jinit(classpath="/home/yenlow/software/pharmacoepi.jar", 
       force.init=TRUE, parameters="-Xmx2g");
.jclassPath()

# required functions
source("../hdpsFunctions.R")
source("../rfFunctions.R")

# treatment effect is 5 in F, 10 in M
# tr      m   f
# 0      10   10
# 1      15   20 

id <- "id"
exposed <- "exposure" # name of exposure variable (must be string)
outcome <- "outcome"
expertvariables <- c('gen', 'x1')
outcomemodelvariables_exclmatchingvar <- paste("x", 2:10, sep="")
# empirical variables may be matched depending on HD-PS
empvariables <- c('gen', paste("x", 1:10, sep=""))
empvariables_cat <- c('gen')
empvariables_num <- setdiff(empvariables, empvariables_cat)

# formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS <- as.formula(paste(exposed, " ~ ", paste(expertvariables,
                                              collapse=" + "), sep=""))

# formula for outcome after expertPS
fmod <- paste("outcome ~ exposure +", 
              paste(outcomemodelvariables_exclmatchingvar,
                    collapse=" + ", sep=""))

# desired metrics
desiredOrder <- c("n1", "effect_matched", "se_matched",
                "n0", "effect_adj", "se_adj", "SMD")

# simulation variables
set.seed(10)
n <- 1000
nsim <- 1

# initializations
Nmethods <- 10 # 10 methods
resultsArray <- array(dim=c(Nmethods, 7, nsim))
heterogeneousEffects <- array(dim=c(n, 3, nsim))
smdmatArray <- array(dim=c(Nmethods+1, length(empvariables), nsim))
Noutcomesmat <- matrix(nrow=nsim, ncol=2)
time <- array(dim=c(Nmethods, 3, nsim))
dimnames(time)[[1]] <- c("naive", "expertPS", "hdPS", "lassoPS", "rfPS", 
                         "powerForest", "propForest", "euclidean", "lassoMV", 
                         "rdmSamp1x")
dimnames(time)[[2]] <- c("user", "system", "elapsed")
dimnames(resultsArray)[[1]] <- dimnames(time)[[1]]

# objects to handle various methods's intermediary results
rfImpt <- matrix(nrow=length(empvariables), ncol=nsim)
rownames(rfImpt) <- empvariables
lassoPSBeta <- rfImpt
lassoMVvar=wantedVarList=matchedID=var_corOutcomeList=rdmBag=rdmJK=list()
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=bestlambda=c()
estimators <- 50 # number of estimators for causal forests

# genre variable
gen <- c(rep(0,n/2), rep(1,n/2))

createPSResultsArray <- function(extracted) {
  results <- rep(NA, 7)
  names(results) <- desiredOrder
  results['n0'] <- extracted[[1]]['n0']
  results['n1'] <- extracted[[1]]['n1']
  results['effect_adj'] <- extracted[[1]]['coeff_adj']
  results['effect_matched'] <- extracted[[1]]['coeff_matched']
  results['se_adj'] <- extracted[[1]]['se_adj']
  results['se_matched'] <- extracted[[1]]['se_matched']
  results['SMD'] <- extracted[[1]]['SMD']
  return(list(results=results, matchedID=extracted[[2]]))
}


for(i in 1:nsim){
  cat('\n\n------', i, '-----\n')
  # generate PS
  ps <-rep(NA, n)
  ps[gen == 0] <- .25
  ps[gen == 1] <- .75
  
  # generate treatment
  tr <- rep(NA, n)
  tr[gen == 0] <- 1*(runif(n/2) < .25)
  tr[gen == 1] <- 1*(runif(n/2) <.75)
  
  # generate noise variabless
  x <- matrix(rnorm(n*10), ncol = 10)
  bnoise1 <- c(5, 5, rep(0, 8))
  bnoise2 <- c(0, 0, -3, 2, rep(0, 6))
  
  # generate treatment effect (based on gender)
  y <- extra <- rep(NA, n)
  tre <- 5*(gen == 0) + 10*(gen == 1)
  extra[gen == 0] <- (x[gen == 0, ]%*%bnoise1)
  extra[gen == 1] <- (x[gen == 1, ]%*%bnoise2)
  # generate outcome
  y <- 10 + tr*tre + extra + rnorm(n)
  
  # merge dataset
  ds <- as.data.frame(cbind(y, ps, tr, gen, x))
  ds[, id] <- 1:nrow(ds)
  rownames(ds) <- 1:nrow(ds)
  colnames(ds) <- c(outcome, 'ps', exposed, empvariables, id)
  Nexposed[i] <- sum(ds[, exposed])
  Noutcomesmat[i,] <- sum(ds[, outcome])
  nminor <- min(table(ds[, exposed]))
  
  ############# naive estimation #############
  print("############# naive estimation ############")
  start <- proc.time()[1:3]
  
  results <- rep(NA, 7)
  names(results) <- desiredOrder
  results['n0'] <- nrow(ds)
  results['n1'] <- nrow(ds)
  results['effect_adj'] <- mean(y[tr == 1]) - mean(y[tr == 0])
  naiveResults <- list(results=results, matchedID=NA)
  end <- proc.time()[1:3]
  time["naive",,i] <- end-start
  
  
  ############# expertPS #############
  print("############ expertPS ############")
  start=proc.time()[1:3]
  xmat <- as.matrix(ds[, empvariables])
  expertPsMod <- cv.glmnet(x=xmat[, expertvariables], y=ds[, exposed],
                           family="binomial", standardize=F, alpha=1, nfold=5)
  psl <- unlogit(as.numeric(predict(expertPsMod, xmat[, expertvariables], s=0)))
  names(psl) <- ds[, id]
  tryobjN <- try(expertResults <- extractResults(ps=psl,
                                                 exposurevec=ds[, exposed], 
                                                 fmod=fmod, data=ds, id=id, 
                                                 exposed=exposed, 
                                                 outcome=outcome, 
                                                 logitFlag=TRUE,
                                                 outfile=NULL, verbose=FALSE,
                                                 continuous=TRUE))
  expertResults <- createPSResultsArray(expertResults)
  end=proc.time()[1:3]
  time["expertPS",,i] <- end-start
  
  if (FALSE) {
  ############# hdPS #############
  print("############# hdPS ############")
  start <- proc.time()[1:3]
  
  # prepare data into format required by hdps
  # read data file as single string (to use addPatientsFromBuffer)
  patientheader <- paste(id, exposed, outcome, sep="\t")
  patientstring <- paste(paste(ds[, id], ds[, exposed], ds[,outcome],
                               sep="\t"), collapse="\n")
  # for reading in whole file as single string
  # patientstring=readChar(patientfile, file.info(patientfile)$size)
  # patientstring=gsub("\"","",patientstring)
  # patients=dataset[,c("pid_org",exposed,outcome[1])] #get patient id, exp, out
  # write.table(patients,file=patientfile,sep="\t",
  # col.names=T,row.names=F,na="") # output data to patient.txt;
  datainstring <- paste(patientheader, patientstring, sep="\n")
  variables <- empvariables_cat
  dimdata <- ds[, c(id, variables)]
  
  ### INVOKE pharmacoepi.jar ### (handles cat variables only)
  hdpsobj <- hdps(datainstring, dimdata, outDir=tempDir, Nmostfreq=10, k=10,
                  stratifyDim=FALSE, outfile="output_cohort.txt",
                  FullOutput=TRUE, verbose=T, ZeroCellCorrection=F)
  hdpsobj$selectedvariables # (see Fig2, Schneeweiss 2008)
  wantedvar <- hdpsobj$selectedvariables[grep("1Once$",
                                              hdpsobj$selectedvariables)]
  
  var_corExposed <- empvariables_num[abs(cor(ds[, exposed],
                                             ds[, empvariables_num]))>0.05]
  var_corOutcome <- empvariables_num[abs(cor(ds[, outcome], # include in PS
                                             ds[, empvariables_num]))>0.05]
  IV <- setdiff(var_corExposed, var_corOutcome) # exclude from PS
  
  # Estimate the PS (force numerical variables into PS model)
  dataPS <- cbind(ds[,c(id, exposed, var_corOutcome)],
                  hdpsobj$hdpsdata[,wantedvar])
  hdPsMod <- glm(paste(exposed, "~ . -", id), data=dataPS, family="binomial")
  names(hdPsMod$fitted.values) <- as.character(dataPS[,id])
  hdResults <- extractResults(ps=hdPsMod$fitted, exposurevec=hdPsMod$y,
                              fmod=NULL, data=ds, id=id, exposed=exposed,
                              outcome=outcome, logitFlag=TRUE, outfile=NULL,
                              verbose=FALSE, continuous=TRUE)
  hdResults <- createPSResultsArray(hdResults)
  end <- proc.time()[1:3]
  time["hdPS",,i] <- end-start
  }
  time["hdPS",,i] <- 0
  results = rep(NA, 7)
  names(results) <- desiredOrder
  wantedvar <- NA
  var_corExposed <- NA
  var_corOutcome <- NA
  IV <- NA
  hdResults <- list(results=results, matchedID=NA)
  
  ############### lassoPS ############
  print("############# lassoPS ############")
  start=proc.time()[1:3]
  
  xmat <- as.matrix(ds[, empvariables])
  lassoPsmod <- cv.glmnet(xmat, ds[,exposed], alpha=1, family="binomial",
                          standardize=F, nfold=5)
  bestlambda[i] <- lassoPsmod$lambda.1se
  lassoPSBeta[,i] <- coeffAtlambda(lassoPsmod)[-1]  # exclude intercept
  # get estimated ps in logit form
  psl <- unlogit(as.numeric(predict(lassoPsmod, xmat, s=lassoPsmod$lambda.1se)))
  names(psl) <- ds[,id]
  lassoResults <- extractResults(ps=psl, exposurevec=ds[,exposed], fmod=NULL,
                                 data=ds, id=id, exposed=exposed,
                                 outcome=outcome, logitFlag=TRUE,
                                 outfile=NULL, verbose=FALSE, continuous=TRUE)
  lassoResults <- createPSResultsArray(lassoResults)
  end <- proc.time()[1:3]
  time["lassoPS",,i] <- end-start
  
  
  ############# randomForest #############
  print("############# rfPS ############")
  start <- proc.time()[1:3]
  rfPsMod <- randomForest(xmat, ds[,exposed], ntree=100, importance=T,
                          nodesize=100)
  ps <- rfPsMod$predicted
  names(ps) <- ds[,id]
  rfImpt[,i] <- rfPsMod$importance[,1]
  tryobj <- try(extractResults(ps=ps, exposurevec=ds[,exposed], data=ds,
                               fmod=NULL, id=id,exposed=exposed,
                               outcome=outcome, logitFlag=TRUE,
                               outfile=NULL, verbose=FALSE, continuous=TRUE))
  if(class(tryobj)!="try-error") rfResults<-tryobj else rfResults<-NA
  rfResults <- createPSResultsArray(rfResults)
  end <- proc.time()[1:3]
  time["rfPS",,i] <- end-start
  
  
  ############# Causal Forests #############
  print("############# power forest ############")
  start <- proc.time()[1:3]
  powerEffects <- causalForest(ds, method="power", covariates=empvariables,
                               name=name, estimators=estimators, desired='effect',
                               exposure=exposed, outcome=outcome)
  results <- rep(NA, 7)
  names(results) <- desiredOrder
  results['n0'] <- nrow(ds)
  results['n1'] <- nrow(ds)
  results['coeff_adj'] <- mean(powerEffects, na.rm=T)
  powerResults <- list(results=results, matchedID=NA)
  end <- proc.time()[1:3]
  time["powerForest",,i] <- end-start
  
  
  print("############ propensity forest ############")
  start <- proc.time()[1:3]
  propEffects <- causalForest(ds, covariates=empvariables, method="propensity",
                              name=name, estimators=estimators, desired='effect',
                              exposure=exposed, outcome=outcome)
  results <- rep(NA, 7)
  names(results) <- desiredOrder
  results['n0'] <- nrow(ds)
  results['n1'] <- nrow(ds)
  results['coeff_adj'] <- mean(propEffects, na.rm=T)
  propResults <- list(results=results, matchedID=NA)
  end <- proc.time()[1:3]
  time["propForest",,i] <- end-start 
  
  # reset current directory
  setwd('~/projects/confounderControl/scott/')
  
  
  ############# similarity #############
  print("############ similarity ############")
  
  xmat_ctrl <- xmat[ds[, exposed]==0,]
  xmat_trted <- xmat[ds[, exposed]==1,]
  rownames(xmat_ctrl) <- ds[ds[, exposed]==0, id]
  rownames(xmat_trted) <- ds[ds[, exposed]==1, id]
  
  runSimPS <- function(method="euclidean", caliper=0.7, nsd=3,
                       algorithm="brute") {
    matchedSim <- matchByDist(xmat_ctrl, xmat_trted, method=method,
                              k_neighbors=5, caliper=caliper, nsd=nsd,
                              algorithm=algorithm)
    simResults <- extractResults(ps=matchedSim, exposurevec=NULL,
                                 data=ds, fmod=NULL, id=id, exposed=exposed,
                                 outcome=outcome, logitFlag=F,
                                 outfile=NULL, verbose=FALSE, continuous=TRUE)
    return(simResults)
  }
  
  start <- proc.time()[1:3]
  euclideanResults <- runSimPS(method="euclidean", nsd=3, algorithm="brute")
  euclideanResults <- createPSResultsArray(euclideanResults)
  end=proc.time()[1:3]
  time["euclidean",,i] <- end-start
  
  
  print("############ lassoMV model ############")
  start <- proc.time()[1:3]
  
  xmat <- Matrix(as.matrix(ds[, c(exposed, empvariables)]), sparse=T)
  penalty.factor <- c(0, 0, rep(1, ncol(xmat)-2)) # force W and gen into model
  # tune lambda for glmnet using 5-fold CV
  lassoMvmod <- cv.glmnet(xmat, ds[,outcome], alpha=1, standardize=F, nfold=5,
                          penalty.factor=penalty.factor, family="gaussian")
  coeff <- na.omit(sort(coeffAtlambda(lassoMvmod)))
  
  # bootstrap, repeat 100 times for 95%CI
  # use global lambda (lambda.1se from lassomod.cv)
  lassomodbootstrap <- genCI(xmat, ds[,outcome], ntrials=100,
                             lambda=lassoMvmod$lambda.1se, alpha=1,
                             penalty.factor=penalty.factor, family="gaussian")
  lassomodCI <- setCL(lassomodbootstrap)
  lassoMVvar[[i]] <- lassomodCI[[1]]  # get coeff
  
  # extract results
  ORCInonZero=as.data.frame(exp(lassomodCI$betaCIlim[lassomodCI$beta_nonZero,]))
  Smd <- smd(ds[,c(exposed,outcome)], exposed=exposed, variable=outcome,
             verbose=FALSE, categorical=TRUE)
  se_adj <- (lassomodCI$betaCIlim[exposed,3] - 
               lassomodCI$betaCIlim[exposed,1])/2/1.96
  results <- rep(NA, 7)
  names(results) <- desiredOrder
  results['n0'] <- nrow(ds)
  results['n1'] <- nrow(ds)
  results['SMD'] <- Smd
  results['effect_adj'] <- lassomodCI$betaCIlim[exposed,2]
  results['se_adj'] <- se_adj
  lassoMVResults <- list(results=results, matchedID=NULL)
  
  end=proc.time()[1:3]
  time["lassoMV",,i] <- end-start
  
  
  ####### random sampling #######
  print("############ random downsampling 1x ############")
  start <- proc.time()[1:3]
  matchedID_rdm <- cbind(ds[ds[,exposed]==1,id], 
                         sample(ds[ds[,exposed]==0,id], 
                                nminor, replace=FALSE))
  rdm1xResults <- extractResultsRdm(matchedID_rdm, ds, fmod=NULL,
                                    id=id, exposed=exposed, outcome=outcome,
                                    verbose=FALSE, continuous=TRUE)
  rdm1xResults <- createPSResultsArray(rdm1xResults)
  end <- proc.time()[1:3]
  time["rdmSamp1x",,i] <- end-start
  
  # consolidate matchedIDs
  matchedID[[i]] <- list(naiveResults, 
                         expertResults[[2]], hdResults[[2]], lassoResults[[2]],
                         rfResults[[2]], powerResults[[2]], 
                         propResults[[2]], euclideanResults[[2]],
                         lassoMVResults[[2]], rdm1xResults[[2]])
  names(matchedID[[i]]) <- dimnames(time)[[1]]
  
  # consolidate results
  resultsmat <- as.data.frame(rbind(naiveResults[[1]], expertResults[[1]], 
                                    hdResults[[1]], lassoResults[[1]], 
                                    rfResults[[1]], powerResults[[1]], 
                                    propResults[[1]], euclideanResults[[1]], 
                                    lassoMVResults[[1]], rdm1xResults[[1]]))
  rownames(resultsmat) <- dimnames(time)[[1]]
  print(dim(resultsArray))
  print(dim(as.matrix(resultsmat[, desiredOrder])))
  save(resultsmat, resultsArray, file="test.RData")
  resultsArray[,,i] <- as.matrix(resultsmat[, desiredOrder])
  wantedVarList[[i]] <- wantedvar
  var_corOutcomeList[[i]] <- var_corOutcome
  heterogeneousEffects[,,i] <- as.matrix(data.frame(effect=tre,
                                          power=powerEffects,
                                          propensity=propEffects))

  # save periodidically
  if (i %in% seq(50, 900, by=50)) {
    save(ds, resultsArray, time, wantedVarList, rfImpt, lassoPSBeta, lassoMVvar,
         bestlambda, var_corOutcomeList, matchedID, Noutcomesmat, Nexposed,
         heterogeneousEffects, file=paste(name, 'results.RData', sep="_"))
  }
  
}

# remember to terminate the cores when done
stopCluster(cl)

dimnames(resultsArray)[1:2] <- list(rownames(resultsmat), desiredOrder)
dimnames(heterogeneousEffects)[1:2] <- list(1:n, 
                                            c('effect', 'power', 'propensity'))
save(resultsArray, time, wantedVarList, rfImpt, lassoPSBeta, lassoMVvar,
     bestlambda, var_corOutcomeList, matchedID, Noutcomesmat, Nexposed,
     heterogeneousEffects, file=paste(name, 'results.RData', sep="_"))


## DISREGARD
if (FALSE) {
  # compute absolute deviations from true treatment effect
  err.Naive <- c(matrix(est.Naive, nsim, n) - tre)
  err.NaivePS <- c(est.NaivePS - tre)
  err.Lasso <- c(matrix(est.Lasso, nsim, n) - tre)
  err.LassoPS <- c(est.LassoPS - tre)
  err.PropensityTree <- c(est.propensity.tree - tre)
  err.PropensityForest <- c(est.propensity.forest - tre)
  err.PowersTree <- c(est.powers.tree - tre)
  err.PowersForest <- c(est.powers.forest - tre)
  
  calc.bias <- function(x) {
    return(c(mean(abs(x)), sqrt(var(abs(x))/length(x))))
  }
  
  out <- rbind(calc.bias(err.Naive),
               calc.bias(err.NaivePS),
               calc.bias(err.Lasso),
               calc.bias(err.LassoPS),
               calc.bias(err.PropensityTree),
               calc.bias(err.PropensityForest),
               calc.bias(err.PowersTree),
               calc.bias(err.PowersForest))
  
  dimnames(out) <- list(c('err.Naive',
                          'err.NaivePS',
                          'err.Lasso',
                          'err.LassoPS',
                          'err.PropensityTree',
                          'err.PropensityForest',
                          'err.PowersTree',
                          'err.PowersForest'),
                        c("mean","se"))
  
  # mses
  errs <- c(mean(err.Naive^2), mean(err.NaivePS^2), mean(err.Lasso^2),
            mean(err.LassoPS^2), mean(err.PropensityTree^2),
            mean(err.PropensityForest^2),
            mean(err.PowersTree^2),
            mean(err.PowersForest^2))
  
  save.image('data/sim1-image-4-5-16.rda')
  
  q(save='no')
}
