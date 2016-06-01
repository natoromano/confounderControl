################################################################################
# (Script)
# Run a simulation based on covariates extracted from clinical data, and 
# performs treatment effect estimation methods on it.
# Methods:
# - naive estimation
# - naive lasso
# - direct estimation by machine learning (lassoMV)
# - PS-based: expertPS, lassoPS
# - similarity-based: euclidean
# - causal forests
# - random downsampling
# 
#
# /!\ SHOULD BE RAN FROM COMMAND LINE, ELSE CALLING PYTHON SCRIPTS WITH SYSTEM
# DOES NOT WORK
#
# For the new methods, we require at least 5 samples of each treatment class
################################################################################
rm(list=ls())
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath="/home/yenlow/scripts"
source(paste(RScriptPath, "/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath, "/R/patientChar.R", sep="")) # patientchar

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
setwd('~/projects/confounderControl/realistic/')

baseDir = getwd();
tempDir <- paste(baseDir, "tmp", sep="/") # output directory
patientfile <- paste(baseDir, "patients.txt", sep="")
dir.create(tempDir)

# initialize the Java subsystem. include a jdbc object
.jinit(classpath="/home/yenlow/software/pharmacoepi.jar", 
       force.init=TRUE, parameters="-Xmx2g");
.jclassPath()

# required functions
source("../hdpsFunctions.R")
source("../rfFunctions.R")

id <- "id"
exposed <- "exposure" # name of exposure variable (must be string)
outcome <- "outcome"

# desired metrics
desiredOrder <- c("n1", "effect_matched", "bias_matched", "mse_matched",
                  "se_matched", "n0", "effect_adj", "bias_adj", "mse_adj",
                  "se_adj", "SMD")
preOrder <- c('n0', 'n1', 'effect_adj', 'effect_matched', 'se_adj', 
              'se_matched', 'SMD')

# simulation variables
set.seed(10)
n <- 1000
if (!exists('name')) name <- 'sim_tests'
if (!exists('nsim')) nsim = 2
if (!exists('n_EHR_var')) n_EHR_var = 10
if (!exists('alp')) alp = 2

# initializations
Nmethods <- 8
resultsArray <- array(NA, dim=c(Nmethods, 11, nsim))
heterogeneousEffects <- array(dim=c(n, 5, nsim))
Noutcomesmat <- matrix(nrow=nsim, ncol=2)
time <- array(dim=c(Nmethods, 3, nsim))
dimnames(time)[[1]] <- c("naive", "expertPS", "lassoPS", "powersForest", 
                         "propForest", "euclidean", "lassoMV", 
                         "rdmSamp1x")
dimnames(time)[[2]] <- c("user", "system", "elapsed")

# objects to handle various methods's intermediary results
lassoMVvar=matchedID=list()
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=bestlambda=c()
estimators <- 800 # number of estimators for causal forests

createPSResultsArray <- function(extracted) {
  results <- c(n0=extracted[[1]]['n0'], n1=extracted[[1]]['n1'],
               effect_adj=extracted[[1]]['coeff_adj'],
               effect_matched=extracted[[1]]['coeff_matched'],
               se_adj=extracted[[1]]['se_adj'],
               se_matched=extracted[[1]]['se_matched'],
               SMD=extracted[[1]]['SMD'])
  return(list(results=results, matchedID=extracted[[2]]))
}

# load EHR data
load('EHR_source.RData')

for(i in 1:nsim){
  cat('\n\n------', i, '-----\n')
  
  #### SIMULATE DATASET ###
  # generate genre variable
  gen <- c(rep(0, n/2), rep(1, n/2))
  
  # pool from EHR data
  patients <- sample(1:nrow(d), n)
  x_ehr <- d[patients, sample(1:ncol(d), n_EHR_var)]
  x_ehr <- matrix(rnorm(n * n_EHR_var, 0, 1), n, n_EHR_var)
  # scale necessary variables
  x_ehr[, colMeans(x_ehr, na.rm=T) > 1] <- scale(x_ehr[, colMeans(x_ehr, 
                                                                  na.rm=T) > 1])
  rownames(x_ehr) <- 1:nrow(x_ehr)
  age <- scale(d[patients, 'age'])
  
  # merge together
  covariates <- cbind(gen, age, x_ehr)
  colnames(covariates) <- c('gen', 'age', paste('x', 1:ncol(x_ehr), sep=""))
  
  # set empirical variables names
  empvariables <- colnames(covariates)
  
  # generate PS
  ps_var <- c('gen', 'age', paste('x', sample(1:ncol(x_ehr), 3), sep=""))
  b_ps <- c(0, alp * - 1.5, alp * 1.43, -0.55, 0.22, 0.21)
  ps <- c(1/(1 + exp(-(b_ps[1] + b_ps[-1] %*% t(covariates[, ps_var])))))
  
  # generate treatment
  tr <- rep(NA, n)
  tr <- 1 * (runif(n) < ps)
  
  # generate treatment effect
  # b_tre <- rep(0, ncol(covariates))
  # names(b_tre) <- colnames(covariates)
  # hetero_var <- paste('x', sample(1:ncol(x_ehr), 3), sep='')
  # b_tre[hetero_var] <- rnorm(3, 0, 10)
  tre <- 10 # + as.matrix(covariates) %*% b_tre
  
  # set relevant expert variables
  expertvariables <- ps_var
  
  # generate outcome
  b_risk <- rep(0, ncol(covariates))
  names(b_risk) <- colnames(covariates)
  risk_factors <- setdiff(paste('x', sample(1:ncol(x_ehr), 3), sep=''), ps_var)
  b_risk[risk_factors] <- rnorm(length(risk_factors), 0, 10)
  y <- -4.5 + alp * 1.2 * gen - alp * 0.3 * age + 
    tr*tre + as.matrix(covariates) %*% b_risk + rnorm(n)
  
  # consolidate data
  ds <- as.data.frame(cbind(y, ps, tr, covariates))
  ds[, id] <- 1:nrow(ds)
  rownames(ds) <- 1:nrow(ds)
  colnames(ds) <- c(outcome, 'ps', exposed, empvariables, id)
  Nexposed[i] <- sum(ds[, exposed])
  Noutcomesmat[i,] <- sum(ds[, outcome])
  nminor <- min(table(ds[, exposed]))
  
  # set formulas for linear models
  # formula for expertPS model: exposure ~ var1 + var 2 + ...
  expertvariables = setdiff(expertvariables, c('age'))
  fPS <- as.formula(paste(exposed, " ~ ", paste(expertvariables,
                                                collapse=" + "), sep=""))
  
  # formula for outcome after expertPS
  if (length(risk_factors) > 0) {
    fmod <- paste("outcome ~ exposure +", paste(risk_factors,
                                                collapse=" + ", sep=""))
  } else {
    fmod <- "outcome ~ exposure"
  }
  
  
  ############# naive estimation #############
  print("############# naive estimation ############")
  start <- proc.time()[1:3]
  results <- c(n0=nrow(ds), n1=nrow(ds), 
               effect_adj=(mean(y[tr == 1]) - mean(y[tr == 0])),
               effect_matched=NA, se_adj=NA, se_matched=NA, SMD=NA)
  naiveResults <- list(results=results[preOrder], matchedID=NA)
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
  tryobje <- try(expertResults <- extractResults(ps=psl, 
                                                 exposurevec=ds[, exposed], 
                                                 fmod=fmod, data=ds, id=id, 
                                                 exposed=exposed, 
                                                 outcome=outcome, 
                                                 logitFlag=TRUE, outfile=NULL, 
                                                 verbose=FALSE, 
                                                 continuous=TRUE))
  if (class(tryobje) != 'try-error') {
    expertResults <- createPSResultsArray(expertResults)
  } else {
    expertResults <- c(n0=NA, n1=NA, effect_adj=NA, effect_matched=NA,
                       se_adj=NA, se_matched=NA, SMD=NA)
  }  
  end <- proc.time()[1:3]
  time["expertPS",,i] <- end-start
  
  
  ############### lassoPS ############
  print("############# lassoPS ############")
  start=proc.time()[1:3]
  
  xmat <- as.matrix(ds[, empvariables])
  lassoPsmod <- cv.glmnet(xmat, ds[, exposed], alpha=1, family="binomial",
                          standardize=F, nfold=5)
  bestlambda[i] <- lassoPsmod$lambda.1se
  # get estimated ps in logit form
  psl <- unlogit(as.numeric(predict(lassoPsmod, xmat, s=lassoPsmod$lambda.1se)))
  names(psl) <- ds[,id]
  tryobjl <- try(lassoResults <- extractResults(ps=psl, exposurevec=ds[, exposed], 
                                                fmod=NULL, data=ds, id=id, exposed=exposed,
                                                outcome=outcome, logitFlag=TRUE,
                                                outfile=NULL, verbose=FALSE, 
                                                continuous=TRUE))
  if (class(tryobjl) != 'try-error') {
    lassoResults <- createPSResultsArray(lassoResults)
  } else {
    lassoResults <- c(n0=NA, n1=NA, effect_adj=NA, effect_matched=NA,
                      se_adj=NA, se_matched=NA, SMD=NA)
  }
  end <- proc.time()[1:3]
  time["lassoPS",,i] <- end-start
  
  
  ############# Causal Forests #############
  print("############# powers forest ############")
  start <- proc.time()[1:3]
  powerEffects <- causalForest(ds, method="power", covariates=empvariables,
                               name=name, estimators=estimators, desired='effect',
                               exposure=exposed, outcome=outcome,
                               confidence.interval=T)
  results <- c(n0=nrow(ds), n1=nrow(ds), 
               effect_adj=mean(powerEffects$effect, na.rm=T),
               effect_matched=NA, se_adj=NA, se_matched=NA, SMD=NA)
  powerResults <- list(results=results[preOrder], matchedID=NA)
  end <- proc.time()[1:3]
  time["powersForest",,i] <- end-start
  
  
  print("############ propensity forest ############")
  start <- proc.time()[1:3]
  propEffects <- causalForest(ds, covariates=empvariables, method="propensity",
                              name=name, estimators=estimators, desired='effect',
                              exposure=exposed, outcome=outcome,
                              confidence.interval=T)
  results <- c(n0=nrow(ds), n1=nrow(ds), 
               effect_adj=mean(propEffects$effect, na.rm=T),
               effect_matched=NA, se_adj=NA, se_matched=NA, SMD=NA)
  propResults <- list(results=results[preOrder], matchedID=NA)
  end <- proc.time()[1:3]
  time["propForest",,i] <- end-start 
  
  # reset current directory
  setwd('~/projects/confounderControl/realistic/')
  
  
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
  
  
  print("############ lassoMV model ###########")
  start <- proc.time()[1:3]
  
  xmat <- Matrix(as.matrix(ds[, c(exposed, empvariables)]), sparse=T)
  penalty.factor <- c(0, 0, rep(1, ncol(xmat)-2)) # force W and gen into model
  # tune lambda for glmnet using 5-fold CV
  lassoMvmod <- cv.glmnet(xmat, ds[, outcome], alpha=1, standardize=F, nfold=5,
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
  results <- c(n0=nrow(ds), n1=nrow(ds), SMD=Smd, se_adj=se_adj,
               effect_adj=lassomodCI$betaCIlim[exposed, 2],
               effect_matched=NA, se_matched=NA)
  lassoMVResults <- list(results=results[preOrder], matchedID=NULL)
  
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
  matchedID[[i]] <- list(naiveResults[[2]], expertResults[[2]], lassoResults[[2]], 
                         powerResults[[2]], propResults[[2]], 
                         euclideanResults[[2]], lassoMVResults[[2]], 
                         rdm1xResults[[2]])
  names(matchedID[[i]]) <- dimnames(time)[[1]]
  
  # consolidate results
  resultsmat <- as.data.frame(rbind(naiveResults[[1]], expertResults[[1]], 
                                    lassoResults[[1]], 
                                    powerResults[[1]], propResults[[1]], 
                                    euclideanResults[[1]], lassoMVResults[[1]], 
                                    rdm1xResults[[1]]))
  rownames(resultsmat) <- dimnames(time)[[1]]
  
  # computes MSE for a homogeneous estimate of treatment effect
  constant_error <- function(estimate, effects) {
    return(mean(abs(effects - estimate), na.rm=T))
  }
  
  resultsmat$bias_adj <- resultsmat$effect_adj-mean(tre, na.rm=T)
  resultsmat$bias_matched <- resultsmat$effect_matched-mean(tre, na.rm=T)
  resultsmat$mse_adj <- sapply(resultsmat$effect_adj, constant_error, tre)
  resultsmat$mse_matched <- sapply(resultsmat$effect_matched, constant_error, tre)
  # special case for causal forests: they can model heterogenerous effects
  resultsmat['powersForest', 'bias_adj'] <- mean(powerEffects$effect - tre, na.rm=T)
  resultsmat['powersForest', 'mse_adj'] <- mean(abs(powerEffects$effect - tre), na.rm=T)
  resultsmat['propForest', 'bias_adj'] <- mean(propEffects$effect - tre, na.rm=T)
  resultsmat['propForest', 'mse_adj'] <- mean(abs(propEffects$effect - tre), na.rm=T)
  resultsArray[,,i] <- as.matrix(resultsmat[, desiredOrder])
  # other results of interest
  heterogeneousEffects[,,i] <- as.matrix(data.frame(effect=tre,
                                                    powers=powerEffects$effect,
                                                    propensity=propEffects$effect,
                                                    powers.ci=powerEffects$ci,
                                                    propensity.ci=propEffects$ci))
  
  # save periodidically
  if (i %in% seq(50, 900, by=50)) {
    save(resultsArray, time, lassoMVvar,
         bestlambda, matchedID, Noutcomesmat, Nexposed,
         heterogeneousEffects, file=paste(name, 'results.RData', sep="_"))
  }
  
}

# remember to terminate the cores when done
stopCluster(cl)

save(ds, resultsArray, time, lassoMVvar, bestlambda, 
     matchedID, Noutcomesmat, Nexposed, heterogeneousEffects, 
     file=paste(name, 'results.RData', sep="_"))

dimnames(resultsArray)[1:2] <- list(rownames(resultsmat), desiredOrder)
dimnames(heterogeneousEffects)[1:2]<-list(1:n,c('effect','powers','propensity',
                                                'powers.ci', 'propensity.ci'))
save(ds, resultsArray, time, lassoMVvar, bestlambda, 
     matchedID, Noutcomesmat, Nexposed, heterogeneousEffects, 
     file=paste(name, 'results.RData', sep="_"))
