##############################################################################
# (Script)
# To HTN data, performs expertPS, hdPS, lasso logit, RF, Euclidean dist, 
# Jaccard, Dice,Cosine similarities, Pearson and Spearman correl, causal forests
#
# /!\ SHOULD BE RAN FROM COMMAND LINE, ELSE CALLING PYTHON SCRIPTS WITH SYSTEM
# DOES NOT WORK

# Note: HDPS was commented out as the Java package crashed on this dataset.
#
# Nathanael Romano
##############################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

# parallel core computing for lasso
require(doParallel)
require(foreach)
# specify number of cores req
cl=makeCluster(10) 
registerDoParallel(cl)

###### PARAMS #####
RScriptPath="/home/yenlow/scripts"
RProjectPath="~/projects/confounderControl/"
source(paste(RScriptPath, "/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath, "/R/patientChar.R", sep="")) # patientchar
setwd(paste(RProjectPath, 'HTN', sep=""))

baseDir = getwd();

installnewpackage(c("rJava", "Matching", "Epi", "glmnet", "randomForest", 
                    "vegan", "FNN", "Matrix", "doParallel", "foreach", "lars"))

require(rJava)
require(Matching)
require(Epi) #clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)
require(Matrix)
require(lars)

# initialize the Java subsystem.  include a jdbc object
# set JVM max heap size to to 12GB (set 12GB to get 60% of it, i.e. ~7GB)
# use pharmacoepi_juan.jar for more than 100 variables
#.jinit(classpath="/home/yenlow/software/pharmacoepi_juan.jar", 
# force.init = TRUE, parameters="-Xmx12g")
.jinit(classpath="/home/yenlow/software/pharmacoepi.jar", 
       force.init=TRUE, parameters="-Xmx20g")
.jclassPath()

#required functions
source(paste(RScriptPath,"/R/logit.R", sep=""))
source("../hdpsFunctions.R")
source("../rfFunctions.R")

# load data
load("htn_processed.RData", verbose=T)

id <- "pid"
exposed <- 'W' # name of exposure variable (must be string)
outcome <- c('Y.systolic', 'Y.diastolic')
expertvariables <- c()
forcedVar <- c()
empvariables <- setdiff(colnames(df), c(id, exposed, outcome)) ### TODO
empvariables_cat <- setdiff(empvariables, c("date.first.drug", "DM.age"))
empvariables_num <- c("date.first.drug", "DM.age")
print("Loaded data.")

# formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS = as.formula(paste(exposed, " ~ ", 
                       paste(sample(x=empvariables, size=20),
                             collapse=" + "), sep=""))

# remove bias and inCI (not applicable in this case since no ground truth)
# desiredOrder <- c("n1", "ORmatched", "ORlow_matched", "ORupp_matched",
#                   "coeff_matched", "se_matched", "n0", "ORadj", "ORlow_adj",
#                   "ORupp_adj", "coeff_adj", "se_adj", "SMD","KS","KL","Cstat")
desiredOrder <- c("n1", "coeff_matched", "se_matched", "n0", "coeff_adj",
                  "se_adj", "SMD", "KS", "KL", "Cstat")
# change variable name for consistency
ds <- df
# standardized
ds$date.first.drug <- (ds$date.first.drug - mean(ds$date.first.drug)) /
    sd(ds$date.first.drug)
ds$DM.age <- (ds$DM.age - mean(ds$DM.age)) / sd(ds$DM.age)

Noutcome <- length(outcome)
Ndraw <- 100
Nmethods <- 12 # 17 methods (incl. 3 rdm sampling, 2 causal forests, 1 naive)
nminor <- min(table(ds[,exposed]))
resultsArray <- array(dim=c(Nmethods, length(desiredOrder), Noutcome))
estimators <- 50 # number of estimators for causal forests
name <- 'HTN'
#Noutcomesmat=matrix(nrow=Noutcome,ncol=2)
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=NA
lassoMVvar=matchedID=hdpsobj=rdmBag=rdmJK=heterogeneousEffects=list()

### hdPS values depend on outcome; do within outcomes loop
# generate binary x matrix for hdps
dimdata=ds[,c(id, empvariables_cat)]
dim(dimdata)
time <- array(dim=c(Nmethods, 3, Noutcome))
dimnames(time)[[1]] <- c("naive", "expertPS", "hdPS", "lassoPS", "rfPS", 
                         "powerForest", "propForest", "euclidean",
                         "lassoMV", "rdmSamp1x", "rdmBag", "rdmJK")
dimnames(time)[[2]] <- c("user", "system", "elapsed")

##### expertPS #####
expertPsMod <- glm(fPS, data=ds, family="binomial")
ps_exp <- expertPsMod$fitted.values
names(ps_exp)=ds[,id]
print('Computed expert estimates of PS.')

##### lassoPS #####
xmat <- data.matrix(ds[, empvariables])
# set penalty to 0 to force certain varialbes into model
penalty.factor <- rep(1, ncol(xmat))
penalty.factor[colnames(xmat) %in% forcedVar] <- 0

# tune lambda for glmnet using 5-fold CV
lassoPsmod <- cv.glmnet(xmat, ds[, exposed], alpha=1, family="binomial",
                        standardize=F, nfold=5, penalty.factor=penalty.factor,
                        parallel=TRUE)
plot(lassoPsmod)
psl <- as.vector(predict(lassoPsmod, xmat, s=lassoPsmod$lambda.1se,
                         type="response"))
names(psl) <- ds[,id]
print('Computed lasso estimates of PS.')

##### rfPS #####
# if regression
# rfPsMod=randomForest(xmat,ds[,exposed], ntree=100, importance=F)
# ps_rf=predict(rfPsMod)
# hist(ps_rf)

# if classification
rfPsMod <- randomForest(xmat, as.factor(ds[,exposed]), ntree=10, importance=T)
ps_rf <- predict(rfPsMod, type="prob")[, 2]
range(ps_rf, na.rm=T)
# oob-adjust prob from RF to avoid 100% prob
ps_rf <- padjfromRF(rfPsMod)
hist(ps_rf)
names(ps_rf) <- ds[,id]
ps_rf[is.nan(ps_rf)] = NA
ps_rf[is.na(ps_rf)] = mean(ps_rf, na.rm=T)
rfImpt <- rfPsMod$importance[,1]
print('Computed rf estimates of PS.')

##### random downsampling #####
matchedID_rdm1x <- cbind(ds[ds[,exposed]==1,id],
                         sample(ds[ds[,exposed]==0,id],
                                nminor, replace=FALSE))

##### similarity #####
xmat_ctrl <- xmat[ds[,exposed]==0,]
xmat_trted <- xmat[ds[,exposed]==1,]
rownames(xmat_ctrl) <- ds[ds[,exposed]==0, id]
rownames(xmat_trted) <- ds[ds[,exposed]==1, id]

runSimPS<-function(method="jaccard", caliper=0.7, 
                   nsd=3, algorithm="kd_tree", i){
  matchedSim <- matchByDist(xmat_ctrl, xmat_trted, method=method, k_neighbors=5,
                            caliper=caliper, nsd=nsd, algorithm=algorithm)
  simResults <- extractResults(ps=matchedSim, exposurevec=NULL,
                               data=ds, fmod=NULL, id=id,exposed=exposed,
                               outcome=outcome[i], logitFlag=F,
                               outfile=NULL, verbose=FALSE, continuous=TRUE)
  return(simResults)
}


##### lassoMV #####
smat <- Matrix(cbind(ds[,exposed], xmat), sparse=T)
colnames(smat)[1] <- exposed
# force variables into model
penalty.factorMV <- rep(1,ncol(smat))
penalty.factorMV[colnames(smat) %in% forcedVar] <- 0

# save models
save(ds, rfPsMod, rfImpt, lassoPsmod, ps_exp, psl, 
     ps_rf, file="HTN_models.RData")
print('Saved preliminary models.')

#### Loop through 2 outcomes
for(i in 1:Noutcome){
  #for(i in c(1:9)){
  
  cat("\n\n-------", i, outcome[i],"-------\n")
  #if(i>1) file.remove(paste(tempDir, "output_cohort.txt", sep="/"))
  
  ############# naive estimation #############
  print("############# naive estimation ############")
  start <- proc.time()[1:3]
  
  results <- rep(NA, 16)
  names(results) <- c("n0", "n1", "KS", "KL", "Cstat", "SMD", "ORadj",
                      "ORlow_adj", "ORupp_adj", "coeff_adj", "se_adj",
                      "ORmatched", "ORlow_matched", "ORupp_matched",
                      "coeff_matched", "se_matched")
  results['n0'] <- 2000
  results['n1'] <- 2000
  results['coeff_adj'] <- mean(df[df$W==1, outcome[i]]) - 
                          mean(df[df$W==0, outcome[i]])
  results['ORadj'] <- exp(results['coeff_adj'])
  naiveResults <- list(results=results, matchedID=NA)
  end <- proc.time()[1:3]
  time["naive",,i] <- end-start
  
  
  ###### expertPS ######
  print("############ expertPS ############")
  start <- proc.time()[1:3]
  expertResults=extractResults(ps=ps_exp, exposurevec=ds[,exposed], fmod=NULL,
                               data=ds, id=id, exposed=exposed, 
                               outcome=outcome[i], logitFlag=TRUE,
                               outfile="HTN_expertPS.pdf", continuous=TRUE,
                               verbose=TRUE , printvalues=FALSE, ylim=NULL)
  
  end <- proc.time()[1:3]
  time["expertPS",,i]=end-start
  
  ####### hdPS ####### 
  print("############ hdPS ############")
  #tempDir <- paste(baseDir, "tmp", i, sep="/") # output directory
  #patientfile <- paste(tempDir, "patients.txt", sep="")
  #dir.create(tempDir)
  #
  # prepare data into format required by hdps
  # read data file as single string (to use addPatientsFromBuffer)
  #patientheader <- paste(id, exposed, outcome[i], sep="\t")
  #patientstring <- paste(paste(ds[,id], ds[,exposed], ds[,outcome[i]],
  #                             sep="\t"), collapse="\n")
  #datainstring <- paste(patientheader, patientstring, sep="\n")
  #
  ### INVOKE pharmacoepi.jar ### (handles cat variables only)
  #hdpsobj[[i]]=hdps(datainstring, dimdata, outDir=tempDir, Nmostfreq=80,
  #                  k=200, stratifyDim=FALSE, outfile="output_cohort.txt", 
  #                  FullOutput=FALSE, verbose=T, ZeroCellCorrection=F)
  # hdpsobj$selectedvariables # (see Fig2, Schneeweiss 2008 for definitions)
  #save(hdpsobj, file="hdpsobj.RData")
  # load(file="hdpsobj.RData",verbose=T)
  #wantedvar <- hdpsobj[[i]]$selectedvariables
  
  #var_corExposed <- empvariables_num[abs(cor(ds[,exposed],
  #                                           ds[,empvariables_num]))>0.05]
  #var_corOutcome <- empvariables_num[abs(cor(ds[,outcome[i]],
  #                                           ds[,empvariables_num]))>0.05]
  #IV <- setdiff(var_corExposed, var_corOutcome) # exclude from PS
  
  # Estimate the PS (force numerical variables into PS model)
  #dataPS <- cbind(ds[,c(id,exposed,var_corOutcome,forcedVar)],
  #                hdpsobj[[i]]$hdpsdata[,setdiff(wantedvar)])
  #hdPsMod <- glm(paste(exposed,"~ . -", id), data=dataPS, family="binomial")
  #names(hdPsMod$fitted.values) <- as.character(dataPS[,id])
  #summary(hdPsMod$fitted.values)
  
  #hdResults <- extractResults(ps=hdPsMod$fitted, continuous=TRUE,
                              #exposurevec=hdPsMod$y, fmod=NULL, data=ds, id=id,
                              #exposed=exposed, outcome=outcome[i],
                              #logitFlag=TRUE, outfile="HTN_hdPS.pdf",
                              #verbose=TRUE, printvalues=FALSE, ylim=NULL)
  hdResults <- list(results=rep(NA, 16), matchedID=NA)
  
  
  ####### lassoPS ####### 
  print("############ lassoPS ############")
  lassoResults <- extractResults(ps=psl, exposurevec=ds[,exposed], fmod=NULL,
                              data=ds, id=id, exposed=exposed,
                              outcome=outcome[i], logitFlag=TRUE,
                              outfile="HTN_lassoPS.pdf",  continuous=TRUE,
                              verbose=T, printvalues=FALSE, ylim=90)
  
  
  ####### rfPS ####### 
  print("############ rfPS ############")
  tryobj <- try(extractResults(ps=ps_rf, exposurevec=ds[,exposed],
                            data=ds, fmod=NULL, id=id, exposed=exposed,
                            outcome=outcome[i], logitFlag=TRUE, continuous=TRUE,
                            outfile="HTN_rfPS.pdf", verbose=TRUE,
                            printvalues=FALSE, ylim=5))
  if(class(tryobj)!="try-error") rfResults=tryobj else rfResults=list(NULL,NULL)
  
  
  ############# Causal Forests #############
  print("############# power forest ############")
  start <- proc.time()[1:3]
  
  powerEffects <- causalForest(ds, method="power", covariates=empvariables,
                           name=name, estimators=estimators, output='effect',
                           exposure='W', outcome=outcome[i])
  results <- rep(NA, 16)
  names(results) <- c("n0", "n1", "KS", "KL", "Cstat", "SMD", "ORadj",
                      "ORlow_adj", "ORupp_adj", "coeff_adj", "se_adj",
                      "ORmatched", "ORlow_matched", "ORupp_matched",
                      "coeff_matched", "se_matched")
  results['n0'] <- 2000
  results['n1'] <- 2000
  results['coeff_adj'] <- mean(powerEffects, na.rm=T)
  results['ORadj'] <- exp(results['coeff_adj'])
  powerResults <- list(results=results, matchedID=NA)
  
  end <- proc.time()[1:3]
  time["powerForest",,i] <- end-start
  
  
  print("############ propensity forest ############")
  start <- proc.time()[1:3]
  
  propEffects <- causalForest(ds, covariates=empvariables, method="propensity",
                           name=name, estimators=estimators, output='effect',
                           exposure='W', outcome=outcome[i])
  results <- rep(NA, 16)
  names(results) <- c("n0", "n1", "KS", "KL", "Cstat", "SMD", "ORadj",
                      "ORlow_adj", "ORupp_adj", "coeff_adj", "se_adj",
                      "ORmatched", "ORlow_matched", "ORupp_matched",
                      "coeff_matched", "se_matched")
  results['n0'] <- 2000
  results['n1'] <- 2000
  results['coeff_adj'] <- mean(propEffects, na.rm=T)
  results['ORadj'] <- exp(results['coeff_adj'])
  propResults <- list(results=results, matchedID=NA)
  
  end <- proc.time()[1:3]
  time["propForest",,i] <- end-start  
 

  ####### similarity ####### 
  print("############ similarity ############")
  euclideanResults=runSimPS(method="euclidean", nsd=3, algorithm="brute", i=i)

  
  ####### lassoMV ####### 
  print("############ lassoMV ############")
  xmat <- Matrix(as.matrix(ds[, empvariables]), sparse=T)
  # tune lambda for glmnet using 5-fold CV
  lassoMvmod <- cv.glmnet(xmat, ds[,outcome], alpha=1, family="gaussian",
                          standardize=F, nfold=5, penalty.factor=penalty.factor)
  coeff <- na.omit(sort(coeffAtlambda(lassoMvmod)))
  
  Smd=smd(ds[,c(exposed, outcome[i])], exposed=exposed, variable=outcome[i],
          verbose=FALSE, categorical=TRUE)
  results <- c(n0=nrow(smat), n1=NA, 
               KS=NA, KL=NA, Cstat=NA, SMD=Smd, 
               ORadj=NA, ORlow_adj=NA, ORupp_adj=NA, 
               coeff_adj=coeff, se_adj=NA, rep(NA, 5))
  lassoMVResults <- list(results=results, matchedID=NULL)
  
  
  ####### random sampling ensemble ####### 
  print("############ random downsampling 1x ############")
  rdm1xResults <- extractResultsRdm(matchedID_rdm1x, ds, fmod=NULL, id=id,
                                    exposed=exposed, outcome=outcome[i],
                                    verbose=FALSE, continuous=TRUE)
  
  for(j in 1:Ndraw){
    ############### random bagging (repeatedly downsample major class Ndraw 
    # times WITH replacement) ############
    matchedID_bag <- cbind(ds[ds[,exposed]==1,id], 
                           sample(ds[ds[,exposed]==0,id], nminor, replace=TRUE))
    rdmBag[[j]] <- extractResultsRdm(matchedID_bag, ds, fmod=NULL, id=id,
                                     exposed=exposed, outcome=outcome[i],
                                     verbose=FALSE, continuous=TRUE)
    
    ############### random jackknife (repeatedly downsample major class Ndraw 
    # times WITHOUT replacement) ############
    matchedID_JK <- cbind(ds[ds[,exposed]==1,id],
                          sample(ds[ds[,exposed]==0,id], nminor))
    rdmJK[[j]] <- extractResultsRdm(matchedID_JK, ds, fmod=NULL, id=id,
                                    exposed=exposed, outcome=outcome[i],
                                    verbose=FALSE, continuous=TRUE)
  }
  
  temp <- matrix(unlist(lapply(rdmBag, `[`,"results")),
                 nrow=Ndraw, ncol=length(rdmBag[[1]][[1]]),byrow=T)
  colnames(temp) <- names(rdmBag[[1]][[1]])
  rdmBagResults<- list(results=colMeans(temp,na.rm=T), matchedID=NULL)
  
  temp <- matrix(unlist(lapply(rdmJK,`[`,"results")),
                 nrow=Ndraw, ncol=length(rdmJK[[1]][[1]]), byrow=T)
  colnames(temp) <- names(rdmJK[[1]][[1]])
  rdmJKResults <- list(results=colMeans(temp,na.rm=T), matchedID=NULL)
  
  # consolidate matchedIDs
  matchedID[[i]] <- list(naiveResults[[2]], expertResults[[2]], hdResults[[2]], 
                      lassoResults[[2]], rfResults[[2]], powerResults[[2]],
                      propResults[[2]], euclideanResults[[2]],
                      lassoMVResults[[2]], rdm1xResults[[2]],
                      rdmBagResults[[2]], rdmJKResults[[2]])
  names(matchedID[[i]]) <- c("naive", "expertPS", "hdPS", "lassoPS", "rfPS",
                          "power", "propensity", "euclidean",
                          "lassoMV", "rdmSamp1x", "rdmBag", "rdmJK")
  heterogeneousEffects[[i]] <- data.frame(power=powerEffects,
                                          propensity=propEffects)
  names(heterogeneousEffects[[i]]) <- c('power', 'propensity')

  ####### Consolidate results #######
  resultsmat <- as.data.frame(rbind(naiveResults[[1]], expertResults[[1]], 
                                    hdResults[[1]], lassoResults[[1]], 
                                    rfResults[[1]], powerResults[[1]],
                                    propResults[[1]],
                                    euclideanResults[[1]],
                                    lassoMVResults[[1]], rdm1xResults[[1]],
                                    rdmBagResults[[1]], rdmJKResults[[1]]))
  rownames(resultsmat) <- names(matchedID[[i]])
  resultsArray[,,i] <- as.matrix(resultsmat[, desiredOrder])
  
}   # end of Noutcomes
dimnames(resultsArray) <- list(rownames(resultsmat), desiredOrder, outcome)
names(matchedID) <- outcome
names(heterogeneousEffects) <- outcome

save(rfPsMod, rfImpt, lassoPsmod, lassoMvmod,
     ps_exp, psl, ps_rf, hdpsobj, file="HTN_models.RData")
save(resultsArray, matchedID, heterogeneousEffects, file="HTN_results.RData")

# remember to terminate the cores when done
stopCluster(cl)
# resultsArray
# 
# round(apply(resultsArray,c(1,2),mean,na.rm=T),3)
# apply(time[,"elapsed",],1,summary)
# apply(rfImpt,1,summary)
