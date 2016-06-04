################################################################################
# (Script)
# Run Scott Powers' simulation 1, only with Causal Forests, to explore their
# results.
#
# /!\ SHOULD BE RAN FROM COMMAND LINE, ELSE CALLING PYTHON SCRIPTS WITH SYSTEM
# DOES NOT WORK
#
################################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath="/home/yenlow/scripts"
source(paste(RScriptPath, "/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath, "/R/patientChar.R", sep="")) # patientchar
name <- 'CF_1000'

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
estimators <- 100 # number of estimators for causal forests

# genre variable
gen <- c(rep(0,n/2), rep(1,n/2))

# generate PS
ps <-rep(NA, n)
ps[gen == 0] <- .25
ps[gen == 1] <- .75

# generate treatment
tr <- rep(NA, n)
tr[gen == 0] <- 1*(runif(n/2) < .25)
tr[gen == 1] <- 1*(runif(n/2) < .75)

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
nminor <- min(table(ds[, exposed]))

############# Causal Forests #############
print("############# power forest ############")
powerEffects <- causalForest(ds, method="power", covariates=empvariables,
                             name=name, estimators=estimators,
                             exposure=exposed, outcome=outcome,
                             confidence.interval=TRUE)


print("############ propensity forest ############")
propEffects <- causalForest(ds, covariates=empvariables, method="propensity",
                            name=name, estimators=estimators,
                            exposure=exposed, outcome=outcome,
                            confidence.interval=TRUE)

# reset current directory
setwd('~/projects/confounderControl/scott/')
heterogeneousEffects <- as.matrix(data.frame(effect=tre,
                                            power=powerEffects$effect,
                                            power.ci=powerEffects$ci,
                                            propensity=propEffects$effect,
                                            propensity.ci=propEffects$ci))

# remember to terminate the cores when done
stopCluster(cl)
dimnames(heterogeneousEffects) <- list(1:n, c('effect', 'power', 'power.ci', 
                                              'propensity', 'propensity.ci'))
save(ds, heterogeneousEffects, file=paste(name, 'results.RData', sep="_"))
