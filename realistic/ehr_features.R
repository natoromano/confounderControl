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
# Plots bias versus confounding strength
#
# /!\ SHOULD BE RAN FROM COMMAND LINE, ELSE CALLING PYTHON SCRIPTS WITH SYSTEM
# DOES NOT WORK
#
# For the new methods, we require at least 5 samples of each treatment class
################################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath <- "/home/yenlow/scripts"
RProjectPath <- "~/projects/confounderControl/"
source(paste(RScriptPath,"/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath,"/R/patientChar.R", sep="")) # patientchar

require("ggplot2")

# inputs and parameters
name <- 'masked_age_observed_09'
nsim <- 100
n_EHR_var <- 100
alps <- c(0, 1, 1, 2, 5) # degree of confounding
desired.match.methods <- c('lassoPS', 'expertPS', 'euclidean')
desired.adj.methods <- c('naive', 'expertPS', 'lassoPS') #'propForest', 'powersForest')
Nmethods <- 8 # change after adding the two extra methods
resultsTotal <- array(dim=c(Nmethods, 4, length(alps)))

for (k in 1:length(alps)) {
  
  alp <- alps[k]
  tryobj <- try(source('~/projects/confounderControl/realistic/sim.R'))
  
  if (class(tryobj)!='try-error') {
    # remove extreme values for expertPS
    invalidID <- which(resultsArray["expertPS", "se_matched",]>10)
    resultsArray["expertPS", , invalidID] <- NA
    
    # fetch relevant values
    resultsTotal[, 1, k] <- abs(apply(resultsArray, 
                                  c(1,2), mean, na.rm=T)[,'bias_matched'])
    resultsTotal[, 2, k] <- abs(apply(resultsArray, 
                                  c(1,2), mean, na.rm=T)[,'bias_adj'])
    resultsTotal[, 3, k] <- apply(resultsArray, 
                             c(1,2), sd, na.rm=T)[,'bias_matched'] / sqrt(nsim)
    resultsTotal[, 4, k] <- apply(resultsArray, 
                             c(1,2), sd, na.rm=T)[,'bias_adj'] / sqrt(nsim)
    # save data
    file <- paste(name, ".RData", sep="")
    save(resultsTotal, file=paste(RProjectPath, 
                                  paste('realistic', file, sep="/"), sep=""))
  }
  
}

resultsTotal <- resultsTotal[,, -2]
alps = alps[-2]
# change column names
resultsNames <- c('bias_matched', 'bias_adj', 'se_matched', 'se_adj')
dimnames(resultsTotal) <- list(rownames(resultsArray), resultsNames, alps)

# save again
file <- paste(name, ".RData", sep="")
save(resultsTotal, file=paste(RProjectPath, 
                              paste('realistic', file, sep="/"), sep=""))

# process results to get plot
n_var = rep(alps, length(desired.match.methods) + length(desired.adj.methods))
biases = c()
methods = c()
se = c()
i <- 0
for (m in desired.match.methods) {
  for (al in alps) {
    i <- i + 1
    biases[i] <- abs(resultsTotal[m, 'bias_matched', as.character(al)])
    se[i] <- resultsTotal[m, 'se_matched', as.character(al)]
    if (grepl('PS', m)) methods[i] <- paste(m, 'match', sep="_")
    else methods[i] <- m
  }
}

for (m in desired.adj.methods) {
  for (al in alps) {
    i <- i + 1
    biases[i] <- abs(resultsTotal[m, 'bias_adj', as.character(al)])
    se[i] <- resultsTotal[m, 'se_adj', as.character(al)]
    if (grepl('PS', m)) methods[i] <- paste(m, 'adj', sep="_")
    else methods[i] <- m  }
}

# combine results
df = data.frame(nvar=n_var, bias=biases, se=se, method=methods)

# plot results
bias.plot <- ggplot(data = df, mapping = aes(x = nvar, y = bias, 
                                             group = method, color = method)) +
  geom_line() + geom_point() + xlab("confounding strength (alpha)") +
  theme_bw() + theme(legend.key = element_blank()) +
   geom_errorbar(aes(ymin=bias-se, ymax=bias+se), width=.9,
               position=position_dodge(0.05)) + 
  ggtitle("Evolution of bias with confounding strength")

# save plot
file <- paste(name, ".png", sep="")
filename <- paste(RProjectPath, paste('realistic/', file, sep="/"), sep="")
ggsave(filename=filename, plot=bias.plot)
