###############################################################################
# (Script)
# To simulated data 2000x10, do expertPS, hdPS, lassoPS, rfPS (adj and match)
# Euclidean dist, Jaccard, Dice, Cosine, Pearson and Spearman correllation.
# Does so for various confounding strengths (can be specified in the alps
# vector).
# Plot bias vs confounding strength for specified methods and saves plot.
#
# /!\ SHOULD BE RAN FROM COMMAND LINE, ELSE CALLING PYTHON SCRIPTS WITH SYSTEM
# DOES NOT WORK
#
# Nathanael Romano
##############################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath <- "/home/yenlow/scripts"
RProjectPath <- "~/projects/confounderControl/"
source(paste(RScriptPath,"/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath,"/R/patientChar.R", sep="")) # patientchar

require("ggplot2")

# inputs and parameters
name <- 'full_correlation_not_masked'
alps <- c(0, 1, 2, 5) # alpha values for confounding strength
mask <- 0 # (optionally) set masked covariate
desired.match.methods <- c('expertPS','euclidean', 'rdmSamp1x')
desired.adj.methods <- c('lassoMV', 'expertPS', 'propForest', 'powerForest')
Nmethods <- 16 # change after adding the two extra methods
resultsTotal <- array(dim=c(Nmethods, 4, length(alps)))

for (k in 1:length(alps)) {
  
  alp <- alps[k]
  tryobj <- try(source('~/projects/confounderControl/models/sim2000x10/ps.R'))
  if (class(tryobj)!='try-error') {
    # remove extreme values for expertPS
    invalidID <- which(resultsArray["expertPS", "se_matched",]>10)
    resultsArray["expertPS", , invalidID] <- NA
    
    # fetch relevant values
    resultsTotal[, 1, k] <- apply(resultsArray, 
                                        c(1,2), mean, na.rm=T)[,'bias_matched']
    resultsTotal[, 2, k] <- apply(resultsArray, 
                                        c(1,2), mean, na.rm=T)[,'bias_adj']
    resultsTotal[, 3, k] <- apply(resultsArray, 
                                        c(1,2), sd, na.rm=T)[,'bias_matched']
    resultsTotal[, 4, k] <- apply(resultsArray, 
                                        c(1,2), sd, na.rm=T)[,'bias_adj']
    # save data
    file <- paste(name, ".RData", sep="")
    save(resultsTotal, file=paste(RProjectPath, 
                                  paste('experiments', file, sep="/"), sep=""))
    }

}

# change column names
resultsNames <- c('bias_matched', 'bias_adj', 'se_matched', 'se_adj')
dimnames(resultsTotal) <- list(rownames(resultsArray), resultsNames, alps)

# save again
file <- paste(name, ".RData", sep="")
save(resultsTotal, file=paste(RProjectPath, 
                              paste('experiments', file, sep="/"), sep=""))

# process results to get plot
alphas = rep(alps, length(desired.match.methods) + length(desired.adj.methods))
biases = c()
methods = c()
se = c()
i <- 0
for (m in desired.match.methods) {
  for (al in alps) {
    i <- i + 1
    biases[i] <- abs(resultsTotal[m, 'bias_matched', as.character(al)])
    se[i] <- resultsTotal[m, 'se_matched', as.character(al)]
    methods[i] <- paste(m, 'match', sep="_")
  }
}
for (m in desired.adj.methods) {
  for (al in alps) {
    i <- i + 1
    biases[i] <- abs(resultsTotal[m, 'bias_adj', as.character(al)])
    se[i] <- resultsTotal[m, 'se_adj', as.character(al)]
    methods[i] <- paste(m, 'adj', sep="_")
  }
}

# combine results
df = data.frame(alpha=alphas, bias=biases, se=se, method=methods)

# plot results
bias.plot <- ggplot(data = df, mapping = aes(x = alpha, y = bias, 
                                group = method, color = method)) +
  geom_line() + geom_point() + xlab("confounding factor (alpha)") +
  theme_bw() + theme(legend.key = element_blank()) +
  # geom_errorbar(aes(ymin=bias-se, ymax=bias+se), width=.9,
  #            position=position_dodge(0.05)) + 
  ggtitle("Evolution of bias with confounding strength")

# save plot
file <- paste(name, ".png", sep="")
filename <- paste(RProjectPath, 
                  paste('experiments/plots/biases', file, sep="/"), sep="")
ggsave(filename=filename, plot=bias.plot)
