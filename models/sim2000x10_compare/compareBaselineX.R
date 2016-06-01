##############################################################################
# To sim2000x10 data, compare how well baseline characteristics match
# after expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman 
# correl and rdm sampling
#
# 23-Oct-14 fix IV to just x7
# 15-Oct-14 use rfPS in reg mode (set min node size to 5% of sample size, i.e. 
# 100, Malley 2012) 
# 17-Sep-14 moved rdmSamp1x to last row
# 26-Aug-14 Yen Low
##############################################################################

############# SET-UP #############
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath <- "/home/yenlow/scripts"
source(paste(RScriptPath,"/R/utils.R", sep="")) # installnewpackage
ProjectPath <- "~/projects/confounderControl"
setwd(paste(ProjectPath, "/models/sim2000x10_compare", sep=""))
getwd()

installnewpackage(c('ggplot2', 'reshape2', 'grid'))

require(ggplot2)
require(reshape2)
require("grid")

source("../../charts.R")
source("../../hdpsFunctions.R")
source(paste(RScriptPath, "/R/utils.R", sep=""))
source(paste(RScriptPath, "/R/patientChar.R", sep=""))
source(paste(RScriptPath, "/R/logit.R", sep=""))

############# LOAD AND PROCESS RESULTS #############
# load results from running various methods on PAD_10sparseCid
load(file="../sim2000x10/sim2000x10_1000MC.RData", verbose=T)

catVar <- c("x1", "x3", "x5", "x6", "x8", "x9")
numVar <- c("x2", "x4", "x7", "x10")
expertvariables <- c("x1", "x2", "x3", "x4", "x5", "x6")

modelOrder <- c("expertPS", "hdPS", "lassoPS", "rfPS",
             "euclidean", "jaccard", "dice", "cosine", "pearson", "spearman",
             "lassoMV", "bootstrap", "jackknife", "rdmSamp1x")


# check for large values (possibly due to quasi separation, esp in expertPS)
tail(sort(resultsArray["expertPS", "se_matched",]))
invalidID <- which(resultsArray["expertPS", "se_matched",]>10)
resultsArray["expertPS", , invalidID]
save(resultsArray, file="../sim2000x10/resultsArray_untrimmed.RData")
# trim invalid models
resultsArray["expertPS", , invalidID] <- NA
save(resultsArray, file="../sim2000x10/resultsArray_trimmed.RData")

# get results from all methods
apply(resultsArray, c(1,2), mean, na.rm=T)


############# SELECTED VARIABLES #############
# number of times variables were used in hdPS model
hdVarID <- sapply(wantedVarList,substr, start=2, stop=3)
mode(hdVarID) <- "numeric"
timesSel_hdps <- c(table(catVar[unlist(hdVarID)]),
                   table(unlist(var_corOutcomeList)))
timesSel_hdps <- timesSel_hdps[paste("x", 1:10, sep="")]

png("sim2000x10_hdPSImpt.png", width=11, height=7, 
    units="in", res=300, bg="transparent")
barplot(timesSel_hdps[10:1], horiz=TRUE, las=1, 
        main="Number of times selected by hdPS")
dev.off()


# rf variable importance of rfPS model
rownames(rfImpt) <- paste("x", 1:10, sep="")
varImp <- rowMeans(rfImpt)
sort(varImp,dec=T)

png("sim2000x10_rfImpt.png", width=11, height=7,
    units="in", res=300, bg="white")
dotplot(rfImpt[nrow(rfImpt):1,], pch="|", col=1, cex=0.8, 
        main="RF variable importance")
dev.off()

# number of times variables were used in lassoMV model
sellassoMVvar <- lapply(lassoMVvar, 
              function(x) setdiff(rownames(x)[x[,2]!=0 & (x[,3]<=0|x[,1]>=0)],
                                              c("(Intercept)", "exposure")))

# extract pvalues and smd of baseline variables between groups
pmat <- apply(-log10(pmatArray), c(1,2), mean, na.rm=T)
pmat[is.na(pmat)] <- NA
pmat <- pmat[c("before",modelOrder),]

smdmat <- apply(smdmatArray, c(1,2), mean, na.rm=T)
smdmat <- smdmat[c("before", modelOrder),]

# number of times variabls are selected
varselmat <- matrix(NA,nrow=length(modelOrder), ncol=ncol(pmat),
                    dimnames=list(modelOrder, colnames(pmat)))
varselmat["expertPS",] <- 0
varselmat["expertPS",expertvariables] <- 1000
varselmat["hdPS",] <- timesSel_hdps
varselmat["lassoPS",] <- rowSums(!is.na(lassoPSBeta))
varselmat["rfPS",] <- rowSums(rfImpt>0)
varselmat[c(5:10,12:14),] <- 1000
varselmat["lassoMV",] <- table(unlist(sellassoMVvar))[paste("x", 1:10, sep="")]
varselmat["lassoMV", is.na(varselmat["lassoMV",])] <- 0


############# HEATMAPS #############
colscale=colorpanel(10, "white", "black")

png("sim2000x10_heatmap_pval.png", width=6.5, height=7, units="in",
    res=150, bg="transparent")
heatmap.2(pmat, Rowv=F, Colv=F, scale="none", col=colscale, trace="none",
          keysize=1, density.info="none", na.rm=T, na.col="transparent")
# heatmap.2(t(1-pmat),Rowv=F,Colv=F,col=colscale,scale="none",
#          trace="none",keysize=1,density.info="none")
dev.off()

png("sim2000x10_heatmap_smd.png", width=6.5, height=7, units="in", 
    res=150, bg="transparent")
heatmap.2(abs(smdmat), Rowv=F, Colv=F, scale="none", col=colscale, trace="none",
          keysize=1, density.info="none")
dev.off()

png("sim2000x10_heatmap_varsel.png", width=6.5, height=7, units="in", 
    res=150, bg="transparent")
heatmap.2(varselmat/1000, Rowv=F, Colv=F, scale="none", col=colscale,
          trace="none", keysize=1, density.info="none")
dev.off()

