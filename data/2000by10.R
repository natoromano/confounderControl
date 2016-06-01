###############################################################################
# (Script) 
# Generate simulated data set 2000 patients x 100 variables (additional IV, 
# risk factors, colliders, interrcorrel noise)
# and 2000 patients x 10 variables (Setoguchi 2008, Scenario E)
# 
# 30-Mar-14 Yen Low
# 22-Feb-16 refactoring, Nathanael Romano
###############################################################################

############# SETUP #############
# rm(list=ls(all=T))
source("~/projects/confounderControl/hdpsFunctions.R")

############# SET PARAMETERS #############
n <- 2000 # number of patients
p <- 10 # number of variables
# set logit(exposure) coefficients; includes b0
b <- c(0, alp * 0.8, alp * -0.25, alp * 0.6, alp * -0.4,
       -0.8, -0.5, 0.7, 0, 0, 0) # change 2 and 3
# set logit(outcome) coefficients; includes a0
a <- c(-3.85, alp * 0.3, alp * -0.36, alp * -0.73, alp * -0.2,
       0, 0, 0, 0.71, -0.19, 0.26) # change 3 and 4
gamma <- -0.4 # coeff for PS covariate adjustment
minSum <- 40 # ensure minimum of 40 outcomes

############# SIMULATE #############
while (minSum<=40) { # ensure minimum of 40 outcomes
  # generate 10 variables based on standard normal distribution
  # introduce intercorrelation between 0.2 and 0.9
  v1 <- rnorm(n,0,1)
  v3 <- rnorm(n,0,1)
  
  x1 <- as.integer(v1>=0)
  x2 <- rnorm(n,0,1)
  x3 <- as.integer(v3>=0)
  x4 <- rnorm(n,0,1)
  
  v5 <- 1*((1-sqrt(0.7))*v1+sqrt(0.7)*rnorm(n,0,1))
  v6 <- -1*((1-sqrt(0.0))*x2+sqrt(0.0)*rnorm(n,0,1))
  v8 <- -1*((1-sqrt(0.7))*v3+sqrt(0.7)*rnorm(n,0,1))
  v9 <- 1*((1-sqrt(0.1))*x4+sqrt(0.1)*rnorm(n,0,1))
  
  x5 <- as.integer(v5>=0)
  x6 <- as.integer(v6>=0)
  x7 <- rnorm(n,0,1) # strong IV
  x8 <- as.integer(v8>=0)
  x9 <- as.integer(v9>=0)
  x10 <- rnorm(n,0,1) # weak risk factor
  
  xmat <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  # dim(xmat)
  
  # view intercorrelations
  cormat=cor(xmat)
  write.table(cormat, file="100x100cormat.txt", sep="\t",
              col.names=F, row.names=F)
  
  png("100x100cormat.png", width=7, height=7, res=100, units="in",
      bg="transparent")
  heatmap(cormat, Rowv=NA, symm=T, col=cm.colors(10), scale="none")
  dev.off()
  
  cormat[abs(cormat)<0.1] <- NA # suppress for easy viewing
  save(xmat, file="xmat2000by10.RData")
  
  # simulate propensity scores based on scenario E from eq in Appendices 1-2
  # handle b0 (i.e. b[1]) differently
  Pexp <- c(1/(1+exp(-(b[1] + b[-1] %*% t(xmat)))))
                     #+b[3]*x2*x2
                     #+0.5*b[2]*x1*x3
                     #+0.7*b[3]*x2*x4
                     #+0.5*b[5]*x4*x5
                     #+0.5*b[6]*x5*x6))))
  # sum(mapply(function(x,y) x*y, b[-1], xmat[2,]))
  
  # simulate true exposure
  rdmNum <- runif(n)
  temp <- as.data.frame(cbind(rdmNum, Pexp))
  temp$exposure <- 0
  temp$exposure[temp$rdmNum<temp$Pexp] <- 1
  
  # simulate estimated prob of outcome based on eq in Appendices 1-2
  Poutcome=c(1/(1+exp(-(a[1]+a[-1] %*% t(xmat) + gamma*temp$exposure))))
  
  # simulate true outcome
  temp$rdmNum2 <- runif(n)
  temp$Poutcome <- Poutcome
  temp$outcome <- 0
  temp$outcome[temp$rdmNum2<temp$Poutcome] <- 1
  
  # check that there are at least 40 outcomes
  tab <- table(temp$outcome,temp$exposure)
  Noutcomes <- tab[2,]
  minSum <- sum(Noutcomes)
  cat("minSum: ", minSum,"\n")
}

############# PUTTING MATRICES TOGETHER #############
ds10 <- cbind(outcome=temp$outcome, exposure=temp$exposure, ps=temp$Pexp,
              Poutcome=temp$Poutcome, as.data.frame(xmat))
		
############# SAVE #############
# get true OR (0.67)
cat("Pre-defined OR: ", exp(gamma), '\n') # pre-defined OR (0.67)
mod <- glm(paste("outcome ~ exposure +", paste("x", c(1:4, 8:10), 
                                               collapse=" + ", sep="")),
				data=ds10, family="binomial")
trueCoeff <- coefficients(mod)["exposure"]
trueOR <- exp(trueCoeff)

# get crude OR (0.58)
# print(chisq.test(tab))
print(tab)
OR <- tab[1,1]*tab[2,2]/tab[1,2]/tab[2,1]
cat("True estimated OR:", trueOR, "\n") # actual re-derived OR from data = 0.45
cat("Crude OR: ", OR, "\n")

# get true standardized mean difference (-0.09)
trueSmd <- smd(ds10, exposed="exposure", variable="outcome",
               verbose=FALSE, categorical=TRUE)
cat("SMD: ", trueSmd, "\n\n")

if (mask != 0) {
  masked.var = paste('x', mask, sep="")
  ds10$masked = ds10[, masked.var]
  ds10[, masked.var] = NA
}
save(ds10, temp, file="2000x10.RData")

# check underlying distributions
png("2000x100histograms.png", width=7, height=11, res=100,
    units="in", bg="transparent")
par(mfrow=c(3, 2))
hist(Pexp)
h <- hist(logit(Pexp), prob=T)
lines(density(logit(Pexp)))

hist(Poutcome)
hist(logit(Poutcome), prob=T)
lines(density(logit(Poutcome)))

hist(temp$exposure)
hist(temp$outcome)
dev.off()
