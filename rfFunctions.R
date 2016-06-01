###############################################################################
# (Functions)
# Functions to generate ORs from random forest, as well as functions
# implementing direct effect estimation through Causal Forests.
#
# 29-Feb-16 refactoring, Nathanael Romano
# 08-Sep-14 Yen Low
##############################################################################
RProjectPath <- "~/projects/confounderControl/"

############# SET-UP #############
require(randomForest)
source("~/projects/confounderControl/hdpsFunctions.R")  # for smd and logit

############# FUNCTIONS FOR POWERS AND PROPENSITY FOREST #############

# causalForests
# -------------
# Arguments:
# - data: dataset used, in a dataframe format
# - name: name of the experiment
# - method: either "power" or "propensity", for Power Forests or Propensity
# Forests
# - desired: 'effect' or 'outcome', typical case is effect
# - estimators: number of base learners
# - exposure, outcome: exposure and outcome covariate name
# - confidence.interval: if CI is desired, will be added to returned dataframe
# - train: separate training set or NULL
causalForest <- function(data, name, covariates=NULL, 
                         method="power", estimators=100,
                         exposure='exposure', desired='effect',
                         outcome='outcome', confidence.interval=FALSE,
                         train=NULL) {
  # Outputs a data frame with columns mean(treated outcome) and mean(control 
  # outcome), one row per patient
  if (is.null(covariates)) covariates = paste('x', 1:10, sep="")
  setwd(paste(RProjectPath, 'experiments', sep=""))
  orderedData = data[c(exposure, outcome, covariates)]
  filename <- paste(paste(name, outcome, sep="_"), '.txt', sep="")
  filepath <- paste(RProjectPath, 'experiments/', filename, sep="")
  write.table(orderedData, file=filename, col.names=FALSE, row.names=FALSE,
              quote=FALSE)
  
  if (!is.null(train)) {
    train_ordered <- train[c(exposure, outcome, covariates)]
    filename_train <- paste(paste(name, outcome, 'train', sep="_"), 
                            '.txt', sep="")
    filepath_train <- paste(RProjectPath, 'experiments/', filename_train, sep="")
    write.table(train_ordered, file=filename_train, 
                col.names=FALSE, row.names=FALSE,
                quote=FALSE)
  } else {
    filename_train <- filename
  }
  
  if (method=='power') fct <- 'fit_power_forest.py'
  else fct <- 'fit_propensity_forest.py'
  output_filename <- paste('output', outcome, method, filename, sep="_")
  if (!confidence.interval) ci <- 'no'
  else ci <- paste('output', outcome, method, 'ci', filename, sep="_")
  command <- paste('python', fct, filename_train, filename, 
                   output_filename, estimators, desired, ci)
  system(command)
  
  if (desired == 'effect') {
    output <- read.delim(file=paste(RProjectPath, 'experiments/', 
                                    output_filename, sep=""), 
                         col.names=c('effect'), header=FALSE)
    output <- data.frame(effect = as.numeric(as.vector(output$effect)))
  } else {
    output <- read.delim(file=paste(RProjectPath, 'experiments/', 
                                    output_filename, sep=""), 
                         col.names=c('treated', 'control'), header=FALSE)
    output <- data.frame(treated = as.numeric(as.vector(output$treated)),
                         control = as.numeric(as.vector(output$control)))
  }
  
  if (ci != 'no') {
    confint <- read.delim(file=paste(RProjectPath, 'experiments/', ci, sep=""),
                        col.names=c('ci'), header=FALSE)
    output$ci <- as.numeric(as.vector(confint$ci))
  }
  return(output)
}


############# FUNCTIONS FOR RFMV #############
padjfromRF <- function(rfmod, p=NULL, class=NULL){
  if (is.null(class)) class <- rfmod$classes[2]
  if (is.null(p)) p <- rfmod$votes[,class]
  oob <- sum(rfmod$confusion[,3] * rowSums(rfmod$confusion[,-3]) / 
               sum(rfmod$confusion[,-3]))
  padj <- p - (2*p - 1) * oob
  return(padj)
}

ORfromRF <- function(xmat, yvec, xvar, ntree=100, class=NULL, balance=TRUE,
                     verbose=TRUE){
  # rf set-up
  if (balance==TRUE) {
    rfMV <- randomForest(xmat, as.factor(yvec), ntree=ntree, importance=T,
                         proximity=T, sampsize=rep(min(table(yvec)), 2))
  } else if(balance==FALSE) {
    rfMV <- randomForest(xmat, as.factor(yvec), ntree=ntree, importance=T,
                         proximity=T)
  } else if (balance.is.null()) { # regression mode
    rfMV <- randomForest(xmat, as.numeric(as.character(yvec)), ntree=ntree, 
                         importance=T, proximity=T)  
  } else {
    stop("balance must be: \n
          TRUE (recommended, randomly bootstrap majority class to match 
          minority class), \n
          FALSE (to leave classes unbalanced), \n
          or NULL (for RF in regression mode)")
  }
  
  if (is.null(class)) class <- rfMV$classes[2]
  if (is.factor(yvec)) { 
    numLevels <- nlevels(yvec)
    # returns logit(probabilities of class "1" at each level of "cilostazol")
    # ESL/Hastie section 10.13.2 partial dependence
    # http://stats.stackexchange.com/questions/93202/odds-ratio-from-decision-
    # tree-and-random-forest
    pdep <- unlogit(partialPlot(rfMV, xmat, as.character(xvar), 
                                which.class=class)$y)
    if (numLevels==2) {
      OR <- (pdep[2]/(1-pdep[2]))/(pdep[1]/(1-pdep[1]))
      beta <- log(OR)
    
      # get marginal probability of exposure at different levels of exposure
      xmatpred <- xmat
      xmatpred[,xvar] <- levels(yvec)[1]
      p0=predict(rfMV,xmatpred,type="prob")[,class]
      xmatpred[,xvar] <- levels(yvec)[2]
      p1 <- predict(rfMV, xmatpred, type="prob")[,class]
      p0adj <- padjfromRF(rfMV, p0, class=class)
      p1adj <- padjfromRF(rfMV, p1, class=class)
  
      # calculate paired OR and get their quantiles      
      ORvec <- mapply(function(x,y) y/(1-y)/(x/(1-x)), p0adj, p1adj)
      kstest <- ks.test(log(ORvec), "pnorm")
        
      # consider +/-1.96 SE if normal
      if (kstest$p.value>0.1) {
        SE_beta <- sd(log(ORvec),na.rm=T) / sqrt(nrow(xmat))
        betaCI <- beta + c(-1.96,1.96)*SE_beta
      } else {
        betaCI <- quantile(log(ORvec), c(0.025,0.975), na.rm=T)
      }
      ORCI <- exp(betaCI)
      if (verbose==TRUE) {
        print(rfMV)
        hist(log(ORvec))  # log OR cos OR is non-negative
        summary(ORvec)
        qqnorm(log(ORvec))
        abline(1, 1)
        print(kstest)
        cat("beta:", beta, "[",betaCI,"]\n")
        cat("OR  :", OR, "[",ORCI,"]\n")
      }
      
    } else {
      stop("yvec must be a binary factor.")
    }  # end of if(binary classification) section
  }  # end of if(classification) section
  list(rfmod=rfMV, betaest=beta, betaCI=betaCI, ORest=OR, ORCI=ORCI,ORvec=ORvec)
}

# adapted from http://jayyonamine.com/?p=762
# https://gist.github.com/jayyonamine/5d9d802476a5c8e8b0cd#file-gistfile1-r

# library('ggplot2')
# library('randomForest')
# set.seed(2014)
# 
# rf_predict<-function(rf_object, data){
#   if (rf_object$type=="classification"){
#     p <-predict(rf_object, data, type="prob")
#     p<-as.vector(p[,2])
#   } else {
#     p <-predict(rf_object, data)
#   }
#   return (p)
# }
# 
# plot_partial<- 
#   function(rf, X, dv, iv, class, conf_int_lb=.025, 
#            conf_int_ub=.975, range_low=NULL, 
#            range_high=NULL, delta=FALSE, num_sample=NULL)
#   {
#     iv_name<-iv
#     if (is.factor(dv)==TRUE){
#       factor_var<-unique(dv)
#       #the test set needs all factor levels.  so, we build them and will drop them before we plot
#       factor_names=attributes(factor_var)$levels
#       fix_factor_df=data.frame(X[1:length(factor_names),])
#       colnames(fix_factor_df)=colnames(X)
#       fix_factor_df[, iv_name]=factor_names
#       y_hat_df=data.frame(matrix(NA,0, 2))
#       y_temp=data.frame(matrix(NA, nrow(X), 2))
#       y<-predict(rf, X)[,2]
#       for (i in 1:length(factor_names)){
#         X[, iv_name]=as.factor(factor_names[i])
#         X_temp=rbind(X, fix_factor_df)
#         p=as.vector(predict(rf, X, type="prob")[,class])
#         y_temp[,1]=p[1:nrow(X)] #drop the fix_factor_df rows
#         if (delta==TRUE){
#           y_temp[,1]<-y_temp[,1]-y
#         }
#         y_temp[,2]<-factor_names[i]
#         y_hat_df<-rbind(y_hat_df, y_temp)
#       }
#       plot<- qplot(y_hat_df[,2], y_hat_df[,1], 
#                    data = y_hat_df, 
#                    geom="boxplot", 
#                    main = paste("Partial Dependence of", (iv_name), "on", (dv_name))) +
#         ylab(bquote("Predicted values of" ~ .(dv_name))) +
#         xlab(iv_name)
#       return (plot)
#     } else {
#       conf_int <-(conf_int_ub-conf_int_lb)*100
#       temp<-sort(X[, iv_name])
#       if (is.null(num_sample)==FALSE){
#         temp<-sample(temp, num_sample)
#       }
#       if (is.null(range_low)==FALSE & is.null(range_high)==FALSE){
#         low_value<-quantile(temp, range_low)
#         high_value<-quantile(temp, range_high)
#         temp<-temp[temp<high_value & temp>low_value]
#       }
#       y_hat_mean<-vector()
#       y_hat_lb<-vector()
#       y_hat_ub<-vector()
#       y<-rf_predict(rf, X)
#       for (i in 1:length(temp)){
#         X[, iv_name] <- temp[i]
#         y_hat<-rf_predict(rf, X)
#         if (delta==TRUE){
#           y_hat<-y_hat-y
#         }
#         y_hat_mean[i]<-weighted.mean(y_hat)
#         y_hat_lb[i]<-quantile(y_hat, conf_int_lb)
#         y_hat_ub[i]<-quantile(y_hat, conf_int_ub)
#       }
#       df_new<-as.data.frame(cbind(temp, y_hat_mean, y_hat_lb, y_hat_ub))
#       plot<- ggplot(df_new, aes(temp)) + 
#         geom_line(aes(y=y_hat_mean), colour="blue") + 
#         geom_ribbon(aes(ymin=y_hat_lb, ymax=y_hat_ub), alpha=0.2) +
#         geom_rug(aes()) +
#         xlab(iv_name) +
#         ylab(bquote("Predicted values of" ~ .(dv_name))) +
#         ggtitle(paste("Partial Dependence of", (iv_name), "on", (dv_name), "\n with", (conf_int), "% Confidence Intervals"))
#       return (plot)
#     }
#   }
