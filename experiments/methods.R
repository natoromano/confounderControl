###############################################################################
# (Functions)
# Helper methods for various experiments. Most of them output plots.
#
# Nathanael Romano
###############################################################################
.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.1", .libPaths()))

RScriptPath="/home/yenlow/scripts"
RProjectPath="~/projects/confounderControl/"
source(paste(RScriptPath,"/R/utils.R", sep="")) # installnewpackage
source(paste(RScriptPath,"/R/logit.R", sep="")) # coeffAtlambda

installnewpackage(c("reshape2", "ggplot2", "gridExtra"))
require("reshape2")
require("ggplot2")


# OUTLIER ANALYSIS
# ----------------
# plot.PCA.outliers
# -----------------
# Given a dataset, and a set of errors, plots points from the dataset along
# with the color-coded errors (matched by columns), in a 2-dimensional space
# created using Principal Component Analysis.
# Covariates to be used for PCA can be specified as an extra argument, as
# well as a method name for the plot.
plot.PCA.outliers <- function(data, errors, pc1=1, pc2=2, covariates=NULL,
                              method='powers forest') {
  if (is.null(covariates)) covariates = setdiff(colnames(data), "id")
  pca <- prcomp(data[, covariates], scale=T)
  projected <- data.frame(pca$x)
  projected$error = errors
  PC1 = paste('PC', pc1, sep="")
  PC2 = paste('PC', pc2, sep="")
  pcplot <- qplot(x=projected[, PC1], y=projected[, PC2], 
                  colour=projected$error) +
    scale_colour_gradient(name="error (%)", low="green", high="red") +
    xlab(PC1) + ylab(PC2) + ggtitle('Relative error from') +
    ggtitle(paste('Error from', method))
  return(pcplot)
}


# plot.outliers
# -----------------
# Given a dataset, and a set of errors, plots points from the dataset along
# with the color-coded errors (matched by columns), using the two provided
# covariates.
plot.outliers <- function(data, errors, cov1, cov2) {
  pcplot <- qplot(x=data[, cov1], y=data[, cov2], 
                  colour=errors) +
    scale_colour_gradient(name="error", low="green", high="red") +
    xlab(cov1) + ylab(cov2) + ggtitle('Error on treatment effect')
  return(pcplot)
}


# plot.PCA.CI
# -----------------
# Given a dataset, and a set of effects, plots points from the dataset along
# with the color-coded CI (matched by columns), in a 2-dimensional space
# created using Principal Component Analysis. The column on which to find the
# the confidence intervals in the effects dataframe should be specified.
# Covariates to be used for PCA can be specified as an extra argument, as
# well as a method name for the plot.
plot.PCA.CI <- function(data, effects, pc1=1, pc2=2, covariates=NULL,
                        ci.column='propensity.ci') {
  if (is.null(covariates)) covariates = setdiff(colnames(data), "id")
  pca <- prcomp(data[, covariates], scale=T)
  projected <- data.frame(pca$x)
  projected$CI = effects[, ci.column] / sqrt(1000)
  PC1 = paste('PC', pc1, sep="")
  PC2 = paste('PC', pc2, sep="")
  pcplot <- qplot(x=projected[, PC1], y=projected[, PC2], 
                  colour=projected$CI) +
    scale_colour_gradient(name="95% CI width", low="green", high="red") +
    xlab(PC1) + ylab(PC2) + ggtitle('Width of 95% CI') +
    ggtitle(paste('95% CI width for', ci.column))
  return(pcplot)
}


# plot.power.CI
# -----------------
# Given a dataset, and a set of powers forest effects, plots the effects in
# sorted order along with the confience intervals (which shoudl be in the
# effects data frame as well)
plot.power.CI <- function(effects, data, order.by='powers', method="powers") {
  ci <- paste(method, 'ci', sep='.')
  df <- as.data.frame(effects)
  if (order.by %in% c('effect', 'propensity', 'powers')) {
    order <- order(effects[, order.by])
  } else {
    order <- order(data[, order.by])
  }
  df[, ci] = df[, ci] / sqrt(1000)
  df <- df[order,]
  df$id <- 1:nrow(df)
  p1 <- ggplot(df, aes(id, powers)) +
    geom_line() + 
    geom_ribbon(data=df, aes(ymin = powers-powers.ci, 
                             ymax = powers+powers.ci), alpha=0.3) + 
    xlab('Patient ID') + ylab('Estimated effect (+95% CI)') +
    ggtitle('Bias and 95% CI from power forest')
  return(p1)
}


# plot.propensity.CI
# -----------------
# Given a dataset, and a set of propensity forest effects, plots the effects in
# sorted order along with the confience intervals (which shoudl be in the
# effects data frame as well)
plot.propensity.CI <- function(effects, data, order.by='propensity', 
                               method="propensity") {
  ci <- paste(method, 'ci', sep='.')
  df <- as.data.frame(effects)
  if (order.by %in% c('effect', 'propensity', 'power')) {
    order <- order(effects[, order.by])
  } else {
    order <- order(data[, order.by])
  }
  df[, ci] = df[, ci] / sqrt(1000)
  df <- df[order,]
  df$id <- 1:nrow(df)
  p1 <- ggplot(df, aes(id, propensity)) +
    geom_line() + 
    geom_ribbon(data=df, aes(ymin = propensity-propensity.ci, 
                             ymax = propensity+propensity.ci), alpha=0.3) + 
    xlab('Patient ID') + ylab('Estimated effect (+95% CI)') +
    ggtitle('Bias and 95% CI from propensity forest')
  return(p1)
}


# plot.error.density
# ------------------
# Plots the density of errors on estimated effects, using R's default KDE.
plot.error.density <- function(effects, column) {
  effects <- as.data.frame(effects)
  p <- plot(density((effects[, column] - effects$effect)/effects$effect),
    main=paste('Relative biases from', column, 'forest'),
    xlab='Relative biases')
  return(p)
}


# plot.error.density
# ------------------
# Plots the density of estimated effects' CIs, using R's default KDE.
plot.CI.density <- function(effects, column) {
  ci <- paste(column, 'ci', sep='.')
  effects <- as.data.frame(effects)
  print(ci)
  p <- plot(density(effects[, ci] / sqrt(1000)),
            main=paste('95% CI widths from', column, 'forest'),
            xlab='Confidence interval width')
  return(p)
}


# apply.methods
# -------------
# Apply all previous methods and saves a pdf file with all plots, for a given
# set of effects and dataset.
apply.methods <- function(effects, data, name) {
  pdf(file=paste(name, 'densities.pdf', sep="_"))
  par(mfrow=c(2, 2))
  plot.CI.density(effects, 'powers')
  plot.CI.density(effects, 'propensity')
  plot.error.density(effects, 'powers')
  plot.error.density(effects, 'propensity')
  dev.off()
  pdf(file=paste(name, 'gg_CI.pdf', sep="_"))
  p1 <- plot.propensity.CI(effects, data)
  p2 <- plot.power.CI(effects, data)
  grid.arrange(p1, p2, ncol=1, nrow=2)
  dev.off()
  pdf(file=paste(name, 'gg_viz.pdf', sep="_"))
  p3 <- plot.PCA.CI(data, effects,
              covariates=c('gen', paste('x', 1:10, sep=''), 'exposure'),
              ci.column="powers.ci")
  p4 <- plot.PCA.CI(data, effects, 
              covariates=c('gen', paste('x', 1:10, sep=''), 'exposure'),
              ci.column="propensity.ci")
  powerError <- abs(effects[,'powers'] 
                          - effects[,'effect'])/effects[,'effect']
  propError <- abs(effects[,'propensity'] 
                   - effects[,'effect'])/effects[,'effect']
  p5 <- plot.PCA.outliers(data, powerError, method="powers forest",
                covariates=c('gen', paste('x', 1:10, sep=''), 'exposure'))
  p6 <- plot.PCA.outliers(data, propError, method="propensity forest",
                covariates=c('gen', paste('x', 1:10, sep=''), 'exposure'))
  grid.arrange(p3, p4, p5, p6, ncol=2, nrow=2)
  dev.off()
}


# HTN DATA
# --------
# Methods for applying the methods for HTN dataset. They may be obsolete, but
# can be used to create new methods for real clinical data.
plot.htn <- function(molten, method='density', outcome) {
  if (method == 'density') {
  htn.plot <- ggplot(data = molten, mapping = aes(effect, color = method)) +
    geom_density(alpha=0.1, adjust=1.,
                 trim=TRUE) + xlab("Effect") +
    theme_bw() + theme(legend.key = element_blank()) + 
    xlim(-10, 10) + coord_cartesian(ylim=c(0, 0.5))
    # geom_errorbar(aes(ymin=bias-se, ymax=bias+se), width=.9,
    #            position=position_dodge(0.05)) + 
    ggtitle(paste("Density of treatment effects for", outcome))
  } else {
    htn.plot <- ggplot(data = molten, mapping = aes(x=pid, y=effect, 
                                      group = method, color = method)) +
      geom_line() + xlab("Patient ID") +
      theme_bw() + theme(legend.key = element_blank()) +
      # geom_errorbar(aes(ymin=bias-se, ymax=bias+se), width=.9,
      #            position=position_dodge(0.05)) + 
      ggtitle(paste("Treatment effects for", outcome))
  }
  
  # save plot
  file <- paste('htn', method, outcome, 'results.png', sep="_")
  filename <- paste(RProjectPath, paste('HTN', file, sep="/"), sep="")
  ggsave(filename=filename, plot=htn.plot)
  return(htn.plot)
}

process.htn.results <- function(data, results, effects,
                                adj_methods, match_methods, outcome=1) {
  outcomes = c('Y.systolic', 'Y.diastolic')
  data = data.frame(pid=as.numeric(rownames(data)))
  for (m in adj_methods) {
    if (m %in% c('power', 'propensity')) {
      data[, m] <- effects[[outcome]][, m]
    } else {
      data[, paste(m, 'adj', sep='_')] <- rep(results[m, 'coeff_adj',
                                                outcomes[outcome]], nrow(data))
    }
  }
  for (m in match_methods) {
    if (m %in% c('power', 'propensity')) {
      data[, paste(m, 'adj', sep='_')] <- effects[[outcome]][, m]
    } else {
      data[, paste(m, 'match', sep='_')] <- rep(results[m, 'coeff_matched', 
                                                outcomes[outcome]], nrow(data))
    }
  }
  molten <- melt(data, id.vars='pid')
  colnames(molten) <- c('pid', 'method', 'effect')
  return(molten)
}


# BALANCE ASSESSMENT FUNCTIONS
# ----------------------------
# SMD.frame
# ---------
# Helper function to built a SMD plot. Creates an intermediate dataframe.
# SHOULD BE USED AFTER RUNNING METHODS, WHEN ENV VARIABLES STILL DEFINED
SMD.frame <- function(method) {
  results.after = afterMatching[[method]]
  smds.after = c()
  for (var in empvariables_num) {
    smds.after[var]<-abs(as.numeric(results.after[['numvar']][var, 'smd']))
  }
  for (var in empvariables_cat) {
    smds.after[var]<-abs(as.numeric(results.after[['catvar']][[var]][['smd']]))
  }
  df = data.frame(variable=names(smds.after), after=smds.after)
  colnames(df) <- c('variable', method)
  return(df)
}


# plot.SMD
# --------
# Given a list of methods, plots the SMDs for these methods, before and after
# matching
# SHOULD BE USED AFTER RUNNING METHODS, WHEN ENV VARIABLES STILL DEFINED
plotSMD <- function(methods) {
  # computing SMDs before matching
  results.before = beforeMatching
  smds.before = c()
  for (var in empvariables_num) {
    smds.before[var]<-abs(as.numeric(results.before[['numvar']][var, 'smd']))
  }
  for (var in empvariables_cat) {
    smds.before[var]<-abs(as.numeric(results.before[['catvar']][[var]][['smd']]))
  }
  dataplot <- data.frame(variable=names(smds.before), before=smds.before)
  
  # addind the methods in method
  for (method in methods) {
    dataplot <- merge(dataplot, SMD.frame(method), by='variable')
  }
  datamelt <- melt.data.frame(data=dataplot, id.vars=c('variable'), 
                              variable_name='Method',
                              value.name='SMD')
  varNames <- as.character(dataplot$variable)[order(dataplot$before)]
  datamelt$variable <- factor(datamelt$variable, levels=varNames)
  colnames(datamelt) <- c('variable', 'Method', 'SMD')
  ggplot(data = datamelt, mapping = aes(x = variable, y = SMD,
                                        group = Method, color = Method)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
    coord_flip() +
    theme_bw() + theme(legend.key = element_blank())
}


# plot2dPCA
# ---------
# Given a dataset and matchedIDs, performs a 2-d visualization of matched
# groups, using Principal Components Analysis.
plot2dPCA <- function(data, matchedID, method, pc1=1, pc2=2, cov=1:10) {
  if (method != 'before') {
      matchedIDs <- as.vector(matchedID[[1]][[method]])
      matcheddata = data[matchedIDs,]
    } else {
      matcheddata = data
    }
    covariates <- paste('x', cov, sep="")
    pca <- prcomp(matcheddata[, covariates], scale=T)
    projected <- data.frame(matcheddata$exposure, pca$x)
    PC1 = paste('PC', pc1, sep="")
    PC2 = paste('PC', pc2, sep="")
    pcplot <- qplot(x=projected[,PC1], y=projected[,PC2], 
                    colour=factor(matcheddata$exposure)) +
              scale_colour_hue(name="exposure") + xlab(PC1) + ylab(PC2) +
              ggtitle(paste('Matched data for', method)) + 
              xlim(-5, 5) + ylim(-5, 5) # for easy comparison, remove if nec.
    name = paste(paste(c('pc', pc1, pc2, method), collapse="_"), '.png', sep="")
    filename = paste(RProjectPath, paste('experiments/plots', name, sep="/"), 
                     sep="")
    ggsave(filename=filename, plot=pcplot)
    return(pcplot)
}

# resultsTables
# -------------
resultsTables <- function(resultsArray) {
  # store resultsArray separately to compute absolute deviations (old method)
  absolute.results <- resultsArray
  bias <- c('bias_matched', 'bias_adj')
  absolute.results[,bias,] <- abs(resultsArray[,bias,])
  
  matched.columns <- c('ORmatched', 'se_matched', 'bias_matched', 'n1')
  adj.columns <- c('ORadj', 'se_adj', 'bias_adj', 'n0')
  
  # compute results array
  apply(absolute.results[, matched.columns,], c(1, 2), mean, na.rm=T)
  apply(absolute.results[, adj.columns,], c(1, 2), mean, na.rm=T)
}
