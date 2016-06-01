# confounderControl
This code was adapated from Yen Low's code for the paper "Comparing high-dimensional confounder control methods for rapid cohort studies from electronic health records" (published in the Journal of Comparative Effectiveness Research).
However it has been widely refactored and updated (to add handling of continous outcomes, and diagnostics methods, increasing confounding strength in simulations, and causal forests). 

# Organisation
Its main parts are the following.

### OLD SIMULATIONS
The *data* directory reproduces the simulations from the initial paper, with a few modifications (removing non-linearities in the PS model, and sampled signs for the inter-covariates correlations). The *model/sim2000x10* directory implements 17 methods (in ps.R) that can be applied to the previous simulation.

### EXPERIMENTS
In the *experiments* directory, methods.R has various methods that can be applied on the results of experiments, to produce informative plots on covariate balance, lower-dimensional representation of errors / confidence intervals, density plots, etc.
confounding.R implements an experiment where the previous simulation is ran for increasing confounding strength, and the evolution of bias with confounding strength is plotted for various methods.

### HTN DATA
*Note: the data is not on the online repository.*
The *HTN* directory has a pre-processed version of the hypertension (HTN) dataset, and the code to apply the 17 methods to it (in ps.R). The main differences with the previous simulations are the absence of ground truth effect (which is probably heterogeneous), and the fact that the outcome is continous (there are in fact two outcomes: decrease in diastolic and systolic blood pressure).

### SCOTT POWERS'S SIMULATIONS
The *scott* directory reproduces the simulations built by Scott Powers to demonstrate the effectiveness of causal forests, and applies multiple of the previous methods to them.

### SEMI-CLINICAL SIMULATION
*Note: the source of EHR covariates is not on the online repository.*
The *realistic* directory holds experiments where a causal structure is built on top of actual EHR data (samples from the HTN dataset), to test various treatment effect estimation methods. 
The file sim.R holds the actual simulation and methods, whereas ehr-features.R offers the same simulation as before, where the confounding strength is increased to test the methods' robustness to high confounding. This directory also holds a source of EHR covariates which can be used for simulations

# References
Low Y.S., Gallego B., Shah N.H., Comparing high-dimensional confounder control methods for rapid cohort studies from electronic health records.. J Comp Eff Res, 2015.

Rosenbaum P.R., Rubin D.B., The Central Role of the Propensity Score in Observational Studies for Causal Effects. Biometrika, 1983.

Austin P.C., An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behav Research, 2011.

Stuart E.A., Matching methods for causal inference: A review and a look forward. Stat Sci, 2010.

Austin P.C., Type I Error Rates, Coverage of Confidence Intervals, and Variance Estimation in Propensity-Score Matched Analyses. Int J Biostat, 2009.

Rosenbaum P.R., Rubin D.B., The Bias Due to Incomplete Matching. Biometrics, 1985.

Imbens G. W., Nonparametric Estimation of the Average Treatment Effects Under Exonogeneity: A Review. The Review of Economics and Statistics, 2004.

Imbens G.W., Abadie A., Bias-Corrected Matching Estimators for Average Treatment Effects. Journal of Business Economics and Statistics, 2011.

Hade M. E., Bo L., Bias associated with using the estimated propensity score as a regression covariate. Stat Med., 2014.

Lunt M., Selecting an Appropriate Caliper Can Be Essential for Achieving Good Balance With Propensity Score Matching., Am. J. Epidemiol, 2014.

Pike M.C., Anderson J., Day N., Some insights into Miettinen's multivariate confounder score approach to case-control study analysis., Epidemiology and Community Health, 1979.

Austin P.C., Balance diagnostics for comparing the distribution of baseline covariates between treatment groups in propensity-score matched samples., Stat. Med., 2009.

Wooldridge J., Should Instrumental Variables be Used as Matching Variables?, Michigan State University, MI, 2009.

Webb M., Wilson J., Chong J., An Analysis of Quasi-complete Binary Data with Logistic Models:  Applications to Alcohol Abuse Data, Journal of Data Science, 2004.

Athey S., Imbens G., Machine Learning Methods for Estimating Heterogeneous Causal Effects, arXiv, 2015.

Wager S., Athey S., Estimation and Inference of Heterogeneous Treatment Effects using Random Forests, arXiv, 2014.

Wager S., Hastie T., Efron B.,  Confidence Intervals for Random Forests: The Jackknife and the Infinitesimal Jackknife, arXiv, 2014.

#### Nathanael Romano, Yen Low, Stanford Center for Biomedical Informatics Research

