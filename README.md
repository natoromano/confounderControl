# confounderControl
This code was adapated from Yen Low's code for the paper "Comparing high-dimensional confounder control methods for rapid cohort studies from electronic health records" (published in the Journal of Comparative Effectiveness Research).
However it has been widely refactored and updated (to add handling of continous outcomes, and diagnostics methods, increasing confounding strength in simulations, and causal forests). 


# Organisation
Its main parts are the following.

OLD SIMULATIONS
- the data directory reproduces the simulations from the initial paper, with a few modifications (removing non-linearities in the PS model, and sampled signs for the inter-covariates correlations)
- the model/sim2000x10 directory implements 17 methods (in ps.R) that can be applied to the previous simulation

EXPERIMENTS
- methods.R has various methods that can be applied on the results of experiments, to produce informative plots on covariate balance, lower-dimensional representation of errors / confidence intervals, density plots, etc.
- confounding.R implements an experiment where the previous simulation is ran for increasing confounding strength, and the evolution of bias with confounding strength is plotted for various methods

HTN DATA
- the HTN directory has a pre-processed version of the hypertension (HTN) dataset, and the code to apply the 17 methods to it (in ps.R). The main differences with the previous simulations are the absence of ground truth effect (which is probably heterogeneous), and the fact that the outcome is continous (there are in fact two outcomes: decrease in diastolic and systolic blood pressure).

SCOTT POWERS'S SIMULATION
- the scott directory reproduces the simulations built by Scott Powers to demonstrate the effectiveness of causal forests, and applies multiple of the previous methods to them.

SEMI-CLINICAL SIMULATION
- the realistic directory holds an experiment where a causal structure is built on top of actual EHR data (sample from the HTN dataset), to test various treatment effect estimation methods
- the file sim.R holds the actual simulation and methods, whereas ehr-features.R offers the same simulation as before, where the confounding strength is increased to test the methods' robustness to high confounding
- this directory also holds a source of EHR covariates which can be used for simulations

# Nathanael Romano, Yen Low, Stanford Center for Biomedical Informatics Research

