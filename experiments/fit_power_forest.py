#! /usr/bin/env python

################################################################################
# Usage: fit_powers_forest.py <training file> <test file> <output file> 
# <n_estimators> <output: effect or outcomes>
# Fit a Powers forest to the data in <training file>, which has tab
# delimited columns, no header row, and columns <tx>\t<outcome><features...>
# Then apply the model to the data in <test file> (same format).
# Output the estimates on the test data to <output file>
################################################################################

import sys
import numpy as np
from sklearn.ensemble.forest import PowersForest

training_file = sys.argv[1]
test_file = sys.argv[2]
output_file = sys.argv[3]
n_estimators = int(sys.argv[4])
confidence_interval = sys.argv[6] != 'no'

training_data = np.loadtxt(fname=training_file)
train_w = training_data[:,0]
train_y = training_data[:,1]
train_X = training_data[:,2:]

test_data = np.loadtxt(fname=test_file)
test_w = test_data[:,0]
test_y = test_data[:,1]
test_X = test_data[:,2:]

pf = PowersForest(n_estimators=n_estimators,
                  min_samples_leaf=3)
print "Fitting Powers Trees..."
pf.fit(X=train_X, y=train_y, w=train_w)


print "Applying to test set..."
method = sys.argv[5]
if method == 'effect':
    # Returns vector of mean((treated outcomes) - (control outcomes))
    if confidence_interval:
        estimates, ci = pf.predict_effect(X=test_X, 
                        confidence_interval=True)
    else:
        estimates = pf.predict_effect(X=test_X, 
                        confidence_interval=False)
else:
    # Returns matrix of [mean(treated outcomes), mean(control outcomes)]
    estimates = pf.predict_outcomes(X=test_X)

np.savetxt(fname=output_file,
           X=estimates,
           delimiter='\t')

if confidence_interval:
    np.savetxt(fname=sys.argv[6],
               X=ci,
               delimiter='\t')

exit
