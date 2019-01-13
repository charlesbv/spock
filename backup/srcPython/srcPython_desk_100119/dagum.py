import numpy as np

# PARAMETERS TO SET
delta = 0.05 # confidence level
pt = 0.00001   # Pt
epsilon_relative = 0.01  # relative error on Pt

# DAGUM'S EQUATION 
epsilon_absolute = epsilon_relative * pt # in Dagum's equation, the absolute error is used, which is the relative error times the probability
num_f1 = 4 * ( np.exp(1) - 2 )
num_f2 = ( 1 - pt ) * pt
num_f3 = np.log( 2 / delta )
num = num_f1 * num_f2 * num_f3
den = epsilon_absolute**2

nlim_samples = num / den 
nlim_ensembles = np.sqrt( nlim_samples )

print "- Error: " + '{0:.0e}'.format(epsilon_relative)
print "- Confidence level: " + '{0:.0f}'.format( ( 1 - delta ) * 100 ) + '%'
print "- Probability: ", pt
print "--> Minimum number of samples computed by Dagum's method: nlim_samples = " + '{0:.2e}'.format(nlim_samples) + "."
print "--> This corresponds to a minimum number of ensembles for each of the 2 spaceacreft of: nlim_ensembles = sqrt(nlim) = " + '{0:.2e}'.format(nlim_ensembles) + "."
