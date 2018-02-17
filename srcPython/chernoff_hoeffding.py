import numpy as np

# PARAMETERS TO SET
delta = 0.05 # confidence level
epsilon_relative = 0.01 # relative error on Pt
pt = 0.216602000 # probability, does not directly appear in the final equation but used to compute the absolute error

# CHERNOOF-HOEFFDING'S EQUATION 
epsilon_absolute = epsilon_relative * pt # in the Chernoof-Hoeffding's equation, the absolute error is used, which is the relative error times the probability
num = np.log( 2 / delta )
den = 2 * epsilon_absolute**2
nlim_samples = num / den 
nlim_ensembles = np.sqrt( nlim_samples )

print "- Error: " + '{0:.0e}'.format(epsilon_relative)
print "- Confidence level: " + '{0:.0f}'.format( ( 1 - delta ) * 100 ) + '%'
print "- Probability: ", pt
print "--> Minimum number of samples computed by Chernoof-Hoeffding's method: nlim_samples = " + '{0:.2e}'.format(nlim_samples) + "."
print "--> This corresponds to a minimum number of ensembles for each of the 2 spaceacreft of: nlim_ensembles = sqrt(nlim) = " + '{0:.2e}'.format(nlim_ensembles) + "."
