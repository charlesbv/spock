# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

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
