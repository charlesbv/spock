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
