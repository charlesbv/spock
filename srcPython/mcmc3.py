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
import pymc as pm
from matplotlib import pyplot as plt

#import MyModule
#M = Model(MyModule)
    


# def make_model(x):
#     a = pm.Exponential('a',beta=x,value=0.5)

#     @pm.deterministic
#     def b(a=a):
#         return 100-a

#     @pm.stochastic
#     def c(value=0.5, a=a, b=b):
#         return (value-a)**2/b

#     return locals()

# M = pm.Model(make_model(3))
    


# from pymc.examples import gelman_bioassay
# M = pm.MAP(gelman_bioassay)
# M.fit()



# M.alpha.value
# #array(0.8465892309923545)
# M.beta.value
# #array(7.7488499785334168)



# M.AIC
# #7.9648372671389458
# M.BIC
# #6.7374259893787265



# N = pm.NormApprox(gelman_bioassay)
# N.fit()



# N.mu[N.alpha]
# #array([ 0.84658923])
# N.mu[N.alpha, N.beta]
# #array([ 0.84658923,  7.74884998])
# N.C[N.alpha]
# #matrix([[ 1.03854093]])
# N.C[N.alpha, N.beta]
# #matrix([[  1.03854093,   3.54601911],
# #        [  3.54601911,  23.74406919]])



# N.sample(100)
# N.trace('alpha')[::10]
# #array([-0.85001278,  1.58982854,  1.0388088 ,  0.07626688,  1.15359581,
# #       -0.25211939,  1.39264616,  0.22551586,  2.69729987,  1.21722872])
# N.trace('beta')[::10]
# #array([  2.50203663,  14.73815047,  11.32166303,   0.43115426,
# #        10.1182532 ,   7.4063525 ,  11.58584317,   8.99331152,
# #        11.04720439,   9.5084239 ])


### SETUP ###
x = pm.Normal('x', 0,1)
y = pm.Normal('y', 3,.2)
z = pm.Normal('z', x, y, value=[-1,2,4], observed=True)
M = pm.MCMC([x,y,z])
###

M.use_step_method(pm.Metropolis, x, proposal_sd=1., proposal_distribution='Normal')


"""
   if all(self.stochastic.value != 0.):
       self.proposal_sd = ones(shape(self.stochastic.value)) * \
                           abs(self.stochastic.value) * scale
   else:
       self.proposal_sd = ones(shape(self.stochastic.value)) * scale
"""   



M.use_step_method(pm.AdaptiveMetropolis, [x,y,z], \
                      scales={x:1, y:2, z:.5}, delay=10000)



A = pm.Normal('A', value=np.zeros(100), mu=0., tau=1.)



A = [pm.Normal('A_%i'%i, value=0., mu=0., tau=1.) for i in xrange(100)]
