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
# This script calculates the accommodation coefficient of a surface using Goodman (1976) model. The equation is found in moe05 (eq 5)
# !!!!!!!!! figure 4 of moe05 shows  a bad agreement with orbit data < 300 km. they expect better agreeement at higher altitudes but has this agreement been actually verified?...

import numpy as np

mass_incoming_molecule  = 15. # in amu (An atomic mass unit (symbolized AMU or amu) is defined as precisely 1/12 the mass of an atom of carbon-12)
mass_surface_molecule =  56.# see array_mass_surface_molecule in amu (An atomic mass unit (symbolized AMU or amu) is defined as precisely 1/12 the mass of an atom of carbon-12)
u = mass_incoming_molecule / mass_surface_molecule
factor = 3.6 * u  / ( 1 + u )**2
print 'factor', factor
# theta = # angle between the incident velocity vector and the TANGENT to the surface.
# alpha = factor * np.sin(theta) 
