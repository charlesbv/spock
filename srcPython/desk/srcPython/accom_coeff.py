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
