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

# This script implements equation (3) in Moe04

import numpy as np
import math

RAD2DEG = 180. / np.pi
DEG2RAD = np.pi / 180.
earth_radius        = 6378.137; # mean equatorial radius (km)
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)


# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
alpha = 0.95 # accomodation coefficient # 0.93 in sutton09a
Tw = 300. # temperature of the satellite surface # assumed at 300 K in Moe04, 273 in Sutton07 and Bruisma03
h_sat = 350. # satellite altitude in km
v_sat = np.sqrt( earth_mu / ( earth_radius + h_sat ) ) * 1000. # np.array([7.7,0 ,0]) * 1000. # in m/s should be a vector but ok if it's a magnitude
Ta = 500 # atmospheric temperature in K
A = np.array([0.1, 0.1, 0.2, 0.2, 0.2, 0.2]) # ram, wake, zenith, nadir, starboard, port. In m^2
flying = 'foreward' # sideways, foreward, hyperthermal
if flying == 'hyperthermal': # Ta -> 0 so that s->+infinity. also assume only one plate
    #Ta = 0
    A = np.array([0.1])
    angle_normal_to_v_sat = np.array([0.]) * DEG2RAD # np.arccos( v_sat_norm_dot_normal )
elif flying == 'foreward':
    angle_normal_to_v_sat = np.array([0., 180., 90., 90., 90., 90.]) * DEG2RAD # np.arccos( v_sat_norm_dot_normal )
elif flying == 'sideways':
    angle_normal_to_v_sat = np.array([90., 90., 90., 90., 0., 180.]) * DEG2RAD # np.arccos( v_sat_norm_dot_normal )
else: 
    print "***! Flying needs to be 'foreward' or 'sideways' or 'hyperthermal'. The program will stop !***"; raise Exception
gamma = np.cos( angle_normal_to_v_sat )
nb_surf = len(A)
A_ref =  A * gamma # 0.2*np.ones([nb_surf]) # in m^2
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT


v_sat_norm = np.linalg.norm(v_sat)
# v_sat_normalized = v_sat / v_sat_norm
# normal = []
# v_sat_norm_dot_normal = np.dot( v_sat_normalized, normal )
R = 8.31 # universal gas constant in J/mol.K
Ma = 18. / 1000 # in kg/mol
s = v_sat_norm / np.sqrt ( ( 2 * R * Ta / Ma) )
Q = 1 + 1 / ( 2 * s**2 )
T = Ma * v_sat_norm**2 / ( 3 * R )  # kinetic temperature
Vr = v_sat * np.sqrt( 2./3 * ( 1 + alpha * ( Tw / T - 1 )  ) )
Vr_norm = np.linalg.norm(Vr)
cd = []
for isurf in range(nb_surf):
    P = np.exp( -( gamma[isurf]**2 ) * ( s**2 ) ) / s
    Z = 1 + math.erf( gamma[isurf] * s )

    term1 = P / np.sqrt(np.pi)
    term2 = gamma[isurf] * Q * Z
    fac1 = gamma[isurf] * Vr_norm / (2 * v_sat_norm)
    fac2 = gamma[isurf] * np.sqrt(np.pi) * Z + P
    term3 = fac1 * fac2
    term = term1 + term2 + term3
    cd.append( A[isurf]/A_ref[isurf]*term )

cd = np.array(cd)
#print cd


# TOTAL NORMALIZED DRAG COEFFICIENT
A_ref_tot = 0
for isurf in range(nb_surf):
    if gamma[isurf] >= 0:
        A_ref_tot = A_ref_tot + A_ref[isurf]
cd_tot_norm = np.sum(cd * A_ref) / A_ref_tot # see sutton09a equation 9
print angle_normal_to_v_sat * RAD2DEG
print cd
print cd * A_ref
print "XXXXXXXXXXXX"
print cd_tot_norm
print A_ref_tot
print "XXXXXXXXXXXX"
