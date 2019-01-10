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

# Converts ECI r, v to osculating elements at  given time
# Inputs:
# - r, v: postiion in ECI coordinates (in km, 3d np array)
# - time: UTC time. format YYYY-MM-DDTHH:MM:SS.microseconds or
# YYYY/MM/DD HH:MM:SS.microseconds (.microseconds is optional)
# - load_spice: load (1) or not (0) spice. load only if this
# script is called for the first time
# Output:orbitl elemnts in km, degrees
# TEST:
# rtest= np.array([-5516.4523254560, -3998.9893496060, 1032.9420683098])
# vtest = np.array([4.2320846839, -4.7300967309, 4.2046535015])
# timetest = '2018/10/23 00:00:00.000000'
# -> output:
# (6905.925593273137,
#  0.0026926214051435087,
#  34.92424067214202,
#  203.3990527168298,
#  53.452885630137864,
#  321.72524725164567,
#  5711.416379576319,
#  15.17813288178356)

import spiceypy as spice
import numpy as np
def cart2kep(r, v, time, load_spice):
    if load_spice == 1:
        path_spice = "/Users/cbv/cspice/data/"
        eop = path_spice + "pck00010.tpc"
        planet_ephem = path_spice + "de432s.bsp"
        earth_binary_pck = path_spice + "earth_000101_190217_181126.bpc"
        #earth_000101_181125_180903.bpc"

        leap_sec = path_spice + "naif0012.tls"
        spice.furnsh(eop)
        spice.furnsh(planet_ephem)
        spice.furnsh(earth_binary_pck)
        spice.furnsh(leap_sec)

    mu_earth = 398600.4418 # km^3/s^2
    rrss = np.linalg.norm(r)
    vrss = np.linalg.norm(v)
    unit_r = r / rrss
    rdotv = np.dot(r, v)

    # Solve for angular momentum (r x v)
    h = np.cross(r, v)
    mag_h = np.linalg.norm(h)

    # Solve for node vector
    K = np.zeros([3])
    K[2] = 1.0
    node = np.cross(K, h)
    mag_node = np.linalg.norm(node)
    unit_node = node / mag_node

    #  Solve for semi-major axis, sma
    sma = 1.0 / ( (2.0/rrss) - ( (vrss*vrss)/mu_earth ) )
    
  
    #  Solve for ecentricity, e and the e vector
    specific_energy = -1.0*(mu_earth/(2.0*sma))
    tempv1 = v * rdotv
    coeff = vrss*vrss - (mu_earth/rrss)

    tempv2 = r*coeff
    e_vector = tempv2 - tempv1

    coeff = 1.0/mu_earth
    e_vector = e_vector * coeff

    eccentricity = np.linalg.norm(e_vector)
    unit_e = e_vector / eccentricity

    #  Solve for inclination, i
    inclination = np.arccos(h[2]/mag_h)

    # Solve for longitude of ascending node
    if (mag_node == 0.0):
        long_an = 0.0
    elif (node[1] >= 0):
        long_an = np.arccos(node[0]/mag_node)
    elif (node[1] < 0):
        long_an = 2*np.pi - np.arccos(node[0]/mag_node)
        
    # Solve for argument of periapse  
    if (mag_node != 0.0):
        coeff = np.dot(unit_node, unit_e)
        w = np.arccos( coeff )
        if (e_vector[2] < 0.0):
            w = (2.0*np.pi - w)
        # else:
        #     w = 0
  

    # Solve for true anomaly
    e_vector_dot_r = np.dot(e_vector, r)
    if (eccentricity != 0):
        f = np.arccos(e_vector_dot_r/(eccentricity*rrss))
        if (rdotv < 0):
          f = (2*np.pi) - np.abs(f)       
    else:
        f = 0
  
    # orbital period
    period = sma**3.0
    period = period / mu_earth
    period = 2.0 * np.pi * np.sqrt( period )

    # AN to sc (phase_angle)
    phase_angle = np.mod(w + f, 2*np.pi)
    rad2deg = 180./np.pi
    return sma, eccentricity, inclination*rad2deg, long_an*rad2deg,\
        w*rad2deg, f*rad2deg, period, phase_angle*rad2deg

