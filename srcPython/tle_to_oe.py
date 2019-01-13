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

# This script reads a TLE file (with 1 or more TLE(s)), prints the inclination,
# eccentricity, RAAN, argument of perigee and mean anomaly directly from the TLE,
# and compute the sma from the # revolution/day (mean motion read from the TLE).
# Important note: these are not keplerian or osculating elements. The elements
# in the TLEs are different from the keplerian and osculating elements.
# Inputs:
# - tle_filename: name of the file that includes the TLE(s)
# - is_cygnss: if set to 1, assumes tle_filename has 8 TLEs from NORAD ID 41884,
# 41885, ..., to 41891 (in the ascending order of NORAD ID)

import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
import norad_id_to_cygnss_id; reload(norad_id_to_cygnss_id)
from norad_id_to_cygnss_id import *
from datetime import datetime, timedelta
import numpy as np

def get_date(read_tle):
    date_raw = read_tle[itle*2][18:32]
    yy = '20' + date_raw[:2]
    doy = 	date_raw[2:5]
    day_decimal = np.float('0.' + date_raw.split('.')[1])
    seconds = day_decimal * 24 * 3600
    date_str = yy + '-' + doy 
    date = datetime.strptime(date_str, "%Y-%j") + timedelta(seconds = seconds)
    date_str = datetime.strftime(date, "%Y-%m-%dT%H:%M:%S.%f")
    return date_str

tle_filename = 'tletemp.txt'
is_cygnss = 1

mu_earth = 398600.4418 # km^3/s^2
earth_radius = 6378.137 # mean equatorial radius (km)   
tle_file = tle_filename = open(tle_filename)
read_tle = tle_file.readlines()
ntle = len(read_tle) / 2
for itle in range(ntle):
    # NORAD ID
    norad = read_tle[itle*2][2:7]
    if is_cygnss == 1:
        cygnss = norad_id_to_cygnss_id(norad)
    
    # date
    date_str = get_date(read_tle)

    # inclination (degrees)
    inc = np.float(read_tle[itle*2+1][8:16])

    # RAAN (degrees)
    raan = np.float(read_tle[itle*2+1][17:25])

    # eccentricity
    ecc = np.float('0.' + read_tle[itle*2+1][26:33])

    # argument of perigee (degrees)
    arg_per = np.float(read_tle[itle*2+1][34:42])

    # mean anomaly (degrees)
    mean_ano = np.float(read_tle[itle*2+1][43:51])

    # mean motion (revs per day)
    mean_motion = np.float(read_tle[itle*2+1][52:63])

    # sma (km)
    sma = ( mu_earth / (mean_motion/24/3600*2*np.pi)**2 )**(1./3)
    sma_earth_radius = sma - earth_radius

    # print results
    if is_cygnss == 1:
        print 'CYGFM(NORAD)  inc     RAAN      ecc    arg_per mean_ano\
 mean_motion    SMA-Re'
        print cygnss + '(' + norad + ') ' + str(inc) + ' ' + str(raan) + ' ' + \
            str(ecc) + ' ' + str(arg_per) + ' ' + str(mean_ano) + ' ' + \
            str(mean_motion) + ' ' + str(sma_earth_radius) 
    else:
        print 'NORAD   inc     RAAN      ecc    arg_per mean_ano\
 mean_motion    SMA-Re'
        print  norad + ' ' + str(inc) + ' ' + str(raan) + ' ' + \
            str(ecc) + ' ' + str(arg_per) + ' ' + str(mean_ano) + ' ' + \
            str(mean_motion) + ' ' + str(sma_earth_radius) 
    
    print 'epoch: ' + date_str + '\n'
