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
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import time
from datetime import datetime
from matplotlib import pyplot as plt
import os



# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
i = 97 # inclination of orbit
hp = 500. # altitude of perigee
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

# Constants
re        = 6378.137; # mean equatorial radius (km)
mu    = 398600.4418; # gravitational parameter (km^3/s^2)
j2    = 1.081874e-3; # J2 zonal harmonic coefficient (no unit)
rad2deg = 180./np.pi
deg2rad = 1 / rad2deg
second2year = 3600 * 24 * 365.25
second2day = 3600 * 24

# apsidal precession rate
rp = re + hp
ra_min = 1 * re
ra_max = 10 * re
dra = 0.2 * re
ra_arr = np.arange(ra_min, ra_max + dra, dra)
ra_arr_re = ra_arr/re
nb_ra  = len(ra_arr)
dwdt = np.zeros([nb_ra])
for ira in range(nb_ra):
    ra = ra_arr[ira]
    sma = ( ra + rp ) / 2.
    e = ( ra - rp ) / ( ra + rp )
    T = 2 * np.pi * np.sqrt( sma**3 / mu )
    n0 = 2 * np.pi / T
    k = re**2 / ( sma**2 * ( 1 - e**2 )**2 )
    dwdt[ira] = 3./4 * n0 * k * j2 * ( 5 * (np.cos(i*deg2rad))**2  - 1 )

# Plot
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

fig_title = 'Rate of apsidal precession VS radius of apogee'
y_label = 'Absolute rate (' + u'\N{DEGREE SIGN}/2 yr)'
x_label = '$r_a$ (units of $R_E$)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

x_axis = ra_arr_re
y_axis = np.abs(dwdt * rad2deg * second2year * 2) # * 2 because want degree per 2 years

ax.plot(x_axis, y_axis, linewidth = 2, color = 'k', label = '')

where_dwdt_below_90 = np.where(np.abs(dwdt * rad2deg * second2year * 2) < 90)[0][0]
ra_dwdt_below_90 = ra_arr_re[where_dwdt_below_90]
ax.plot([min(x_axis), max(x_axis)], [90,90], linewidth = 2, color = 'r', label = '')
ax.plot([ra_dwdt_below_90, ra_dwdt_below_90], [0, np.abs(dwdt[where_dwdt_below_90] * rad2deg * second2year * 2)], linewidth = 2, color = 'b', label = '')
ax.margins(0,0)
ymin = 0 # min of y axis
#ax.set_ylim([ymin, np.max(y_axis)])
ax.set_ylim([0,200])
ax.text( (max(x_axis) + min(x_axis))/2, 90 - ( max(y_axis) - ymin)/250., '90' + u'\N{DEGREE SIGN}', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold', color = 'r')
ax.text( ra_dwdt_below_90, ymin  - ( max(y_axis) - ymin)/38., format(ra_dwdt_below_90, ".1f"), horizontalalignment = 'center', verticalalignment = 'top', fontsize = fontsize_plot, weight = 'bold', color = 'b')



fig_save_name = 'apsidal_precession.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

