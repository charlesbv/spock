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

# plot the lon lat from the files in toshare

import numpy as np
import matplotlib
import pickle
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import pickle
import fileinput
import time
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import subprocess
from matplotlib.ticker import FixedLocator
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap, shiftgrid

filename = 'toshare/20180910T000000_to_20180911T000000.txt'

filespec = open(filename)
read_filespec = filespec.readlines()
nheader = 1
n = len(read_filespec) - nheader
lon = []
lat = []
gain = []
date = []
for i in range(n):
    lon_sub = []
    lat_sub = []
    gain_sub = []
    date_temp  = read_filespec[i + nheader].split()[0]
    date_temp_date = datetime.strptime(date_temp, "%Y-%m-%dT%H:%M:%S")
    date.append(date_temp_date)
    nb_spec_now = (int)((len(read_filespec[i + nheader].split()) - 1) / 3)
    color_spec_sub = ['','','','']
    for ispec in range(nb_spec_now):
        lon_sub.append( np.float( read_filespec[i + nheader].split()[1 + 3*ispec] ) )
        lat_sub.append( np.float( read_filespec[i + nheader].split()[2 + 3*ispec] ) )
        gain_sub.append( np.float( read_filespec[i + nheader].split()[3 + 3*ispec] ) )
    lon.append(lon_sub)
    lat.append(lat_sub)
    gain.append(gain_sub)
            

date = np.array(date)
date_start = date[0]
n_first_hours_val = 10 # in min
n_first_hours_val_sec = n_first_hours_val * 60.
index_n_first_hours = np.where(date < date[0] + timedelta(seconds = n_first_hours_val_sec))[0]
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 16
ratio_fig_size = 4./3
width_fig = 25

dindex = 10
icount = -1

min_lat_range = 0.
max_lat_range = 40.
min_lon_range = -98.
max_lon_range = -45.
step_lon = ( (max_lon_range - min_lon_range)  / 10. )
step_lat =  ( (max_lat_range - min_lat_range)/ 10. )
array_lon_temp = np.arange(min_lon_range, max_lon_range + step_lon,step_lon) 
array_lat_temp = np.arange(min_lat_range, max_lat_range + step_lat,step_lat) 
array_lon = []
array_lat = []
for i in range(len(array_lon_temp)):
    if (array_lon_temp[i] < 0):
        array_lon.append( format(np.abs(array_lon_temp[i]), ".1f") + 'W' )
    else:
        array_lon.append( format(array_lon_temp[i], ".1f") + 'S' )
for i in range(len(array_lat_temp)):
    if (array_lat_temp[i] < 0):
        array_lat.append( format(np.abs(array_lat_temp[i]), ".1f") + 'S' )
    else:
        array_lat.append( format(array_lat_temp[i], ".1f") + 'N' )


for iii in np.arange(0,n,10):#np.arange(index_n_first_hours[0], index_n_first_hours[-1], dindex):
    icount = icount + 1
    print date[iii]
    fig_title = ''#Probability per $P_C$ bin'
    fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in normal
    gs = gridspec.GridSpec(1, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

    ax = fig.add_subplot(gs[0, 0])


    ax.xaxis.set_major_locator(FixedLocator(np.arange(min_lon_range, max_lon_range+1, step_lon)))
    ax.yaxis.set_major_locator(FixedLocator(np.arange(min_lat_range, max_lat_range+1, step_lat)))

    ax.set_xticklabels(array_lon, color = 'k', fontsize = fontsize_plot)
    ax.set_yticklabels(array_lat, color = 'k', fontsize = fontsize_plot)



    plotMap = Basemap( projection       = 'cyl',
                            llcrnrlon        = min_lon_range,#-180,#min_lon , #Lower Left  CoRNeR Longitude
                            urcrnrlon        = max_lon_range  ,#180,#max_lon  , #Upper Right CoRNeR Longitude
                            llcrnrlat        = min_lat_range,#-40,#,min_lat  , #Lower Left  CoRNeR Latitude
                            urcrnrlat        = max_lat_range,#40,#max_lat,   #Upper Right CoRNeR Latitude
             resolution       = 'l'  ,
             suppress_ticks   = False,
             ax = ax,
             )

    alpha = 0.5
    color_continents = [65,105,225,alpha*256]
    color_continents = np.array(color_continents) / 256.
    color_water  = [100,149,237,alpha*256]
    color_water = np.array(color_water) / 256.

    plotMap.fillcontinents(color=tuple(color_continents),lake_color=tuple(color_water))
    plotMap.drawmapboundary(fill_color=tuple(color_water))

    meridians = plotMap.drawmeridians(np.arange(min_lon_range, max_lon_range,step_lon))
    parallels = plotMap.drawparallels(np.arange(min_lat_range, max_lat_range,step_lat))

    #mapRef = plotMap.drawcoastlines(linewidth=0.7, color='darkblue')

    y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
    x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'        
    ax_title = str(date[iii]) + ' UTC'
    ax.set_title(ax_title, weight = 'normal', fontsize  = (int)(fontsize_plot*1.1), y = 1.001)
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='normal') ## make the labels of the ticks in normal   

    nb_spec_now = len(lon[iii])
    for ispec in range(nb_spec_now):
        x, y =  plotMap(lon[iii][ispec], lat[iii][ispec])
        print lon[iii][ispec], lat[iii][ispec]
        point_plot = plotMap.scatter(x, y,  color = 'red', s = 20, zorder = 5)

    fig_save_name = 'ani/' + filename.split('/')[-1].replace(".txt", '_lon_lat_' + str(icount) + '.png')
    fig.set_figheight(height_fig)
    fig.set_figwidth(width_fig)

    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    
os.system('ffmpeg -y -r 10 -i ani/' + filename.split('/')[-1].replace(".txt", '_lon_lat_%d.png') + ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ani/' +  filename.split('/')[-1].replace(".txt", '_lon_lat.mp4'))
