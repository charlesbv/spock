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

# This script is used to post process the data for the analysis with Chris and Scott about sampling with all GNSS systems (GPS, Galileo, Glonass, etc)

from datetime import datetime, timedelta
import numpy as np
import os
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
#sys.path.append("/home/cbv/spock_development_new_structure_kalman_dev/srcPython")
from os import listdir
from read_input_file import *
from read_output_file import *
from gnss_read_spock_spec_bin import *
from netCDF4 import Dataset
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
#from ecef2eci import *
#from eci_to_lvlh import *
from ecef_to_lvlh import *
import pickle
from cygnss_name_to_norad_id import *
import os.path
import pdb
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


input_filename = 'sp1w_gps.txt'
os.system("mpirun -np 4 spock_dev " + input_filename)
os.system("mpirun -np 4 spec_gnss " + input_filename + " -lon=0 -rot=0 -min")
print 'Post processing...'
grid_spec_sat_list = []
grid_gnss_sat_list = []
for cygfm in range(1, 9):
#cygfm = 1
    print 'FM0' + str(cygfm)
    var_in, var_in_order = read_input_file(input_filename)
    dt_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
    gnss_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
    nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')];
    cygfm_to_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
    isc =  cygfm_to_nb[cygfm-1] - 1
    spec_filename = output_file_path_list[isc] + "specular_" + output_file_name_list[isc].replace(".txt",".bin")

    lonGridInfo = [0., 360., 1.] # min, max, size
    latGridInfo = [-45., 45., 1.]
    data_spec = gnss_read_spock_spec_bin(spec_filename.replace('.txt','.bin'), gnss_name_list, dt_output, lonGridInfo, latGridInfo, 0) 
    date = data_spec[0]; lon = data_spec[1]; lat = data_spec[2]; gain = data_spec[3]; gnss = data_spec[4]; grid_spec_sat = data_spec[5]; grid_gnss_sat = data_spec[6]
    if cygfm == 1:
        grid_spec = np.array(grid_spec_sat)
    else:
        grid_spec = grid_spec + np.array(grid_spec_sat)
    if cygfm == 1:
        grid_gnss = np.array(grid_gnss_sat)
    else:
        grid_gnss = grid_gnss + np.array(grid_gnss_sat)
    grid_spec_sat_list.append(grid_spec_sat)
    grid_gnss_sat_list.append(grid_gnss_sat)

lonFilledInd = np.where(grid_spec > 0)[1]
latFilledInd = np.where(grid_spec > 0)[0]


 
minLonGrid = lonGridInfo[0]
maxLonGrid = lonGridInfo[1]
dlon = lonGridInfo[2] # in degree
minLatGrid = latGridInfo[0]
maxLatGrid = latGridInfo[1]
dlat = latGridInfo[2] # in degree
nlon = (int)((maxLonGrid - minLonGrid)/dlon)
nlat = (int)((maxLatGrid - minLatGrid)/dlat)
lonGrid = np.arange(minLonGrid, maxLonGrid, dlon)
latGrid = np.arange(minLatGrid, maxLatGrid, dlat)

lonFilled = lonGrid[lonFilledInd]
latFilled = latGrid[latFilledInd]


# PLOTS

height_fig = 18.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Latitude'
x_label = 'Longitude'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')


fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(2, 1, width_ratios=[1],
                          height_ratios=[ 0.1,1])

#                  gridspec_kw={"height_ratios":[1, 0.05]})
ax = fig.add_subplot(gs[1, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           


# for i in range(len(lonGrid)):
#     grid[3, i] = 100.


#cmap = colors.ListedColormap(['azure', 'lightblue', 'blue', 'red'])  

m = Basemap( projection       = 'cyl',
             llcrnrlon        =  lonGridInfo[0], #Lower Left  CoRNeR Longitude
             urcrnrlon        =  lonGridInfo[1] , #Upper Right CoRNeR Longitude
             llcrnrlat        = latGridInfo[0]  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = latGridInfo[1],   #Upper Right CoRNeR Latitude
             resolution       = 'c'  , #'c' 'l' 'i', 'h' 'f'
             suppress_ticks   = False,
             ax = ax,
             )
m.drawcoastlines(linewidth=0.7, color='blue')


levels = MaxNLocator(nbins=15).tick_values(grid_gnss.min(), grid_gnss.max())
cmap = plt.get_cmap('jet')
# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels. (https://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)


lon_arr = np.zeros([nlon+1])
for ilon in range(nlon+1):
    lon_arr[ilon]  = minLonGrid + ilon * dlon
lat_arr = np.zeros([nlat+1])
for ilat in range(nlat+1):
    lat_arr[ilat]  = minLatGrid + ilat * dlat


im = m.pcolormesh(lon_arr, lat_arr, grid_gnss, cmap=cmap, norm=norm)#linewidth = 0.01,  edgecolors = 'lightgrey')  # !!!!!! be careful: X and Y need to be the number of element of data + 1 (),  otherwise one cell
                                                                                              # will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)
im.set_rasterized(True) # to hide the edges of the bins in the grid
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

# cax = fig.add_subplot(gs[0, 0])
# #cax.axis('off')
# fig.colorbar(im, ax=cax,  orientation='horizontal')
# gs.update(left = 0.11, right=0.87, top = 1,bottom = 0.12, hspace = -0.3)

#ax.margins(0,0)
fig_save_name = 'grid_gnss_' + input_filename.replace('.txt', '.pdf')
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')


height_fig = 18.  # the width is calculated as height_fig * 4/3.                                                                             
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'Latitude'
x_label = 'Longitude'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')


fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(2, 1, width_ratios=[1],
                          height_ratios=[ 0.1,1])

#                  gridspec_kw={"height_ratios":[1, 0.05]})
ax = fig.add_subplot(gs[1, 0])
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           


# for i in range(len(lonGrid)):
#     grid[3, i] = 100.


#cmap = colors.ListedColormap(['azure', 'lightblue', 'blue', 'red'])  

m = Basemap( projection       = 'cyl',
             llcrnrlon        =  lonGridInfo[0], #Lower Left  CoRNeR Longitude
             urcrnrlon        =  lonGridInfo[1] , #Upper Right CoRNeR Longitude
             llcrnrlat        = latGridInfo[0]  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = latGridInfo[1],   #Upper Right CoRNeR Latitude
             resolution       = 'c'  , #'c' 'l' 'i', 'h' 'f'
             suppress_ticks   = False,
             ax = ax,
             )
m.drawcoastlines(linewidth=0.7, color='blue')


levels = MaxNLocator(nbins=15).tick_values(grid_spec.min(), grid_spec.max())
cmap = plt.get_cmap('jet')
# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels. (https://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)


lon_arr = np.zeros([nlon+1])
for ilon in range(nlon+1):
    lon_arr[ilon]  = minLonGrid + ilon * dlon
lat_arr = np.zeros([nlat+1])
for ilat in range(nlat+1):
    lat_arr[ilat]  = minLatGrid + ilat * dlat


im = m.pcolormesh(lon_arr, lat_arr, grid_spec, cmap=cmap, norm=norm)#linewidth = 0.01,  edgecolors = 'lightgrey')  # !!!!!! be careful: X and Y need to be the number of element of data + 1 (),  otherwise one cell
                                                                                              # will be missing (https://matplotlib.org/basemap/api/basemap_api.html#mpl_toolkits.basemap.Basemap.pcolormesh)
im.set_rasterized(True) # to hide the edges of the bins in the grid
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax)

# cax = fig.add_subplot(gs[0, 0])
# #cax.axis('off')
# fig.colorbar(im, ax=cax,  orientation='horizontal')
# gs.update(left = 0.11, right=0.87, top = 1,bottom = 0.12, hspace = -0.3)

#ax.margins(0,0)
fig_save_name = 'grid_sp_' + input_filename.replace('.txt', '.pdf')
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')






