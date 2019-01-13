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

# copy of  collision_2dmap.py on aug 5 2018. modified for cygnss collisioona voidance paper (CA3).
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FixedLocator
from mpl_toolkits.basemap import Basemap, shiftgrid
from read_output_file import *
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from seconds_since_epoch_start_to_utc import *

from read_input_file import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta
from read_collision_file_new_with_dca import *

show_plot = 0
save_results = 1
root_save_fig_name = './'
name_mission = 'other'
if show_plot == 1:
    plt.ion()





run_list = ['oct31_h06_f107_ap_actual_2dgeo.txt']

nb_runs = len(run_list)
# GENERATE DISCTINCT COLORS
NCURVES = 100#nb_runs
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
# gs_with_map = gridspec.GridSpec(1, 1)
# gs_with_map.update(left=0.15, right=0.82, top = 0.90,bottom = 0.1)
# ax1 = fig.add_subplot(gs_with_map[0, 0])
min_lat = 400
min_lon = 400
max_lat = -400
max_lon = -400

# Read collision report to get TCA
for irun in range(nb_runs):
    input_filename = run_list[irun]
    param_in, param_in_name = read_input_file(input_filename)
    output_file_path_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_path_list')]
    output_file_name_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_name_list')]
    dt = param_in[find_in_read_input_order_variables(param_in_name, 'dt_output')]
    nb_satellites = param_in[find_in_read_input_order_variables(param_in_name, 'nb_sc')]
    nb_steps = param_in[find_in_read_input_order_variables(param_in_name, 'nb_steps')]
    date_start = param_in[find_in_read_input_order_variables(param_in_name, 'date_start')]
    output_run_dir_name = output_file_path_list[0].split('/')[-3]
    collision_filename = output_run_dir_name + '/' + output_run_dir_name + '_collision.txt'
    date_collision, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll = read_collision_file_new_with_dca( collision_filename )
    dt_coll = ( date_collision[0][1] - date_collision[0][0] ).total_seconds()
    nb_ca = cpc.shape[0]
    if nb_ca > 1:
        print "For now script works only if there is one single close approach. Here there are " + str(nb_ca) + " close approaches. The program will stop."; raise Exception

    # Figure out which time step in propagation the TCA corresponds to
    tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
    index_tca_in_propagation = (int) ( ( tca - date_start ).total_seconds() / dt )
    #!!!!!!!!!!!!!! comment block below if you don't understand it. this is here because for these 2 runs the time step is one below the time step for all other runs
    if ( ( irun == 7 ) | ( irun == 8 ) ):
        index_tca_in_propagation = index_tca_in_propagation + 1
    #!!!!!!!!!!!!!! end of comment block below if you don't understand it
    # Read position of 2 sc 
    list_var = ["position", "latitude", "longitude"]
    eci_r = np.zeros([nb_satellites, nb_steps, 3])
    lat = np.zeros([nb_satellites, nb_steps])
    lon = np.zeros([nb_satellites, nb_steps])
    for isat in range(nb_satellites):
        out_var, order_var = read_output_file( output_file_path_list[isat] + output_file_name_list[isat], list_var  )
        date = out_var[find_in_read_input_order_variables(order_var, 'date')]
        latitude = out_var[find_in_read_input_order_variables(order_var, 'latitude')]
        longitude = out_var[find_in_read_input_order_variables(order_var, 'longitude')] 
        position = out_var[find_in_read_input_order_variables(order_var, 'position')]
        eci_r[isat, :, :] = position
        lat[isat, :] = latitude
        lon[isat, :] = longitude
        # print lon[isat, index_tca_in_propagation-1] , lat[isat, index_tca_in_propagation-1]
        print lon[isat, index_tca_in_propagation] , lat[isat, index_tca_in_propagation]
        # print lon[isat, index_tca_in_propagation+1] , lat[isat, index_tca_in_propagation+1]
        if lon[isat, index_tca_in_propagation] < min_lon:
            min_lon = lon[isat, index_tca_in_propagation] 
        if lat[isat, index_tca_in_propagation] < min_lat:
            min_lat = lat[isat, index_tca_in_propagation] 
        if lon[isat, index_tca_in_propagation] > max_lon:
            max_lon = lon[isat, index_tca_in_propagation] 
        if lat[isat, index_tca_in_propagation] > max_lat:
            max_lat = lat[isat, index_tca_in_propagation] 

range_lon = max_lon - min_lon
range_lat = max_lat - min_lat
min_lon_array = [-180, min_lon - range_lon / 10, min_lon - 30]# 152.08] #0
min_lat_array = [-90, min_lat -  range_lat / 10, min_lat - 15]#21.17] #-10
max_lon_array = [180, max_lon + range_lon / 10,  min_lon + 30]#153.35] #60
max_lat_array = [90, max_lat + range_lat / 10, min_lat + 15]#22.32]#40
step_lon = [60,range_lon / 10,5]#10
step_lat = [30, range_lat / 10,5]#10
# array_lon = [['180W', '120W', '60W', '0', '60E', '120E', '180E'],['90W', '80W', '70W', '60W', '50W']]
# array_lat = [['90S', '60S', '30S', 'EQ', '30N', '60N', '90N'],['10N', '20N', '30N','40N','50N']]
zoom = 2
array_lon = [ str(ss) for ss in np.arange(min_lon_array[zoom], max_lon_array[zoom],step_lon[zoom]) ]
array_lat = [ str(ss) for ss in np.arange(min_lat_array[zoom], max_lat_array[zoom],step_lat[zoom]) ]

# Set up 2d map

width_fig = 13.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9

fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
        #ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
        #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure


# ax1.xaxis.set_major_locator(FixedLocator(np.arange(min_lon_array[zoom], max_lon_array[zoom]+1, step_lon[zoom])))
# ax1.yaxis.set_major_locator(FixedLocator(np.arange(min_lat_array[zoom], max_lat_array[zoom]+1, step_lat[zoom])))
# ax1.set_xticklabels(array_lon, color = 'k', fontsize = fontsize_plot, weight = 'bold')
# ax1.set_yticklabels(array_lat, color = 'k', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Longitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Latitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')

m = Basemap( projection       = 'cyl',
	     llcrnrlon        = min_lon_array[zoom] , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = max_lon_array[zoom]  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = min_lat_array[zoom]  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = max_lat_array[zoom],   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax1 ,
	     )

m.drawcoastlines(linewidth=1.7, color='black')
#m.bluemarble()                                                                                                               

# ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
marker_array = ['o', 's']
color_array = ['b', 'r']

irun = 0

input_filename = run_list[irun]
param_in, param_in_name = read_input_file(input_filename)
output_file_path_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_path_list')]
output_file_name_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_name_list')]
dt = param_in[find_in_read_input_order_variables(param_in_name, 'dt_output')]
nb_satellites = param_in[find_in_read_input_order_variables(param_in_name, 'nb_sc')]
nb_steps = param_in[find_in_read_input_order_variables(param_in_name, 'nb_steps')]
date_start = param_in[find_in_read_input_order_variables(param_in_name, 'date_start')]
output_run_dir_name = output_file_path_list[0].split('/')[-3]
collision_filename = output_run_dir_name + '/' + output_run_dir_name + '_collision.txt'
date_collision, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll = read_collision_file_new_with_dca( collision_filename )
dt_coll = ( date_collision[0][1] - date_collision[0][0] ).total_seconds()
nb_ca = cpc.shape[0]
if nb_ca > 1:
    print "For now script works only if there is one single close approach. Here there are " + str(nb_ca) + " close approaches. The program will stop."; raise Exception

# Figure out which time step in propagation the TCA corresponds to
tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
index_tca_in_propagation = (int) ( ( tca - date_start ).total_seconds() / dt )
#!!!!!!!!!!!!!! comment block below if you don't understand it. this is here because for these 2 runs the time step is one below the time step for all other runs
if ( ( irun == 7 ) | ( irun == 8 ) ):
    index_tca_in_propagation = index_tca_in_propagation + 1
#!!!!!!!!!!!!!! end of comment block below if you don't understand it

# Read position of 2 sc 
list_var = ["position", "latitude", "longitude"]
eci_r = np.zeros([nb_satellites, nb_steps, 3])
lat = np.zeros([nb_satellites, nb_steps])
lon = np.zeros([nb_satellites, nb_steps])
for isat in range(nb_satellites):
    out_var, order_var = read_output_file( output_file_path_list[isat] + output_file_name_list[isat], list_var  )
    date = out_var[find_in_read_input_order_variables(order_var, 'date')]
    latitude = out_var[find_in_read_input_order_variables(order_var, 'latitude')]
    longitude = out_var[find_in_read_input_order_variables(order_var, 'longitude')] 
    position = out_var[find_in_read_input_order_variables(order_var, 'position')]
    eci_r[isat, :, :] = position
    lat[isat, :] = latitude
    lon[isat, :] = longitude
    for istep in range(len(lon[isat, :])):
        if lon[isat, istep] > 180:
            lon[isat, istep] = lon[isat, istep] - 360

    colorVal = scalarMap.to_rgba(irun)
    # if irun == index_median_run:
    #     # x, y =  m(lon[isat, index_tca_in_propagation-1:index_tca_in_propagation+2] , lat[isat, index_tca_in_propagation-1:index_tca_in_propagation+2])
    #     # m.plot(x, y,  marker='s', markersize=10,color = colorVal,markeredgecolor = 'none',linestyle="None", label = run_list[irun].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%"  )[0]
    # else:
    #        x, y =  m(lon[isat, index_tca_in_propagation: index_tca_in_propagation+2] , lat[isat,index_tca_in_propagation: index_tca_in_propagation+2])
    nb_steps_in_one_orbit = (int)(95 * 60. / dt)
    x, y =  m(lon[isat, index_tca_in_propagation - nb_steps_in_one_orbit / 2:index_tca_in_propagation + nb_steps_in_one_orbit / 2 ] , lat[isat,index_tca_in_propagation  - nb_steps_in_one_orbit / 2: index_tca_in_propagation + nb_steps_in_one_orbit / 2 ])
    print isat, lon[isat, index_tca_in_propagation], lat[isat, index_tca_in_propagation]
    if isat == 0:
        m.plot(x, y,   markersize=10,color = 'blue',markeredgecolor = 'none', label = run_list[irun].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%", linewidth = 3)[0]
    else:
        m.plot(x, y,   markersize=10,color = 'red',markeredgecolor = 'none', linewidth = 3)[0]
# legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quartile of\nF10.7/Ap dist.", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))



if save_results == 1:
    fig_save_name = 'visualization_collision_python'
    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )


        # x, y =  m(lon[isat, index_tca_in_propagation] , lat[isat, index_tca_in_propagation])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'r')[0]
        # x, y =  m(lon[isat, index_tca_in_propagation+1] , lat[isat, index_tca_in_propagation+1])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'b')[0]

raise Exception


############ ZOOM ON COLLISION POINT
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
        #ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
        #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure


# ax1.xaxis.set_major_locator(FixedLocator(np.arange(min_lon_array[zoom], max_lon_array[zoom]+1, step_lon[zoom])))
# ax1.yaxis.set_major_locator(FixedLocator(np.arange(min_lat_array[zoom], max_lat_array[zoom]+1, step_lat[zoom])))
# ax1.set_xticklabels(array_lon, color = 'k', fontsize = fontsize_plot, weight = 'bold')
# ax1.set_yticklabels(array_lat, color = 'k', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Longitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Latitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')

# add_lat = 2.175e1
# add_lon = 1.56659e2 
# min_lon_collision = 0.003 + add_lon
# max_lon_collision = 0.006 + add_lon
# min_lat_collision = 0.002 + add_lat
# max_lat_collision = 0.0035 + add_lat

add_lat = 2.1751e1
add_lon = 1.56662e2
min_lon_collision = 0.0007 + add_lon 
max_lon_collision = 0.0026 + add_lon
min_lat_collision = 0.0020 + add_lat
max_lat_collision = 0.0024 + add_lat

m = Basemap( projection       = 'cyl',
	     llcrnrlon        = min_lon_collision , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = max_lon_collision  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = min_lat_collision  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = max_lat_collision,   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax1 ,
	     )

# ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
marker_array = ['o', 's']
color_array = ['b', 'r']

for irun in range(nb_runs):
    input_filename =  run_list[irun]
    param_in, param_in_name = read_input_file(input_filename)
    output_file_path_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_path_list')]
    output_file_name_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_name_list')]
    dt = param_in[find_in_read_input_order_variables(param_in_name, 'dt')]
    nb_satellites = param_in[find_in_read_input_order_variables(param_in_name, 'nb_satellites')]
    nb_steps = param_in[find_in_read_input_order_variables(param_in_name, 'nb_steps')]
    date_start = param_in[find_in_read_input_order_variables(param_in_name, 'date_start')]
    output_run_dir_name = output_file_path_list[0].split('/')[-3]
    collision_filename =  output_run_dir_name + '/' + output_run_dir_name + '_collision.txt'
    date_collision, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll = read_collision_file_new_with_dca( collision_filename )

    dt_coll = ( date_collision[0][1] - date_collision[0][0] ).total_seconds()
    nb_ca = cpc.shape[0]
    if nb_ca > 1:
        print "For now script works only if there is one single close approach. Here there are " + str(nb_ca) + " close approaches. The program will stop."; raise Exception

    # Figure out which time step in propagation the TCA corresponds to
    tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
    index_tca_in_propagation = (int) ( ( tca - date_start ).total_seconds() / dt )
    #!!!!!!!!!!!!!! comment block below if you don't understand it. this is here because for these 2 runs the time step is one below the time step for all other runs
    if ( ( irun == 7 ) | ( irun == 8 ) ):
        index_tca_in_propagation = index_tca_in_propagation + 1
    #!!!!!!!!!!!!!! end of comment block below if you don't understand it

    # Read position of 2 sc 
    list_var = ["position", "latitude", "longitude"]
    eci_r = np.zeros([nb_satellites, nb_steps, 3])
    lat = np.zeros([nb_satellites, nb_steps])
    lon = np.zeros([nb_satellites, nb_steps])
    for isat in range(nb_satellites):
        out_var, order_var = read_output_file( output_file_path_list[isat] + output_file_name_list[isat], list_var  )
        date = out_var[find_in_read_input_order_variables(order_var, 'date')]
        latitude = out_var[find_in_read_input_order_variables(order_var, 'latitude')]
        longitude = out_var[find_in_read_input_order_variables(order_var, 'longitude')] 
        position = out_var[find_in_read_input_order_variables(order_var, 'position')]
        eci_r[isat, :, :] = position
        lat[isat, :] = latitude
        lon[isat, :] = longitude

        colorVal = scalarMap.to_rgba(irun)
        # if irun == index_median_run:
        #     # x, y =  m(lon[isat, index_tca_in_propagation-1:index_tca_in_propagation+2] , lat[isat, index_tca_in_propagation-1:index_tca_in_propagation+2])
        #     # m.plot(x, y,  marker='s', markersize=10,color = colorVal,markeredgecolor = 'none',linestyle="None", label = run_list[irun].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%"  )[0]
        # else:
    #!!!!!!!!!!!!!! comment block below if you don't understand it. this is here because for these 2 runs the time step is one below the time step for all other runs
        if ( ( irun == 7 ) | ( irun == 8 ) ):
            x, y =  m(lon[isat, index_tca_in_propagation-1: index_tca_in_propagation+1] , lat[isat,index_tca_in_propagation-1: index_tca_in_propagation+1])
    #!!!!!!!!!!!!!! end of comment block below if you don't understand it
        else:
            x, y =  m(lon[isat, index_tca_in_propagation: index_tca_in_propagation+2] , lat[isat,index_tca_in_propagation: index_tca_in_propagation+2])
        print isat, lon[isat, index_tca_in_propagation], lat[isat, index_tca_in_propagation]
        if isat == 0:
            m.plot(x, y,  marker='s', markersize=10,color = colorVal,markeredgecolor = 'none', label = run_list[irun].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%", linewidth = 2  )[0]
        else:
            m.plot(x, y,  marker='o', markersize=10,color = colorVal,markeredgecolor = 'none' , linewidth = 2 )[0]
legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quartile of\nF10.7/Ap dist.", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

        # x, y =  m(lon[isat, index_tca_in_propagation] , lat[isat, index_tca_in_propagation])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'r')[0]
        # x, y =  m(lon[isat, index_tca_in_propagation+1] , lat[isat, index_tca_in_propagation+1])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'b')[0]


raise Exception

######### JUST MEDIAN COLLISION -  3 DT SHOWN
fig = plt.figure(num=None, figsize=(wxidth_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
        #ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
        #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure


# ax1.xaxis.set_major_locator(FixedLocator(np.arange(min_lon_array[zoom], max_lon_array[zoom]+1, step_lon[zoom])))
# ax1.yaxis.set_major_locator(FixedLocator(np.arange(min_lat_array[zoom], max_lat_array[zoom]+1, step_lat[zoom])))
# ax1.set_xticklabels(array_lon, color = 'k', fontsize = fontsize_plot, weight = 'bold')
# ax1.set_yticklabels(array_lat, color = 'k', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Longitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Latitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')


m = Basemap( projection       = 'cyl',
	     llcrnrlon        =  156., #Lower Left  CoRNeR Longitude
	     urcrnrlon        = 157.3  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = 21.1 , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = 22.35,   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax1 ,
	     )

# ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
marker_array = ['o', 's']
color_array = ['b', 'r']
irun = 4

input_filename =  run_list[irun]
param_in, param_in_name = read_input_file(input_filename)
output_file_path_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_path_list')]
output_file_name_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_name_list')]
dt = param_in[find_in_read_input_order_variables(param_in_name, 'dt')]
nb_satellites = param_in[find_in_read_input_order_variables(param_in_name, 'nb_satellites')]
nb_steps = param_in[find_in_read_input_order_variables(param_in_name, 'nb_steps')]
date_start = param_in[find_in_read_input_order_variables(param_in_name, 'date_start')]
output_run_dir_name = output_file_path_list[0].split('/')[-3]
collision_filename = output_run_dir_name + '/' + output_run_dir_name + '_collision.txt'
date_collision, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll = read_collision_file_new_with_dca( collision_filename )

dt_coll = ( date_collision[0][1] - date_collision[0][0] ).total_seconds()
nb_ca = cpc.shape[0]
if nb_ca > 1:
    print "For now script works only if there is one single close approach. Here there are " + str(nb_ca) + " close approaches. The program will stop."; raise Exception

# Figure out which time step in propagation the TCA corresponds to
tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
index_tca_in_propagation = (int) ( ( tca - date_start ).total_seconds() / dt )

# Read position of 2 sc 
list_var = ["position", "latitude", "longitude"]
eci_r = np.zeros([nb_satellites, nb_steps, 3])
lat = np.zeros([nb_satellites, nb_steps])
lon = np.zeros([nb_satellites, nb_steps])
for isat in range(nb_satellites):
    out_var, order_var = read_output_file( output_file_path_list[isat] + output_file_name_list[isat], list_var  )
    date = out_var[find_in_read_input_order_variables(order_var, 'date')]
    latitude = out_var[find_in_read_input_order_variables(order_var, 'latitude')]
    longitude = out_var[find_in_read_input_order_variables(order_var, 'longitude')] 
    position = out_var[find_in_read_input_order_variables(order_var, 'position')]
    eci_r[isat, :, :] = position
    lat[isat, :] = latitude
    lon[isat, :] = longitude

    colorVal = scalarMap.to_rgba(irun)

    x, y =  m(lon[isat, index_tca_in_propagation-1:index_tca_in_propagation+2] , lat[isat, index_tca_in_propagation-1:index_tca_in_propagation+2])
    if isat == 0:
        m.plot(x, y,  marker='s', markersize=10,color = 'k',markeredgecolor = 'none' )[0]
    if isat == 1:
        m.plot(x, y,  marker='o', markersize=10,color = 'k',markeredgecolor = 'none' )[0]

    print isat, lon[isat, index_tca_in_propagation], lat[isat, index_tca_in_propagation]

# legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quartile of\nF10.7/Ap dist.", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))



if save_results == 1:
    fig_save_name = 'position_only_median_sc_before_and_after_collision'
    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    #os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )


        # x, y =  m(lon[isat, index_tca_in_propagation] , lat[isat, index_tca_in_propagation])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'r')[0]
        # x, y =  m(lon[isat, index_tca_in_propagation+1] , lat[isat, index_tca_in_propagation+1])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'b')[0]

    
if show_plot == 1:
    plt.show(); plt.show()
