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
import matplotlib.gridspec as gridspec
from degree_to_time import *
from get_prop_dir import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.ticker import FixedLocator
import numpy as np
from read_input_file import *
import pickle 
from collections import *
import sys
from read_output_file import *
plt.ion()


input_filename = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2] 
# Read the input file of propagator
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_path_propagator = input_variables[6]; output_file_propagator = input_variables[7]

# Read the output
out_var_list_want = ["latitude", "right_asc", "local_time"]
isat = 0
out_var, out_var_list = read_output_file(output_path_propagator[isat] + output_file_propagator[isat], out_var_list_want)
latitude = out_var[1]; right_asc = out_var[2]; local_time = out_var[3]


# ######################################################
# Set up the 2D map

fig = plt.figure()#plt.figure(figsize=(10, 8))
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
fig.suptitle('Precession of the CYGNSS orbit over 50 days', fontsize = 14, weight = 'bold')
fig.set_facecolor('w')

m = Basemap( projection='npstere',boundinglat=0,lon_0=0 ,suppress_ticks=True )
#m.drawcoastlines(linewidth=0.7, color='blue')
#m.fillcontinents(color='coral',lake_color='aqua', alpha = 0.08)

m.drawparallels(np.arange(-90.,91.,20.))
meridian_ticks =  np.arange(0, 360,45.)
parallel_ticks = np.arange(-90.,91.,20.)

m.drawmeridians(meridian_ticks,labels=[0,0,0,0], weight = 'bold')
meridian_ticks_label_temp = degree_to_time(meridian_ticks, 1)
for ilabel in range(len(meridian_ticks_label_temp)):
    if ilabel%2 == 0:
        xtext, ytext = m(meridian_ticks[ilabel], 0)
#        plt.text( xtext, ytext, meridian_ticks_label_temp[ilabel], horizontalalignment = 'center' )



# # TUPLES
# Build the tuples for the visualization of the satellites
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')
spacecraft_list = []
name_satellite = ["" for x in range(nb_satellites)]
colors = ['b','r','k', 'g', 'm', 'y','c']
label_array = ['M1','M2','L1','M3','M4','I1','I2','L2']
nb_steps_in_one_orbit = (int)(95*60/dt) 


array_date = ['2004-03-21T12:00:00',
'2004-03-22T12:00:00',
'2004-03-23T12:00:00',
'2004-03-24T12:00:00',
'2004-03-25T12:00:00',
'2004-03-26T12:00:00',
'2004-03-27T12:00:00']
# '2004-03-28T12:00:00',
# '2004-03-29T12:00:00',
# '2004-03-30T12:00:00']



# array_date = ['2004-03-21T12:00:00',
#               '2004-03-31T12:00:00',
#               '2004-04-09T12:00:00',
#               '2004-04-19T12:00:00',
#               '2004-04-29T12:00:00']

nb_dates = len(array_date)
i_count = -1
for idate in range(nb_dates):
    i_count = i_count + 1
    nb_seconds_since_start = ( datetime.strptime(array_date[idate], "%Y-%m-%dT%H:%M:%S") - date_start ).total_seconds()
    nb_index_since_start = (int)( nb_seconds_since_start / dt )
    spacecraft = namedtuple('spacecraft',('name',) +  point._fields +  color._fields + ('point_plot',) + ('marker_spacecraft',))
    spacecraft_list.append(spacecraft)
    spacecraft_list[i_count].x, spacecraft_list[i_count].y = m(local_time[nb_index_since_start:nb_index_since_start+nb_steps_in_one_orbit], latitude[nb_index_since_start:nb_index_since_start+nb_steps_in_one_orbit])
    spacecraft_list[i_count].point_plot = m.plot(spacecraft_list[i_count].x, spacecraft_list[i_count].y, markersize=15, color=colors[i_count], linewidth = 2)[0] #colors[i_count]
    print colors[i_count], latitude[nb_index_since_start], latitude[nb_index_since_start+nb_steps_in_one_orbit], local_time[nb_index_since_start], local_time[nb_index_since_start+nb_steps_in_one_orbit]
    # spacecraft_list[i_count].x, spacecraft_list[i_count].y =  m(right_asc[nb_index_since_start:nb_index_since_start+nb_steps_in_one_orbit], latitude[nb_index_since_start:nb_index_since_start+nb_steps_in_one_orbit])
    # spacecraft_list[i_count].point_plot = m.plot(spacecraft_list[i_count].x, spacecraft_list[i_count].y, markersize=15,color = colors[i_count], linewidth = 2, linestyle = 'dotted')[0]
    

fig_save_name = './paper3_stk_results/prcessions_local_time_7days.eps'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none')#, bbox_inches='tight')  
os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")

    
# leg = plt.legend(ncol = idot_count+1, bbox_to_anchor=(0.5, 1.02), loc = 10,borderpad = 0.0001, frameon = False, fontsize = 14, scatterpoints=1,handlelength = 0.01)
# for ileg in range(len(leg.legendHandles)):
#     leg.legendHandles[ileg]._sizes = [40]




