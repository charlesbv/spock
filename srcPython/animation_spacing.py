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

from get_prop_dir import *
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.animation as animation
import pickle
from datetime import datetime, timedelta
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from read_input_file import *
from read_output_file import *
import colorsys
from matplotlib import gridspec

plt.ion()


# Read the input file of propagator
input_filename = "/home/cbv/PropSim/input/main_input/aerie_other_angle.txt"
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_file_propagator = input_variables[6]

name_pickle = '/raid3/Armada/Charles/python/aerie_8sat_ok.pickle'
with open(name_pickle) as f:
    date, position, velocity, longitude, latitude, altitude, true_ano, raan, arg_perigee, right_asc, local_time, angle_asc_node_to_sat, spacing_s1_to_other_sat, spacing, spacing_no_low_plane,spacing_only_low_plane, spacing_only_other_inclination_plane, spacing_only_other_inclination_plane_with_pmpl, spacing_relative_s1_minus180_to_180 = pickle.load(f)




sat_for_spacing = [0,1,2,3,4,7]
sat_for_spacing_no_low_plane = [0,1,3,4]
sat_for_spacing_only_low_plane = [2, 7]
sat_for_spacing_only_other_inclination_plane = [5, 6]
sat_for_spacing_without_reference_sat  = [1,2,3,4,7]

colors = ['b','b','r','b','b','k','k','r'] 
label_array = ['M1','M2','L1','M3','M4','I1','I2','L2']
marker_array = ['s','^','s','o','D','s','^','^']

height_fig = 5
fontsize_plot = 14
fig = plt.figure(num=None, figsize=(height_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

phi_array=np.arange(0,6.28,0.01)
r =1.
ax1.plot( r*np.cos(phi_array), r*np.sin(phi_array), 'k--', linewidth = 2 )
ax1.set_xlim([-r*1.08, r*1.08])
ax1.set_ylim([-r*1.08, r*1.08])
ax1.axis('off')

time_plot = '2017/10/22 00:01:30'   
month_start = 22
index_plot = date.index(time_plot)
fig.suptitle('Month ' + str(month_start), fontsize = 14, weight = 'bold')
for isat in sorted(sat_for_spacing):
    x1 = r*1.0*np.cos((angle_asc_node_to_sat[isat,index_plot] - angle_asc_node_to_sat[0,index_plot])*np.pi/180)                                  
    y1 = r*1.0*np.sin((angle_asc_node_to_sat[isat,index_plot] - angle_asc_node_to_sat[0,index_plot])*np.pi/180)                                                              
    ax1.scatter( [x1], [y1], marker = marker_array[isat], label = label_array[isat], s = 150, color= colors[isat],linewidth = 4, zorder = 5)
spacing_edge = 0.1
gs.update(left=spacing_edge/2, right=1-spacing_edge/2, top = 1-spacing_edge,bottom = 0)
leg = plt.legend(ncol = len(sat_for_spacing), bbox_to_anchor=(0.5, 1.02), loc = 10,borderpad = 0.0001, frameon = False, fontsize = 14, scatterpoints=1,handlelength = 0.01)
for ileg in range(len(leg.legendHandles)):
    leg.legendHandles[ileg]._sizes = [40]
fig.savefig(get_prop_dir(2) + 'output/python_propagator/aerie/spacing_month_' +str(month_start).zfill(2) + '_date_' + datetime.strftime( datetime.strptime(date[index_plot], "%Y/%m/%d %H:%M:%S"), "%y%m%d_%H%M%S") + '.png', facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
raise Exception
index_skip = 180
date_ani = date[0::index_skip]
x1 = r*1.0*np.cos((angle_asc_node_to_sat[:,0::index_skip] - angle_asc_node_to_sat[0,0::index_skip])*np.pi/180)                                  
y1 = r*1.0*np.sin((angle_asc_node_to_sat[:,0::index_skip] - angle_asc_node_to_sat[0,0::index_skip])*np.pi/180)                                                                                             
      

time_text = ax1.text(0.5, 0.95, '', horizontalalignment ='center', transform=ax1.transAxes, fontsize = 15)
line, = ax1.plot([], [], 'o', lw=2)
list_sat = []

for j in range(nb_satellites-2):  
    list_sat.append((ax1.plot([], [], 'p', lw=2, marker = marker_array[j], markersize = 15, color= colors[j])[0],))
    ax1.scatter(10,10, marker = marker_array[j], label = label_array[j], s = 70, color= colors[j])


# Init only required for blitting to give a clean slate.
ax1.legend(scatterpoints=1)
def init():
    tuple_spacecraft = ()
    for j in range(nb_satellites-2):
        list_sat[j][0].set_data([],[])
        tuple_spacecraft = tuple_spacecraft + list_sat[j]

    time_text.set_text('')
    tuple_spacecraft = tuple_spacecraft + (time_text,)
    return tuple_spacecraft

def animate(i):
    tuple_spacecraft = ()
    for j in range(nb_satellites-2):
        thisx1 = [10, x1[j,i]]
        thisy1 = [10, y1[j,i]]
        list_sat[j][0].set_data(thisx1, thisy1)
        tuple_spacecraft = tuple_spacecraft + list_sat[j]

    time_text.set_text(date_ani[i])
    tuple_spacecraft = tuple_spacecraft + (time_text,)
    fig.savefig('./images/animation/'+str(i), facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    return tuple_spacecraft

fps_we_want = 5L
ani = animation.FuncAnimation(fig, animate, init_func=init, frames = x1.shape[1],
                              interval=1000L/fps_we_want, blit=True, repeat = False)
plt.show()

# s1 is M1
# s2 is M2
# s3 is L1
# s4 is M3
# s5 is M4
# s6 is I1
# s7 is I2
# s8 is L2







    # for j in range(nb_satellites-2):
    #     thisx = [0, x1[j,i]]
    #     thisy = [0, y1[j,i]]
    #     line.set_data(thisx, thisy)
    #     list_sat.append((line,))
    # for j in range(len(list_sat)):
    #     tuple_spacecraft = tuple_spacecraft + list_sat[j]


 # 
#    print 'XXXXXXXXXXXXXXXXXXXXXXXXX', len(tuple_spacecraft)
  #  tuple_spacecraft = tuple_spacecraft + (time_text,)
 #   print 'XXXXXXXXXXXXXXXXXXXXXXXXX', len(tuple_spacecraft)
  #  raise Exception
#    fig.savefig('./images/'+str(i), facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')

########WORKS
##IN INI
    # line1.set_data([],[]) 
    # line2.set_data([],[]) 
    # tuple_spacecraft = tuple_spacecraft + (line1,) + (line2,)
##IN ANIME
    # list_sat = []
    # thisx1 = [0, x1[0,i]]
    # thisy1 = [0, y1[0,i]]
    # line1.set_data(thisx1, thisy1)
    # thisx2 = [0, x1[1,i]]
    # thisy2 = [0, y1[1,i]]
    # line2.set_data(thisx2, thisy2)
    # tuple_spacecraft = tuple_spacecraft + (line1,) + (line2,)
