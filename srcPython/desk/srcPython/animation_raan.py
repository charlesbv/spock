from get_prop_dir import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.ticker import FixedLocator
import numpy as np
from read_input_file import *
import pickle 
from collections import *

plt.ion()
plt.isinteractive()

# Read the input file of propagator
input_filename = "/home/cbv/PropSim/input/main_input/aerie_other_angle.txt"
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; output_file_propagator = input_variables[6]

name_pickle = '/raid3/Armada/Charles/python/aerie_8sat_ok.pickle'
with open(name_pickle) as f:
    date, position, velocity, longitude, latitude, altitude, true_ano, raan, arg_perigee, right_asc, local_time, angle_asc_node_to_sat, spacing_s1_to_other_sat, spacing, spacing_no_low_plane,spacing_only_low_plane, spacing_only_other_inclination_plane, spacing_only_other_inclination_plane_with_pmpl, spacing_relative_s1_minus180_to_180 = pickle.load(f)

raise Exception
# ######################################################
# Set up the 2D map

#for month_start in np.arange(0,24,1):
month_start = 10
print month_start
fig = plt.figure()#plt.figure(figsize=(10, 8))
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
fig.suptitle('Month '+str(month_start), fontsize = 14, weight = 'bold')#fig.suptitle('Month '+str(month_start+1), fontsize = 14)
fig.set_facecolor('w')
m = Basemap( projection='npstere',boundinglat=0,lon_0=longitude[0,0]  )
#    m.bluemarble()                 

m.drawparallels(np.arange(-90.,91.,20.))
m.drawmeridians(np.arange(-180.,181.,20.),labels=[True,True,False,True], weight = 'bold')
m.drawcoastlines(linewidth=0.7, color='blue')
m.fillcontinents(color='coral',lake_color='aqua', alpha = 0.08)
# m.drawmapboundary(fill_color='aqua')

# # TUPLES
# Build the tuples for the visualization of the satellites
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')
spacecraft_list = []
name_satellite = ["" for x in range(nb_satellites)]

sat_for_spacing = [0,1,2,3,4,7]
sat_for_spacing_no_low_plane = [0,1,3,4]
sat_for_spacing_only_low_plane = [2, 7]
sat_for_spacing_only_other_inclination_plane = [5, 6]

sat_3_planes = [1, 7, 6] # 7 is 81/475, 6 is 83/500
colors = ['b','r','k']
colors_dot = ['b','b','r','b','b','k','k','r'] 
label_array = ['M1','M2','L1','M3','M4','I1','I2','L2']
marker_array = ['s','^','s','o','D','s','^','^']
# month_start = 12
# month_start= month_start-1
index_month_start = (int)(month_start * 30 *24*3600./dt)
nb_steps_in_one_orbit = (int)(95*60/dt)

i_count = -1
idot_count = -1
index_show_sat = 87
for i in sat_3_planes:
    i_count = i_count + 1
    spacecraft = namedtuple('spacecraft',('name',) +  point._fields +  color._fields + ('point_plot',) + ('marker_spacecraft',))
    spacecraft_list.append(spacecraft)
    spacecraft_list[i_count].x, spacecraft_list[i_count].y =  m(right_asc[i,index_month_start:index_month_start+nb_steps_in_one_orbit], latitude[i,index_month_start:index_month_start+nb_steps_in_one_orbit])
    spacecraft_list[i_count].point_plot = m.plot(spacecraft_list[i_count].x, spacecraft_list[i_count].y, markersize=15,color = colors[i_count], linewidth = 2)[0]
for idot in sorted(sat_for_spacing_no_low_plane + sat_for_spacing_only_low_plane + sat_for_spacing_only_other_inclination_plane):
    spacecraft_list.append(spacecraft)
    idot_count = idot_count + 1
    spacecraft_list[idot_count].x, spacecraft_list[idot_count].y =  m(right_asc[idot,index_month_start+index_show_sat], latitude[idot,index_month_start+index_show_sat])
    if latitude[idot,index_month_start+index_show_sat] > 0:
        spacecraft_list[idot_count].point_plot = m.scatter(spacecraft_list[idot_count].x, spacecraft_list[idot_count].y, color = colors_dot[idot], linewidth = 4, marker=marker_array[idot],label = label_array[idot],s = 200 , zorder= 5)
    
leg = plt.legend(ncol = idot_count+1, bbox_to_anchor=(0.5, 1.02), loc = 10,borderpad = 0.0001, frameon = False, fontsize = 14, scatterpoints=1,handlelength = 0.01)
for ileg in range(len(leg.legendHandles)):
    leg.legendHandles[ileg]._sizes = [40]

fig_save_name = get_prop_dir(2) + 'output/python_propagator/aerie/raan_month_' +str(month_start).zfill(2) + '_date_' + datetime.strftime( datetime.strptime(date[index_month_start+index_show_sat], "%Y/%m/%d %H:%M:%S"), "%y%m%d_%H%M%S") + '.png'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./")

    # for j in range(nb_steps_in_one_orbit):
    #     if ((latitude[i,index_month_start+j+1] > 0) & (latitude[i,index_month_start+j] < 0)):
    #         spacecraft_list[i_count].x, spacecraft_list[i_count].y =  m(longitude[i,index_month_start+j-10:index_month_start+j+nb_steps_in_one_orbit/2+10], latitude[i,index_month_start+j-10:index_month_start+j+nb_steps_in_one_orbit/2.+10])
    #         break

# DON'T ERASE THESE LINES COMMENTED BELOW    
# month_start =  1 / index_show_sat = 50 / date[index_month_start+index_show_sat] =  '2016/01/31 00:25:00'
# month_start =  6 / index_show_sat = 35 / date[index_month_start+index_show_sat] =  '2016/06/29 00:17:30'
# month_start =  10 / index_show_sat = 87 / date[index_month_start+index_show_sat] =  '2016/10/27 00:43:30'
# month_start =  16 / index_show_sat = 75 / date[index_month_start+index_show_sat] =  '2017/04/25 00:37:30'
# month_start =  22 / index_show_sat = 3 / date[index_month_start+index_show_sat] =  '2017/10/22 00:01:30'


latitude[sat_for_spacing_no_low_plane,index_month_start:index_month_start+nb_steps_in_one_orbit].transpose()
latitude[sat_for_spacing_only_other_inclination_plane,index_month_start:index_month_start+nb_steps_in_one_orbit].transpose()





