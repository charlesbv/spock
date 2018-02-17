import numpy as np
from matplotlib import pyplot as plt
from read_input_file import *
from read_output_file import *
import sys
from get_prop_dir import *
from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv

plt.ion()

input_file_prop_name = get_prop_dir(1) + "run_cadre/input/main_input/" + sys.argv[1]
input_file_tle_name = get_prop_dir(1) + "run_cadre/input/tle/" + sys.argv[2]

path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(1) + 'output/python_propagator/'
name_mission = 'CADRE'
root_save_fig_name = path_folder_results + name_mission + '/result/image/' + sys.argv[1].replace(".txt","_") 

# Read tle input file
## Open file
input_file_tle = open(input_file_tle_name, "r")
read_input_file_tle = input_file_tle.readlines()
## Skip header
n_header_tle = 0
while ( read_input_file_tle[n_header_tle].split()[0] != '1' ):
    n_header_tle = n_header_tle + 1
n_tle = len(read_input_file_tle) - n_header_tle - 1 # - 1 to remove last line of file called '<End of file>'
## Read the elements of each tle of the file
nb_tle = n_tle / 2
epoch_tle = [] # Note: the cript propagate of sgp4 only propagate with an accuracy of one second (so the epochs can vary by one second)
r_eci_tle = np.zeros([nb_tle, 3])
v_eci_tle = np.zeros([nb_tle, 3])
for itle in range(nb_tle):
    line1 = read_input_file_tle[itle*2 + n_header_tle].replace("\r", "").replace("\n","") 
    line2 = read_input_file_tle[itle*2 + 1 + n_header_tle].replace("\r", "").replace("\n","")
    sc_tle = twoline2rv(line1, line2, wgs84)
    r_eci_tle[itle, :], v_eci_tle[itle, :] = sc_tle.propagate( sc_tle.epoch.year, sc_tle.epoch.month, sc_tle.epoch.day, sc_tle.epoch.hour, sc_tle.epoch.minute, sc_tle.epoch.second ) 
    epoch_tle.append( sc_tle.epoch )

# Read propagator input file
input_var, order_input_var = read_input_file(input_file_prop_name)
dt = input_var[2]; nb_steps = input_var[3]; satellite_to_plot_path = input_var[6]; satellite_to_plot = input_var[7]; 

# Read propagator output file
isat_prop = 0
r_eci_prop = np.zeros([nb_steps])
var_to_read = ["position"]
var_out, var_out_order = read_output_file( satellite_to_plot_path[isat_prop] + satellite_to_plot[isat_prop], var_to_read )
r_eci_prop = var_out[1]
date_prop = var_out[0]

# Linear interpolate the prop positions at the epochs of the tle
## Compute the number of seconds between each time step of the propagation and a reference time
refce = epoch_tle[0]
seconds_refce_until_date_prop = np.zeros([nb_steps])
for istep in range(nb_steps):
    seconds_refce_until_date_prop[istep] = ( datetime.strptime( date_prop[istep], "%Y/%m/%d %H:%M:%S" ) - refce ).days*24*3600 + ( datetime.strptime( date_prop[istep], "%Y/%m/%d %H:%M:%S" ) - refce ).seconds

## Interpolate at these times
r_eci_prop_interpolated_at_tle_epoch = []
r_eci_tle_that_are_during_prop = []
nb_seconds_since_refce_at_tle_that_are_during_prop = []
for itle in range(nb_tle):
    nb_seconds_since_refce = ( epoch_tle[itle] - refce ).days*24*3600 + ( epoch_tle[itle] - refce ).seconds
    if ( ( nb_seconds_since_refce >= seconds_refce_until_date_prop[0] ) & ( nb_seconds_since_refce <= seconds_refce_until_date_prop[-1] ) ): # only look at tle that have epoch more recent than the epoch start of the propagation and older than the epoch end of the propagation
        r_eci_tle_that_are_during_prop.append( r_eci_tle[itle, :])
        r_eci_prop_interpolated_at_tle_epoch_sublist = []
        nb_seconds_since_refce_at_tle_that_are_during_prop.append( nb_seconds_since_refce )
        where_interval_time = [ np.where( seconds_refce_until_date_prop >= nb_seconds_since_refce )[0][0]-1, np.where( seconds_refce_until_date_prop >= nb_seconds_since_refce )[0][0] ]
        a = ( r_eci_prop[where_interval_time[1], 0] - r_eci_prop[where_interval_time[0], 0] ) / ( seconds_refce_until_date_prop[where_interval_time[1]] - seconds_refce_until_date_prop[where_interval_time[0]])
        b = r_eci_prop[where_interval_time[0], 0] - a * seconds_refce_until_date_prop[where_interval_time[0]] 
        r_eci_prop_interpolated_at_tle_epoch_sublist.append( a*nb_seconds_since_refce + b )
        a = ( r_eci_prop[where_interval_time[1], 1] - r_eci_prop[where_interval_time[0], 1] ) / ( seconds_refce_until_date_prop[where_interval_time[1]] - seconds_refce_until_date_prop[where_interval_time[0]])
        b = r_eci_prop[where_interval_time[0], 1] - a * seconds_refce_until_date_prop[where_interval_time[0]] 
        r_eci_prop_interpolated_at_tle_epoch_sublist.append( a*nb_seconds_since_refce + b )
        a = ( r_eci_prop[where_interval_time[1], 2] - r_eci_prop[where_interval_time[0], 2] ) / ( seconds_refce_until_date_prop[where_interval_time[1]] - seconds_refce_until_date_prop[where_interval_time[0]])
        b = r_eci_prop[where_interval_time[0], 2] - a * seconds_refce_until_date_prop[where_interval_time[0]] 
        r_eci_prop_interpolated_at_tle_epoch_sublist.append( a*nb_seconds_since_refce + b )
        r_eci_prop_interpolated_at_tle_epoch.append(r_eci_prop_interpolated_at_tle_epoch_sublist)

r_eci_prop_interpolated_at_tle_epoch = np.array(r_eci_prop_interpolated_at_tle_epoch)
r_eci_tle_that_are_during_prop = np.array(r_eci_tle_that_are_during_prop)
nb_tle_that_are_during_prop = len(nb_seconds_since_refce_at_tle_that_are_during_prop)
nb_seconds_since_refce_at_tle_that_are_during_prop = np.array(nb_seconds_since_refce_at_tle_that_are_during_prop)

# ## Plot to check the linear interpolation
# fig, ax = plt.subplots(figsize=(20, 14))
# ax.plot( seconds_refce_until_date_prop, r_eci_prop[:, 2], 'b' )
# ax.scatter( nb_seconds_since_refce_at_tle_that_are_during_prop, z_eci_prop_interpolated_at_tle_epoch, color = 'r' )
# ax.set_xlim([nb_seconds_since_refce_at_tle_that_are_during_prop[10], nb_seconds_since_refce_at_tle_that_are_during_prop[-1]])

# Distance between the propagation and the tle positions
dist = np.zeros([nb_tle_that_are_during_prop])
for itle in range(nb_tle_that_are_during_prop):
    dist[itle] = np.linalg.norm( r_eci_tle_that_are_during_prop[itle, :] - r_eci_prop_interpolated_at_tle_epoch[itle, :])

# Plot the distance
fontsize = 16
fig, ax = plt.subplots(figsize = (20, 13), facecolor = 'w', edgecolor = 'k')
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                        
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# time ellapsed since the first tle that is in the interval of the propagation (in days):
x_axis = ( nb_seconds_since_refce_at_tle_that_are_during_prop - nb_seconds_since_refce_at_tle_that_are_during_prop[0] )/3600. / 24
ax.set_title('MinXSS ' + ' '.join(input_file_prop_name.split('/')[-1].split('.')[0].split('_')[1:]).title(), fontsize = fontsize, weight = 'bold')
ax.set_xlabel('Time since first TLE (days)', fontsize = fontsize, weight = 'bold')
ax.set_ylabel('Distance to TLE (km)', fontsize = fontsize, weight = 'bold')
ax.tick_params(axis = 'both', which = 'major', width = 2, pad = 6,labelsize = fontsize, size = 10)

ax.plot(x_axis, dist, linewidth = 2, color = 'k')
ax.margins(0,0)

name_fig_temp = 'distance_propagation_to_tle'
name_fig = root_save_fig_name + name_fig_temp + '.png'
fig.savefig( name_fig, facecolor = fig.get_facecolor(), edgecolor = 'none', bbox_inches = 'tight' )
os.system("rsync -av " + name_fig + " srbwks2014-0008.engin.umich.edu:./CADRE/")
