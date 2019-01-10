# Algorith
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt

from datetime import datetime, timedelta
import numpy as np
from read_input_file import *
from os import listdir
from os.path import isfile, join
deg_sign = u'\xb0'.encode('utf8')


main_input_file_name = sys.argv[1]

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt_output = input_variables[find_in_read_input_order_variables(order_input_variables, 'dt_output')]
dt = input_variables[find_in_read_input_order_variables(order_input_variables, 'dt')]
nb_steps = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_steps')]
satellite_to_plot_path_arr = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_path_list')]
satellite_to_plot_arr = input_variables[find_in_read_input_order_variables(order_input_variables, 'output_file_name_list')]
filename_station = input_variables[find_in_read_input_order_variables(order_input_variables, 'filename_ground_station')]
nb_satellites = input_variables[find_in_read_input_order_variables(order_input_variables, 'nb_sc')]
date_start = input_variables[find_in_read_input_order_variables(order_input_variables, 'date_start')]
date_stop = input_variables[find_in_read_input_order_variables(order_input_variables, 'date_stop')]
# read coverage ground station files
overflight_all_sat_all_station = []


for isat in range(nb_satellites):
    overflight_per_sat_all_station = []
    satellite_to_plot_path = satellite_to_plot_path_arr[isat]
    satellite_to_plot = satellite_to_plot_arr[isat]
    all_station_files = [satellite_to_plot_path + 'coverage/' + f for f in listdir(satellite_to_plot_path + 'coverage/') if ( isfile(join(satellite_to_plot_path + 'coverage/', f)) and f.endswith('.txt') and f.startswith("report") == 0 )]
    if isat == 0:
        nb_stations = len(all_station_files)
        nb_seconds_from_start_to_end = (int)((nb_steps-1)*dt_output)
        time_coverage = np.zeros([nb_seconds_from_start_to_end, nb_satellites, nb_stations])
    station_lon = np.zeros([nb_stations]); station_lat = np.zeros([nb_stations]); station_alt = np.zeros([nb_stations])
    overflight = []
    overflight_all_stations_combined = []
    for istation in range( nb_stations ):
        overflight_per_sat = []
        file_station = open( all_station_files[istation] )
        read_file_station = file_station.readlines()
        print all_station_files[istation]
        # skip header
        skip_header = 0
        while ( read_file_station[skip_header].split()[0] != '#START'):
            if 'longitude' in read_file_station[skip_header]:
                station_lon[istation] = np.float( read_file_station[skip_header].split(':')[-1].split()[0] )
                station_lat[istation] = np.float( read_file_station[skip_header].split(':')[-1].split()[1] )
                station_alt[istation] = np.float( read_file_station[skip_header].split(':')[-1].split()[2] )
            skip_header = skip_header + 1
        skip_header = skip_header + 1
        # read over the file 
        istep = 0
        new_coverage = 0
        while  istep + skip_header < len(read_file_station):
            if  istep + skip_header < len(read_file_station):
                new_coverage = 1
                overflight_per_sat_per_overflight = []
                if read_file_station[istep + skip_header + 1].split()[0].split(':')[2] == '60':
                    read_file_station[istep + skip_header + 1] = read_file_station[istep + skip_header + 1].replace(':60', ':59')
                if read_file_station[istep + skip_header].split()[0].split(':')[2] == '60':
                    read_file_station[istep + skip_header] = read_file_station[istep + skip_header].replace(':60', ':59')

                overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header].split()[0] )
                nb_seconds_from_date_start_until_start_coverage = (int)(( datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S" ) - date_start ).total_seconds())
#                while ( ( istep + skip_header < len(read_file_station) ) & ( np.float( read_file_station[istep + skip_header].split()[4] ) > -998 ) ):
                    
                while ( ( istep + 1 + skip_header < len(read_file_station) ) & ( (datetime.strptime(read_file_station[istep + skip_header + 1].split()[0], "%Y-%m-%dT%H:%M:%S") - datetime.strptime(read_file_station[istep + skip_header].split()[0], "%Y-%m-%dT%H:%M:%S")).total_seconds() < dt*1.001   ) ):
                    #print read_file_station[istep + skip_header ]
                    istep = istep + 1

                    if istep + skip_header + 1 == len(read_file_station):
                        break
                    else:
                        if read_file_station[istep + skip_header + 1].split()[0].split(':')[2] == '60':
                            read_file_station[istep + skip_header + 1] = read_file_station[istep + skip_header + 1].replace(':60', ':59')
                        if read_file_station[istep + skip_header].split()[0].split(':')[2] == '60':
                            read_file_station[istep + skip_header] = read_file_station[istep + skip_header].replace(':60', ':59')

                istep = istep + 1
                overflight_per_sat_per_overflight.append( read_file_station[istep + skip_header - 1].split()[0] )
                nb_seconds_from_date_start_until_end_coverage = (int)( (datetime.strptime(read_file_station[istep + skip_header - 1].split()[0], "%Y-%m-%dT%H:%M:%S" ) - date_start ).total_seconds())
                time_coverage[nb_seconds_from_date_start_until_start_coverage:nb_seconds_from_date_start_until_end_coverage+1, isat, istation] = 1
            if new_coverage == 1:
                overflight_per_sat.append( overflight_per_sat_per_overflight )
                overflight_per_sat_per_overflight.append(all_station_files[istation].split('/')[-1].split('_by')[0])
                overflight_all_stations_combined.append( overflight_per_sat_per_overflight )
                new_coverage = 0
        overflight.append( overflight_per_sat )
    overflight_all_sat_all_station.append(overflight_per_sat_all_station)

raise Exception
nb_station_in_sight_per_sat = np.zeros([nb_satellites, 2,nb_seconds_from_start_to_end]) # 2 for south vs north
for isat in range(nb_satellites):
    for itime in range(nb_seconds_from_start_to_end):
        nb_station_in_sight_per_sat[isat,0, itime] = (int)(len(np.where(time_coverage[itime, isat, :nb_stations/2] == 1)[0])) # north # !!!!! assumes that same number of stations north and south AND all_station_files has first all north then all south
        nb_station_in_sight_per_sat[isat,1, itime] = (int)(len(np.where(time_coverage[itime, isat, nb_stations/2:] == 1)[0])) # south
    
# Plot
# Plot
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

fig_title = 'Number of stations in sight as a function of time'
y_label = '# stations in sight'
x_label = 'Time (days)'
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

color_arr = ['b', 'r','cornflowerblue','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']
x_axis = np.arange(0,nb_seconds_from_start_to_end / 3600. / 24,1/3600. / 24)
#for isat in range(nb_satellites):
isat = 0
ax.scatter(x_axis[::300], nb_station_in_sight_per_sat[isat, 0, ::300], linewidth = 2, color = 'b', label = 'North')
ax.scatter(x_axis[::300], nb_station_in_sight_per_sat[isat, 1, ::300], linewidth = 2, color = 'r', label = 'South')
#ax.scatter(x_axis, nb_station_in_sight_per_sat[isat, 0, :] + nb_station_in_sight_per_sat[isat, 1, :], linewidth = 2, color = 'k', label = 'Both')

legend = ax.legend(loc='upper right', numpoints = 1,  title="", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = main_input_file_name.replace(".txt", "") + '_nb_stations_in_sight.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


print '% of time with 12 North stations on sight:', len(np.where(nb_station_in_sight_per_sat[isat, 0, :] == 12)[0]) * 100./ len(nb_station_in_sight_per_sat[isat, 0, :])
print '% of time with 12 South stations on sight:', len(np.where(nb_station_in_sight_per_sat[isat, 1, :] == 12)[0]) * 100./ len(nb_station_in_sight_per_sat[isat, 1, :])
