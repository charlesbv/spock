from degree_to_time import *
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from read_input_file import *
from read_output_file import *
from get_prop_dir import *
from compute_power import *
from datetime import datetime, timedelta
from orbit_average import *
from dateutil.relativedelta import relativedelta
plt.ion()

input_filename = 'aerie_power.txt'
normalization_coef = 10 # the run is with a area of 0.1 m^2, so multiply by 10 to get a normalized power (= power for an area of 1 m^2)
coeff_real = 0.4*0.3*4*0.28*0.74 # length + width * nb of panels * solar cell eff * packing eff


# READ PROPAGATOR INPUT FILE
input_filename = get_prop_dir(2) + 'input/main_input/' + input_filename
input_variables, order_input_variables = read_input_file(input_filename)
date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; nb_satellites = input_variables[4]; 
output_propagator_path = input_variables[6]; output_propagator_file = input_variables[7]; gps_name = input_variables[5]; 
nb_surfaces = input_variables[8]


# COMPUTE POWER
power = np.zeros([nb_satellites,nb_steps, nb_surfaces])
sun_elevation = np.zeros([nb_satellites,nb_steps])
for isat in range(nb_satellites):
    power[isat, :, :], sun_elevation[isat, :] = compute_power(output_propagator_path[isat] + 'power_'+ output_propagator_file[isat])
total_power = np.zeros([nb_satellites,nb_steps])
for isat in range(nb_satellites):
    for istep in range(nb_steps):
        total_power[isat, istep] = np.sum(power[isat, istep, :] ) *0.74 # !!!!!!!!!!!!!!! 0.74 for packing efficiency

# READ LATITUDE TO THEN COMPUTE ORBIT-AVERAGE POWER
latitude = np.zeros([nb_satellites, nb_steps])
local_time = np.zeros([nb_satellites, nb_steps])
list_var = ["latitude", "local_time"]
lat_var_nb = 1
local_time_var_nb = 2
for isat in range(nb_satellites):
    out_var, order_var = read_output_file(output_propagator_path[isat] + output_propagator_file[isat], list_var)
    if isat == 0:
        date = out_var[0]
    latitude[isat, :] = out_var[lat_var_nb]
    local_time[isat, :] = out_var[local_time_var_nb]
date = out_var[0]

# COMPUTE ORBIT-AVERAGE POWER 
total_power_average = []
for isat in range(nb_satellites):
    total_power_average_temp, date_averaged, index_date_averaged = orbit_average(total_power[isat, :], latitude[isat, :], date)
    total_power_average.append(total_power_average_temp)
index_date_averaged = np.array(index_date_averaged)


# WRITE IN AN OUTPUT FILE THE NORMALIZED POWER (POWER * normalization_coef)
iorbit = 0 # number of the orbit we want to plot the power of
isat = 0 # sat to plot the power of

duration_of_this_orbit_in_minutes = ( index_date_averaged[iorbit,2] - index_date_averaged[iorbit,0] ) * dt / 60.
x_axis = np.arange(0, duration_of_this_orbit_in_minutes, dt / 60.)

inclination = 81.5
filename_out = get_prop_dir(2) + 'output/python_propagator/aerie/aerie_polar_ltan_' + degree_to_time([local_time[isat, index_date_averaged[iorbit,0]]],0)[0].replace(":","_") + '.txt'
file_out = open(filename_out, "w")
print >> file_out, 'This file shows the normalized power (power for a solar panel area of 1 m^2) and the Sun elevation angle over an orbit for AERIE with a gimballed tail. Sun elevation is -999.0 if the spacecraft is in the shade. Time is time in orbit from the ascenging node (ex: 2 means the satellite went north 2 minutes after crossing the Equator).'
print >> file_out, 'The constant for the solar flux is 1358 W/m^2.\nInclination is ' +str(inclination) + ' degrees and local time is ' + degree_to_time([local_time[isat, index_date_averaged[iorbit,0]]],0)[0] +'.'
print >> file_out, 'TIME(min) NORMALIZED_POWER(W) SUN_ELEVATION(deg)'
for i in range(len(total_power[isat, index_date_averaged[iorbit,0]:index_date_averaged[iorbit,2] ])):
    print >> file_out, x_axis[i], total_power[isat, index_date_averaged[iorbit,0]+i ]*normalization_coef, sun_elevation[isat, index_date_averaged[iorbit,0] + i ]
file_out.close()
os.system("rsync -av " + filename_out +" srbwks2014-0008.engin.umich.edu:./")
# PLOT
## POWER OVER ONE ORBIT
height_fig = 11
fontsize_plot = 14
fontsize_title = 16
fig = plt.figure(num=None, figsize=(18, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(2, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = fontsize_title)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('AERIE normalized power * 0.4*0.3*4*0.28*0.74 (top) and Sun elevation (bottom) over one orbit - LTAN = ' + degree_to_time([local_time[isat, index_date_averaged[iorbit,0]]],0)[0], weight = 'bold', fontsize = fontsize_title,  y = 1.01) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, size = 10, labelsize = fontsize_plot)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure


ax1.plot( x_axis, total_power[isat, index_date_averaged[iorbit,0]:index_date_averaged[iorbit,2] ]*normalization_coef*coeff_real, label = 'Power', linewidth = 2, color = 'k')
ax1.plot( [min(x_axis), max(x_axis)], [total_power_average[isat][iorbit]*normalization_coef*coeff_real, total_power_average[isat][iorbit] * normalization_coef*coeff_real], label = 'Orbit-average power (' + '{0:.1f}'.format(total_power_average[isat][iorbit]*normalization_coef*coeff_real) + ' W)', linewidth = 2, color = 'k', linestyle = 'dashdot')
power_requirement = 26.44
ax1.plot( [min(x_axis), max(x_axis)], [power_requirement,power_requirement], label = 'Requirement', linewidth = 2, color = 'r')
ax1.plot( [min(x_axis), max(x_axis)], [power_requirement*1.3,power_requirement*1.3], label = '30 % margin', linewidth = 2, color = 'r', linestyle = 'dashdot')
ax1.set_ylabel('Power (W)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylim([0,ax1.get_ylim()[1]*1.1])
ax1.legend(loc = 3, ncol = 2)
ax1.margins(0,0)

ax2 = fig.add_subplot(gs[1, 0])
ax2.tick_params(axis='both', which='major',  width = 2, pad = 6, size = 10, labelsize = fontsize_plot)
[i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
duration_of_this_orbit_in_minutes = ( index_date_averaged[iorbit,2] - index_date_averaged[iorbit,0] ) * dt / 60.
ax2.plot( x_axis, sun_elevation[isat, index_date_averaged[iorbit,0]:index_date_averaged[iorbit,2] ], label = 'Elevation', linewidth = 2, color = 'k')
ax2.set_xlabel('Time in orbit (min)', fontsize = fontsize_plot, weight = 'bold')
ax2.set_ylabel('Sun elevation (' + u' \N{DEGREE SIGN}' ')', fontsize = fontsize_plot, weight = 'bold')
ax2.set_ylim([-90,90])
ax2.margins(0,1)
ax2.yaxis.set_ticks(np.arange(-90,91,30))
gs.update(left=0.05, right=0.99, top = 0.95,bottom = 0.05)
fig_save_name = filename_out.replace(".txt",".png")# get_prop_dir(2) + 'output/python_propagator/aerie/test_aerie_power_ltan_' + degree_to_time([local_time[isat, index_date_averaged[iorbit,0]]],0)[0].replace(":","_") + '.png'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./")
raise Exception
## ORBIT-AVERAGE POWER
height_fig = 11
fontsize_plot = 14
fig = plt.figure(num=None, figsize=(18, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.99, top = 0.7,bottom = 0.2)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('AERIE orbit-average power over a year', weight = 'bold', fontsize = 20,  y = 1.01) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, size = 10, labelsize = fontsize_plot)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

x_axis = index_date_averaged[:,1] * dt # time ellapsed in seconds for each middle of bin average

ax1.plot( x_axis, total_power_average[0], label = 'Total', linewidth = 2, color = 'k')
ax1.plot( [min(x_axis), max(x_axis)], [20,20], label = 'Requirement', linewidth = 2, color = 'r')
ax1.plot( [min(x_axis), max(x_axis)], [20*1.3,20*1.3], label = '30 % margin', linewidth = 2, color = 'r', linestyle = 'dashdot')
ax1.set_xlabel('Time (months)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Power (W)', fontsize = fontsize_plot, weight = 'bold')

ax1.legend(loc = 3)
date_list = [date_start + relativedelta(month=x) for x in np.arange(1, 13)]
# Aaron wants to go from 0 (1st of Jan) to 12 (31st of Dec)
x_ticks = np.zeros([len(date_list)]) 
x_ticklabels = []
x_ticklabels.append('0')
for imonth in range(1,len(date_list)):
    x_ticks_temp = (date_list[imonth] - date_list[0]).days*24*3600 + (date_list[imonth] - date_list[imonth-1]).seconds
    x_ticks[imonth] = x_ticks_temp
    x_ticklabels.append( str(imonth) )
ax1.xaxis.set_ticks(x_ticks)
ax1.xaxis.set_ticklabels(x_ticklabels, fontsize = fontsize_plot)#, rotation='vertical')
ax1.margins(0,1)
ax1.set_ylim([0, max(total_power_average[0])])


