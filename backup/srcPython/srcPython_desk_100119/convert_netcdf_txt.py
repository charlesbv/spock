from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, shiftgrid
import sys
import os
#sys.path.append('/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython')
sys.path.append('/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython')
from read_input_file import *
from cygnss_read_spock_spec import *
from netCDF4 import Dataset
from datetime import datetime, timedelta
from math import sin, cos, sqrt, atan2, radians


rad2deg = 180./np.pi
# Read positions of specular points from flight data (netCDF files have to be downloaded first (from cygnss-sftp-1 or cygnss-sftp-2))
filename_spec_flight = "/Users/cbv/cygnss/netcdf/2017/197/cyg01.ddmi.s20170716-000000-e20170716-235959.l1.power-brcs.sand037.nc" #"data/cyg01.ddmi.s20170820-000000-e20170820-235959.l1.power-brcs.sand014.nc" #cyg06.ddmi.s20170504-000039-e20170504-141945.l1.power-brcs.sand014.nc"
fh = Dataset(filename_spec_flight, mode='r')
# nc_attrs = fh.ncattrs()
# nc_dims = [dim for dim in fh.dimensions]  # list of nc dimensions
# nc_vars = [var for var in fh.variables]  # list of nc variables
            
sc_pos_x = fh.variables['sc_pos_x'][:]
sc_pos_y = fh.variables['sc_pos_y'][:]
sc_pos_z = fh.variables['sc_pos_z'][:]
sc_vel_x = fh.variables['sc_vel_x'][:]
sc_vel_y = fh.variables['sc_vel_y'][:]
sc_vel_z = fh.variables['sc_vel_z'][:]

sc_pitch = fh.variables['sc_pitch'][:] # Spacecraft pitch angle relative to the orbit frame, in radians at ddm_timestamp_utc
sc_roll = fh.variables['sc_roll'][:] 
sc_yaw = fh.variables['sc_yaw'][:] 
# sc_pitch_att # Spacecraft pitch angle relative to the orbit frame, in radians at att_timestamp_utc

# !!!!!!!!! should we use small_sc_attitude_err and large_sc_attitude_err below?
# small_sc_attitude_err = fh.variables['small_sc_attitude_err'][:] # Set if the absolute value of at least one of the spacecraft roll, pitch or yaw angles is between one and five degrees.
# large_sc_attitude_err = fh.variables['large_sc_attitude_err'][:] # Set if the absolute value of at least one of the spacecraft roll, pitch or yaw angles is greater than or equal to 5 degrees

time_flight = fh.variables['ddm_timestamp_utc'][:]
time_coverage_start = fh.getncattr(fh.ncattrs()[fh.ncattrs().index('time_coverage_start')])
time_coverage_start_datetime = datetime.strptime(time_coverage_start[:-4], "%Y-%m-%dT%H:%M:%S.%f") 
#fh.close()
nb_time_flight = len(sc_pos_x)
date_flight = []
for itime in range(nb_time_flight):
    date_flight.append( datetime.strftime(time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6)), "%Y-%m-%dT%H:%M:%S.%f" ) )

cyg_name = filename_spec_flight.split('/')[-1].split('.')[0]
# # Write time, position, velocity in txt file
# filename_out = "output/" + filename_spec_flight.split('/')[-1].split('.nc')[0] + ".txt"
# file_out = open(filename_out,"w")
# print >> file_out, "#ECEF position and velocity of " + cyg_name + " from " + date_flight[0] + " until " + date_flight[-1]
# print >> file_out, "#Date position(km/s) velocity(km/s) pitch(deg) roll(deg) yaw(deg)"
# print >> file_out, "#START"
# for itime in range(nb_time_flight):
#     if int(date_flight[itime][11:13]) >= 12:
#         print >> file_out, date_flight[itime], format(sc_pos_x[itime]/1000., "e"), format(sc_pos_y[itime]/1000., "e"), format(sc_pos_z[itime]/1000., "e"), format(sc_vel_x[itime]/1000., "e"), format(sc_vel_y[itime]/1000., "e"), format(sc_vel_z[itime]/1000., "e"), sc_pitch[itime]*rad2deg, sc_roll[itime]*rad2deg, sc_yaw[itime]*rad2deg
# file_out.close()

filename_att_out = "output/att_" + filename_spec_flight.split('/')[-1].split('.nc')[0] + ".txt"
file_att_out = open(filename_att_out,"w")
print >> file_att_out, "#BEGINNINGOFHEADER"
print >> file_att_out, "#" + cyg_name + " from " + date_flight[0] + " until " + date_flight[-1]
print >> file_att_out, "#ENDOFHEADER"
for itime in range(nb_time_flight):
    print >> file_att_out, date_flight[itime] +  " (" + str(sc_pitch[itime]*rad2deg) + "; " +  str(sc_roll[itime]*rad2deg) + "; " + str(sc_yaw[itime]*rad2deg) + ")" + " (1; 2; 3)"
print >> file_att_out, "#ENDOFFILE"
file_att_out.close()





# Plot pitch roll yaw
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

fig_title = ''
x_label = 'Time (hours since ' + time_coverage_start[:-4] +  ')' 
fig_att = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
x_axis = time_flight/3600.
fig_att.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(3, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.06)

ax_pitch = fig_att.add_subplot(gs[0, 0])
y_label = 'Pitch ' + u'(\N{DEGREE SIGN})'
ax_pitch.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#ax_pitch.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax_pitch.spines.itervalues()] # change the width of the frame of the figure
ax_pitch.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
y_axis = sc_pitch*rad2deg
ax_pitch.plot(x_axis, y_axis, linewidth = 2, color = 'b')
ax_pitch.xaxis.set_ticklabels("")
ax_pitch.margins(0,0)


ax_roll = fig_att.add_subplot(gs[1, 0])
y_label = 'Roll ' + u'(\N{DEGREE SIGN})'
ax_roll.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#ax_roll.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax_roll.spines.itervalues()] # change the width of the frame of the figure
ax_roll.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
y_axis = sc_roll*rad2deg
ax_roll.plot(x_axis, y_axis, linewidth = 2, color = 'b')
ax_roll.xaxis.set_ticklabels("")
ax_roll.margins(0,0)


ax_yaw = fig_att.add_subplot(gs[2, 0])
y_label = 'Yaw ' + u'(\N{DEGREE SIGN})'
ax_yaw.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_yaw.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax_yaw.spines.itervalues()] # change the width of the frame of the figure
ax_yaw.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
y_axis = sc_yaw*rad2deg
ax_yaw.plot(x_axis, y_axis, linewidth = 2, color = 'b')
ax_yaw.margins(0,0)


fig_save_name =  'attitude_' + filename_spec_flight.split('/')[-1].replace(".nc", ".pdf")
fig_att.savefig(fig_save_name, facecolor=fig_att.get_facecolor(), edgecolor='none', bbox_inches='tight')



