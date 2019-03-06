# This script plots the
# - FM position in polar coord from the beacon station
# - GPS position in polar coord from the FM
# from a csv file output by tds-bop.exe and read by beacon_read_csv.py
import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import csv
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from ecef_to_lvlh import *
from beacon_read_csv import *
from struct import *
from tabulate import tabulate
deg2rad = np.pi/180
filename = 'outputCygnssOct/FM01/pass_3_PRN_21.csv'

date, prn, target_lat, target_lon, target_alt, target_ecef_x, target_ecef_y, target_ecef_z, target_rx_sat_look_angle_az, target_rx_sat_look_angle_el, target_rx_sat_range, sp_lat, sp_lon, sp_ecef_pos_x, sp_ecef_pos_y, sp_ecef_pos_z, sp_gain, rx_sub_sat_lat, rx_sub_sat_lon, rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z, rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z, tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z, rx_power = beacon_read_csv(filename)

rx_to_tx_x = tx_sat_ecef_pos_x - rx_sat_ecef_pos_x
n = len(rx_to_tx_x)
rx = np.array( [rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z]).transpose()
vx = np.array( [rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z]).transpose()
tx = np.array( [tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z]).transpose()
beacon = np.array( [target_ecef_x, target_ecef_y, target_ecef_z]).transpose()
rx_mag = np.linalg.norm(rx, axis = 1)

# FM TO GPS
rx_to_tx = tx - rx
rx_to_tx_lvlh = np.zeros([n, 3])
rx_to_tx_lvlh_norm =  np.zeros([n, 3])
rx_dot_rx_to_tx_lvlh = np.zeros([n])
z_orbit = np.array([0, 0, 1])
for i in range(n):
    rx_to_tx_lvlh[i, :] = ecef_to_lvlh(rx[i, :], vx[i, :], rx_to_tx[i, :])
    rx_to_tx_lvlh_norm[i, :] = rx_to_tx_lvlh[i, :] / np.linalg.norm(rx_to_tx_lvlh[i, :])
    rx_dot_rx_to_tx_lvlh[i] = np.dot(z_orbit, rx_to_tx_lvlh_norm[i, :])
#90. -  acos(cygnss_r_dot_C2G_orbit)/dtor
elev_tx_from_rx = 90. - np.arccos(rx_dot_rx_to_tx_lvlh) / deg2rad
#atan2(C2S_body[1], C2S_body[0])/dtor));
azim_tx_from_rx_pi_to_pi = np.arctan2(rx_to_tx_lvlh[:,1], rx_to_tx_lvlh[:,0]) / deg2rad
azim_tx_from_rx_90port = np.zeros([n])
for i in range(n):
    if (azim_tx_from_rx_pi_to_pi[i] < 0):
        azim_tx_from_rx_90port[i] = 360 + azim_tx_from_rx_pi_to_pi[i]
    else:
        azim_tx_from_rx_90port[i] = azim_tx_from_rx_pi_to_pi[i]
# Until now, azim_tx_from_rx_90port = 90 means tx is port of rx, 270 means tx is startboarf of rx, - means tx is on the ram side of rx, 180 means tx is on the wake side of rx
# Convert so that 0 is ram, 90 is starboard, 180 is wake, 270 is port
azim_tx_from_rx = 360 - azim_tx_from_rx_90port
        
# FM TO BEACON
rx_to_beacon = beacon - rx
rx_to_beacon_lvlh = np.zeros([n, 3])
rx_to_beacon_lvlh_norm =  np.zeros([n, 3])
rx_dot_rx_to_beacon_lvlh = np.zeros([n])
z_orbit = np.array([0, 0, 1])
for i in range(n):
    rx_to_beacon_lvlh[i, :] = ecef_to_lvlh(rx[i, :], vx[i, :], rx_to_beacon[i, :])
    rx_to_beacon_lvlh_norm[i, :] = rx_to_beacon_lvlh[i, :] / np.linalg.norm(rx_to_beacon_lvlh[i, :])
    rx_dot_rx_to_beacon_lvlh[i] = np.dot(z_orbit, rx_to_beacon_lvlh_norm[i, :])
#90. -  acos(cygnss_r_dot_C2G_orbit)/dtor
elev_beacon_from_rx = 90. - np.arccos(rx_dot_rx_to_beacon_lvlh) / deg2rad
elev_beacon_from_rx = -elev_beacon_from_rx # opposite because we want positive elevations to describe the position of the beacon wrt FM
#atan2(C2S_body[1], C2S_body[0])/dtor));
azim_beacon_from_rx_pi_to_pi = np.arctan2(rx_to_beacon_lvlh[:,1], rx_to_beacon_lvlh[:,0]) / deg2rad
azim_beacon_from_rx_90port = np.zeros([n])
for i in range(n):
    if (azim_beacon_from_rx_pi_to_pi[i] < 0):
        azim_beacon_from_rx_90port[i] = 360 + azim_beacon_from_rx_pi_to_pi[i]
    else:
        azim_beacon_from_rx_90port[i] = azim_beacon_from_rx_pi_to_pi[i]

# Until now, azim_beacon_from_rx_90port = 90 means beacon is port of rx, 270 means beacon is startboarf of rx, - means beacon is on the ram side of rx, 180 means beacon is on the wake side of rx
# Convert so that 0 is ram, 90 is starboard, 180 is wake, 270 is port
azim_beacon_from_rx = 360 - azim_beacon_from_rx_90port
azim_beacon_from_rx_body = np.mod(azim_beacon_from_rx + 90., 360)



        
# POSITION OF FM POLAR COORDINATE IN BEACON REFCE FRAME
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']

fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = fontsize_plot)
plt.rc('font', weight='normal') ## make the labels of the ticks in normal

ax = fig.add_subplot(gs[0, 0], projection = 'polar')
ax.plot(target_rx_sat_look_angle_az*deg2rad, 90- target_rx_sat_look_angle_el, linewidth = 4, color = 'blue')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_title('FM polar coord. wrt. beacon', weight = 'normal', fontsize = fontsize_plot*1.1,  y = 1.07) 
ax.tick_params(axis='both', which='major',  width = 2,  labelsize=fontsize_plot, size = 10, right = False)
ax.text(target_rx_sat_look_angle_az[0]*deg2rad, 90- target_rx_sat_look_angle_el[0], date[0] + '\n' + format(target_rx_sat_look_angle_az[0],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(target_rx_sat_look_angle_el[0],'.0f') + u'\N{DEGREE SIGN}' ,fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'right')
ax.text(target_rx_sat_look_angle_az[-1]*deg2rad, 90- target_rx_sat_look_angle_el[-1], date[-1]  + '\n' + format(target_rx_sat_look_angle_az[-1],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(target_rx_sat_look_angle_el[-1],'.0f') + u'\N{DEGREE SIGN}',fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'left')
elev_tick = np.arange(20, 91, 20)
ax.yaxis.set_ticks(elev_tick)
nelev_tick = len(elev_tick)
elev_ticklabel = []
for itick in range(nelev_tick):
    elev_ticklabel.append(format(90-elev_tick[itick], '.0f'))
ax.yaxis.set_ticklabels(elev_ticklabel)
    
#fig.set_size_inches(10, 20)
fig_save_name = filename.replace('.csv', '_fm_wrt_beacon.pdf')
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')


# Table of elveatio vs time every 10 degrees of elvevation
n = len(target_rx_sat_look_angle_el)
i = 1
kel = 10
table_elev = []
table_elev.append([format(np.abs(target_rx_sat_look_angle_el[0]), ".0f")+ u'\N{DEGREE SIGN}',  date[0][11:]])
# eleve from 0 to max, veery 10 degrees
for kel in np.arange(10, 91, 10):
    try:
        ii = np.where(target_rx_sat_look_angle_el >= kel)[0][0]
        table_elev.append([format(target_rx_sat_look_angle_el[ii], ".0f")+ u'\N{DEGREE SIGN}',  date[ii][11:]])
    except IndexError: # means thakel is greater than the max elev
        imax = np.where(target_rx_sat_look_angle_el == np.max(target_rx_sat_look_angle_el))[0][0]
        table_elev.append([format(target_rx_sat_look_angle_el[imax], ".0f")+ u'\N{DEGREE SIGN}',  date[imax][11:]])
        break
kel_max = kel-10
for kel in np.arange(kel_max, 0, -10):
    ii = np.where(target_rx_sat_look_angle_el[imax+1:] <= kel)[0][0]
    table_elev.append([format(target_rx_sat_look_angle_el[imax+1+ii], ".0f")+ u'\N{DEGREE SIGN}',  date[imax+1+ii][11:]])
table_elev.append([format(np.abs(target_rx_sat_look_angle_el[-1]), ".0f")+ u'\N{DEGREE SIGN}',  date[-1][11:]])
print tabulate(table_elev, headers=['FM\nelev.', 'Time\n(UTC)'])
# POSITION OF GPS POLAR COORDINATE IN CYGNSS ORBIT (not body) REFCE FRAME
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']

fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = fontsize_plot)
plt.rc('font', weight='normal') ## make the labels of the ticks in normal

ax = fig.add_subplot(gs[0, 0], projection = 'polar')
ax.plot(azim_tx_from_rx*deg2rad, 90- elev_tx_from_rx, linewidth = 4, color = 'blue')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_title('PRN ' + str((int)(prn[0]))+ ' polar coord. wrt. FM - max elev: ' + format(np.max(elev_tx_from_rx),'.0f') + u'\N{DEGREE SIGN}', weight = 'normal', fontsize = fontsize_plot*1.1,  y = 1.07) 
ax.tick_params(axis='both', which='major',  width = 2,  labelsize=fontsize_plot, size = 10, right = False)
ax.text(azim_tx_from_rx[0]*deg2rad, 90- elev_tx_from_rx[0], date[0] + '\n' + format(azim_tx_from_rx[0],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(elev_tx_from_rx[0],'.0f') + u'\N{DEGREE SIGN}' ,fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'right')
ax.text(azim_tx_from_rx[-1]*deg2rad, 90- elev_tx_from_rx[-1], date[-1]  + '\n' + format(azim_tx_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(elev_tx_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}',fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'left')
#fig.set_size_inches(10, 20)
ax.set_rmax(90)
elev_tick = np.arange(20, 91, 20)
ax.yaxis.set_ticks(elev_tick)
nelev_tick = len(elev_tick)
elev_ticklabel = []
for itick in range(nelev_tick):
    elev_ticklabel.append(format(90-elev_tick[itick], '.0f'))
ax.yaxis.set_ticklabels(elev_ticklabel)


fig_save_name = filename.replace('.csv', '_gps_wrt_fm.pdf')
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  





# POSITION OF BEACON POLAR COORDINATE IN CYGNSS BODY (not orbit) REFCE FRAME
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']

fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = fontsize_plot)
plt.rc('font', weight='normal') ## make the labels of the ticks in normal

ax = fig.add_subplot(gs[0, 0], projection = 'polar')
ax.plot(azim_beacon_from_rx_body*deg2rad, 90- elev_beacon_from_rx, linewidth = 4, color = 'blue')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_title('Beacon polar coord. wrt. FM body plane (0' + u'\N{DEGREE SIGN} ram, 90'+ u'\N{DEGREE SIGN} starboard)', weight = 'normal', fontsize = fontsize_plot*1.1,  y = 1.07) 
ax.tick_params(axis='both', which='major',  width = 2,  labelsize=fontsize_plot, size = 10, right = False)
ax.text(azim_beacon_from_rx_body[0]*deg2rad, 90- elev_beacon_from_rx[0], date[0] + '\n' + format(azim_beacon_from_rx_body[0],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(elev_beacon_from_rx[0],'.0f') + u'\N{DEGREE SIGN}' ,fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'right')
ax.text(azim_beacon_from_rx_body[-1]*deg2rad, 90- elev_beacon_from_rx[-1], date[-1]  + '\n' + format(azim_beacon_from_rx_body[-1],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(elev_beacon_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}',fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'left')
#fig.set_size_inches(10, 20)
ax.set_rmax(90)
elev_tick = np.arange(20, 91, 20)
ax.yaxis.set_ticks(elev_tick)
nelev_tick = len(elev_tick)
elev_ticklabel = []
for itick in range(nelev_tick):
    elev_ticklabel.append(format(90-elev_tick[itick], '.0f'))
ax.yaxis.set_ticklabels(elev_ticklabel)


# ADD THE ANTENNA MAPS
## read the antenna map
filename_gain = 'CygnssAntenna/merged_ant_1_starboard_ddmi_v1_with_ant_1_port_ddmi_v1_test_elevOpposite.bin'
extension = filename_gain.split('.')[-1]
file_gain = open(filename_gain, "rb")
x1 = unpack('i', file_gain.read(4))[0] # don t know what this is ....
x2 = unpack('i', file_gain.read(4))[0] # don t know what this is ....
numAz = unpack('i', file_gain.read(4))[0]
numEl = unpack('i', file_gain.read(4))[0]
az_start_deg = unpack('d', file_gain.read(8))[0]
el_start_deg = unpack('d', file_gain.read(8))[0]
az_inc_deg = unpack('d', file_gain.read(8))[0]
el_inc_deg = unpack('d', file_gain.read(8))[0]
#!!!! for some reason, a negative sign is needed for this file
el_inc_deg = -el_inc_deg
az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
az_deg = np.linspace(az_start_deg, az_stop_deg, numAz)
el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)
gain = np.zeros([numEl,numAz])
for iaz in range(numAz):
    for iel in range(numEl):
        gain[iel, iaz] = unpack('d', file_gain.read(8))[0]
file_gain.close

## ADD THE MAP TO THE POLAR COORD PLOT
x = np.append(az_deg, 360)*deg2rad # need to add 360 so that X is one more dimension than  Z (otherwise the last value of Z is ignored)
y = 90 - np.append(el_deg, 0) # need to add 0 so that Y is one more dimension than  Z (otherwise the last value of Z is ignored)
X, Y = np.meshgrid(x, y)
extension = filename_gain.split('.')[-1]
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

#    if extension == 'bin':
if el_start_deg < el_stop_deg: # file from low to high elevation
    # top of graph is first row of Z,
    # which should be highest elevation (since y axis is
    # low (bottom) to high (top) elevation) but since
    # el_start_deg < el_stop_deg then first row of gain
    # is low elevation so invert it so that first row of Z is
    # highest elevation
    Z = np.zeros([numEl, numAz])
    for iel in range(numEl):
        Z[iel, :] =  gain[numEl-1-iel, :]
else: # first row of gain is highest elevation
    Z = gain
nr, nc = Z.shape

Z = np.ma.array(Z)
max_gain = 15
CS1 = ax.pcolormesh(X, Y, Z, cmap = 'jet',
                    vmin = 0, vmax = max_gain, alpha = 0.3)                    
cbar = plt.colorbar(CS1, ax = ax, pad = 0.075)
cbar.ax.set_ylabel('RCG', fontsize = fontsize_plot, weight = 'normal')

fig_save_name = filename.replace('.csv', '_beacon_wrt_fm_body.pdf')
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


sys.exit()

# POSITION OF BEACON POLAR COORDINATE IN CYGNSS ORBIT (not body) REFCE FRAME
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']

fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.12, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = fontsize_plot)
plt.rc('font', weight='normal') ## make the labels of the ticks in normal

ax = fig.add_subplot(gs[0, 0], projection = 'polar')
ax.plot(azim_beacon_from_rx*deg2rad, 90- elev_beacon_from_rx, linewidth = 4, color = 'blue')
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_title('Beacon polar coord. wrt. FM orbit plane (0' + u'\N{DEGREE SIGN} ram, 90'+ u'\N{DEGREE SIGN} starboard)', weight = 'normal', fontsize = fontsize_plot*1.1,  y = 1.07) 
ax.tick_params(axis='both', which='major',  width = 2,  labelsize=fontsize_plot, size = 10, right = False)
ax.text(azim_beacon_from_rx[0]*deg2rad, 90- elev_beacon_from_rx[0], date[0] + '\n' + format(azim_beacon_from_rx[0],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(elev_beacon_from_rx[0],'.0f') + u'\N{DEGREE SIGN}' ,fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'right')
ax.text(azim_beacon_from_rx[-1]*deg2rad, 90- elev_beacon_from_rx[-1], date[-1]  + '\n' + format(azim_beacon_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(elev_beacon_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}',fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'left')
#fig.set_size_inches(10, 20)
ax.set_rmax(90)
elev_tick = np.arange(20, 91, 20)
ax.yaxis.set_ticks(elev_tick)
nelev_tick = len(elev_tick)
elev_ticklabel = []
for itick in range(nelev_tick):
    elev_ticklabel.append(format(90-elev_tick[itick], '.0f'))
ax.yaxis.set_ticklabels(elev_ticklabel)


fig_save_name = filename.replace('.csv', '_beacon_wrt_fm_orbit.pdf')
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

