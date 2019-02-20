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
deg2rad = np.pi/180
filename = 'outputCygnssOct/FM08/pass_5_PRN_21.csv'

date, prn, target_lat, target_lon, target_alt, target_ecef_x, target_ecef_y, target_ecef_z, target_rx_sat_look_angle_az, target_rx_sat_look_angle_el, target_rx_sat_range, sp_lat, sp_lon, sp_ecef_pos_x, sp_ecef_pos_y, sp_ecef_pos_z, sp_gain, rx_sub_sat_lat, rx_sub_sat_lon, rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z, rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z, tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z, rx_power = beacon_read_csv(filename)

rx_to_tx_x = tx_sat_ecef_pos_x - rx_sat_ecef_pos_x
n = len(rx_to_tx_x)
rx = np.array( [rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z]).transpose()
vx = np.array( [rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z]).transpose()
tx = np.array( [tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z]).transpose()
rx_to_tx = tx - rx
rx_to_tx_lvlh = np.zeros([n, 3])
rx_mag = np.linalg.norm(rx, axis = 1)
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
azim_tx_from_rx = np.zeros([n])
for i in range(n):
    if (azim_tx_from_rx_pi_to_pi[i] < 0):
        azim_tx_from_rx[i] = 360 + azim_tx_from_rx_pi_to_pi[i]
    else:
        azim_tx_from_rx[i] = azim_tx_from_rx_pi_to_pi[i]

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
# POSITION OF GPS POLAR COORDINATE IN CYGNSS REFCE FRAME
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
