# This script reads the csv files output by tds-bop.exe
import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import csv
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from ecef_to_lvlh import *
deg2rad = np.pi/180
filename = 'outputCygnssMar/pass_9_PRN_31.csv'
irow = -1
date = []
with open(filename) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    yy = []; mm = []; dd = []; hr = []; minute = []; sec = []; elapsed_sec = [];
    prn = []; target_lat = []; target_lon = []; target_alt = [];
    target_ecef_x = []; target_ecef_y = []; target_ecef_z = [];
    target_rx_sat_look_angle_az = []; target_rx_sat_look_angle_el = [];
    target_rx_sat_range = []
    sp_lat = []; sp_lon = []; sp_ecef_pos_x = []; sp_ecef_pos_y = [];
    sp_ecef_pos_z = []; sp_gain = []; rx_sub_sat_lat = []; rx_sub_sat_lon = [];
    rx_sat_ecef_pos_x = []; rx_sat_ecef_pos_y = []; rx_sat_ecef_pos_z = [];
    rx_sat_ecef_vel_x = []; rx_sat_ecef_vel_y = []; rx_sat_ecef_vel_z = [];
    tx_sat_ecef_pos_x = []; tx_sat_ecef_pos_y = []; tx_sat_ecef_pos_z = [];
    rx_power= [];
    for row in readCSV:
        irow = irow + 1
        if irow != 0:
            date.append(row[0] + '-' + row[1].zfill(2) + '-' + row[2].zfill(2) + \
                        'T' + row[3].zfill(2) + ':' + row[4].zfill(2) + ':' + row[5].zfill(2))
            yy.append( np.float(row[0]) )
            mm.append( np.float(row[1]) )
            dd.append( np.float(row[2]) )
            hr.append( np.float(row[3]) )
            minute.append( np.float(row[4]) )
            sec.append( np.float(row[5]) )
            elapsed_sec.append( np.float(row[6]) )
            prn.append( np.float(row[7]) )
            target_lat.append( np.float(row[8]) )
            target_lon.append( np.float(row[9]) )
            target_alt.append( np.float(row[10]) )
            target_ecef_x.append( np.float(row[11]) )
            target_ecef_y.append( np.float(row[12]) )
            target_ecef_z.append( np.float(row[13]) )
            target_rx_sat_look_angle_az.append( np.float(row[14]) )
            target_rx_sat_look_angle_el.append( np.float(row[15]) )
            target_rx_sat_range.append( np.float(row[16]) )
            sp_lat.append( np.float(row[19]) )
            sp_lon.append( np.float(row[20]) )
            sp_ecef_pos_x.append( np.float(row[21]) )
            sp_ecef_pos_y.append( np.float(row[22]) )
            sp_ecef_pos_z.append( np.float(row[23]) )
            sp_gain.append( np.float(row[27]) )
            rx_sub_sat_lat.append( np.float(row[28]) )
            rx_sub_sat_lon.append( np.float(row[29]) )
            rx_sat_ecef_pos_x.append( np.float(row[30]) )
            rx_sat_ecef_pos_y.append( np.float(row[31]) )
            rx_sat_ecef_pos_z.append( np.float(row[32]) )
            rx_sat_ecef_vel_x.append( np.float(row[33]) )
            rx_sat_ecef_vel_y.append( np.float(row[34]) )
            rx_sat_ecef_vel_z.append( np.float(row[35]) )
            tx_sat_ecef_pos_x.append( np.float(row[36]) )
            tx_sat_ecef_pos_y.append( np.float(row[37]) )
            tx_sat_ecef_pos_z.append( np.float(row[38]) )
            rx_power.append( np.float(row[42]) )

yy = np.array(yy); mm = np.array(mm); dd = np.array(dd); hr = np.array(hr);
minute = np.array(minute); sec = np.array(sec); elapsed_sec = np.array(elapsed_sec);
prn = np.array(prn); target_lat = np.array(target_lat);
target_lon = np.array(target_lon); target_alt = np.array(target_alt);
target_ecef_x = np.array(target_ecef_x); target_ecef_y = np.array(target_ecef_y);
target_ecef_z = np.array(target_ecef_z);
target_rx_sat_look_angle_az = np.array(target_rx_sat_look_angle_az);
target_rx_sat_look_angle_el = np.array(target_rx_sat_look_angle_el);
target_rx_sat_range = np.array(target_rx_sat_range);
sp_lat = np.array(sp_lat); sp_lon = np.array(sp_lon);
sp_ecef_pos_x = np.array(sp_ecef_pos_x); sp_ecef_pos_y = np.array(sp_ecef_pos_y);
sp_ecef_pos_z = np.array(sp_ecef_pos_z); sp_gain = np.array(sp_gain);
rx_sub_sat_lat = np.array(rx_sub_sat_lat); rx_sub_sat_lon = np.array(rx_sub_sat_lon);
rx_sat_ecef_pos_x = np.array(rx_sat_ecef_pos_x);
rx_sat_ecef_pos_y = np.array(rx_sat_ecef_pos_y); rx_sat_ecef_pos_z = np.array(rx_sat_ecef_pos_z);
rx_sat_ecef_vel_x = np.array(rx_sat_ecef_vel_x); rx_sat_ecef_vel_y = np.array(rx_sat_ecef_vel_y);
rx_sat_ecef_vel_z = np.array(rx_sat_ecef_vel_z);
tx_sat_ecef_pos_x = np.array(tx_sat_ecef_pos_x); tx_sat_ecef_pos_y = np.array(tx_sat_ecef_pos_y);
tx_sat_ecef_pos_z = np.array(tx_sat_ecef_pos_z); rx_power= np.array(rx_power);
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
ax.text(target_rx_sat_look_angle_az[0]*deg2rad, 90- target_rx_sat_look_angle_el[0], date[0] + '\n' + format(target_rx_sat_look_angle_az[0],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(90-target_rx_sat_look_angle_el[0],'.0f') + u'\N{DEGREE SIGN}' ,fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'right')
ax.text(target_rx_sat_look_angle_az[-1]*deg2rad, 90- target_rx_sat_look_angle_el[-1], date[-1]  + '\n' + format(target_rx_sat_look_angle_az[-1],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(90-target_rx_sat_look_angle_el[-1],'.0f') + u'\N{DEGREE SIGN}',fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'left')
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
ax.set_title('GPS polar coord. wrt. FM', weight = 'normal', fontsize = fontsize_plot*1.1,  y = 1.07) 
ax.tick_params(axis='both', which='major',  width = 2,  labelsize=fontsize_plot, size = 10, right = False)
ax.text(azim_tx_from_rx[0]*deg2rad, 90- elev_tx_from_rx[0], date[0] + '\n' + format(azim_tx_from_rx[0],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(90 - elev_tx_from_rx[0],'.0f') + u'\N{DEGREE SIGN}' ,fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'right')
ax.text(azim_tx_from_rx[-1]*deg2rad, 90- elev_tx_from_rx[-1], date[-1]  + '\n' + format(azim_tx_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}'  + '/' + format(90 - elev_tx_from_rx[-1],'.0f') + u'\N{DEGREE SIGN}',fontsize = fontsize_plot, verticalalignment = 'top', horizontalalignment = 'left')
#fig.set_size_inches(10, 20)
ax.set_rmax(90)
fig_save_name = filename.replace('.csv', '_gps_wrt_fm.pdf')
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
