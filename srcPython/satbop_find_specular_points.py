# This script has a similar purpose as find_specular_points.c: at a given time for a given FM, it attempts to select the same 4 specular points as the onboard algo. The main difference is that the orbit predictions (CYGNSS and GPS), as well as the specular point positions, are made with sat-bop.exe, not SpOCK.
# Basically sat-bop.exe outputs the information (SP positions, etc) for ALL existing SPs (not only the top 4. This script selects the top 4 at every second of a pass (a pass is an interval of time for which sat-bop.exe created these output files)

# As such, for a given FM, this script:
# 1. considers all csv files for a given pass, output by sat-bop.exe
# 2. reads each csv file, which contains the SP and CYGNSS position -> store the SP and CYGNSS positions
# 3. once all files are read, compute the elevation and azimuth angles of each SP with respect to the FM in the FM body frame of reference
# 4. read the CYGNSS port and starboard antenna FOM maps
# 5. based on these body elevation and azimuths and the antenna FOM maps, determine the FOM for each SP
# 6. selects the SPs that have the 4 highest FOMs (with a few other small tricks that the onboard algo does)
# Steps 5. and 6. are perfored with the function select_highest_foms

# ASSUMPTIONS:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
# - for a given pass, all csv files have the same length



# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
satbop_output_dir = 'FM03_2018-10-31' #'20180926_fm02' # directory where the csv output files are
satbop_pass = 1#5 # pass to look at (integer or string)
attitude = [0, 0, -90] # attitude of the FM: [pitch, roll, yaw] (in deg) 
order_rotation = [3,2,1] # rder_rotation order 1 means firt you do this rotation. for example: [2,1,3] -> first you roll then pitch then yaw
reproduce_onboard_bug_antenna_selection = 1 # if set to 1 then the PRN and antenna selection algorithm reproduces the onboard bug which doesn't properly switch from fore and aft antenna when the FM is yawed by +-90 deg
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import glob
import ipdb
from datetime import datetime, timedelta
import matplotlib
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import compute_T_sc_to_orbit; reload(compute_T_sc_to_orbit); from compute_T_sc_to_orbit import *
from beacon_read_csv import *
from ecef_to_lvlh_nadir import *
import read_cygnss_agm_antennas; reload(read_cygnss_agm_antennas); from read_cygnss_agm_antennas import *
import select_highest_foms; reload(select_highest_foms); from select_highest_foms import *

if satbop_output_dir[-1] != '/':
    satbop_output_dir = satbop_output_dir + '/'
satbop_pass = str(satbop_pass)
# Determine the name of all csv files for this pass
root_filename = 'pass_' + satbop_pass + '_PRN_'
csv_filename_list = glob.glob(satbop_output_dir + root_filename + '*csv')

# Read the starboard (ant #2) and port (ant #3) antenna FOM maps
fom_map = read_cygnss_agm_antennas()

# Compute the matrix to convert from orbit frame to body frame
T_sc_to_orbit = compute_T_sc_to_orbit(attitude, order_rotation)
T_orbit_to_sc = np.transpose(T_sc_to_orbit)

# Read every csv file
nprn = len(csv_filename_list)
date = []; prn = []; azim_not_int = []; azim = []; elev_not_int = []; elev = []; elev_gps_from_cyg_body = []
first_date = datetime.strptime('2100-04-26T16:50:00', "%Y-%m-%dT%H:%M:%S")
last_date = datetime.strptime('1900-04-26T16:50:00', "%Y-%m-%dT%H:%M:%S")
z_body_here = [0,0,-1];
for iprn in range(nprn):
#iprn = 0 #!!!!!!!! remove this line and use the for loop over all prn: for iprn in range(nprn):
    csv_filename = csv_filename_list[iprn]
    prn_from_name = csv_filename.split('PRN_')[1].split('.')[0]
    print 'iprn '  + str(iprn + 1) + ' (PRN' + prn_from_name.zfill(2) + ')'+ ' out of ' + str(nprn) + ' PRNs'
    
    date_here, prn_here, target_lat, target_lon, target_alt, target_ecef_x, target_ecef_y, target_ecef_z, target_rx_sat_look_angle_az, target_rx_sat_look_angle_el, target_rx_sat_range, sp_lat, sp_lon, sp_ecef_pos_x, sp_ecef_pos_y, sp_ecef_pos_z, sp_gain, rx_sub_sat_lat, rx_sub_sat_lon, rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z, rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z, tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z, rx_power = beacon_read_csv(csv_filename)
    
    
    first_date_here = datetime.strptime(date_here[0], "%Y-%m-%dT%H:%M:%S")
    last_date_here = datetime.strptime(date_here[-1], "%Y-%m-%dT%H:%M:%S")
    if first_date_here < first_date:
        first_date = first_date_here
    if last_date_here > last_date:
        last_date = last_date_here

    prn.append(prn_here[0])
        
    # Compute the elevation and azimuth angles of each SP with respect to the FM in the FM body frame of reference
    dtor = np.pi / 180
    rx = np.array( [rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z]).transpose()
    tx = np.array( [tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z]).transpose()
    vx = np.array( [rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z]).transpose()
    sp = np.array( [sp_ecef_pos_x, sp_ecef_pos_y, sp_ecef_pos_z]).transpose()
    n = len(rx_sat_ecef_pos_x)
    azim_not_int_here = np.zeros([n]); azim_here = np.zeros([n]); elev_not_int_here = np.zeros([n]); elev_here = np.zeros([n]); elev_gps_from_cyg_body_here = np.zeros([n])
    C2S = sp - rx
    C2G = tx - rx
    for i in range(n):
        # CYG TO SP
        C2S_lvlh = ecef_to_lvlh_nadir(rx[i, :], vx[i, :], C2S[i, :])
        C2S_lvlh_norm = C2S_lvlh / np.linalg.norm(C2S_lvlh)
        C2S_body = np.matmul(T_orbit_to_sc, C2S_lvlh_norm)
        C2S_body_norm = C2S_body / np.linalg.norm(C2S_body)
        azim_not_int_here[i] = np.arctan2(C2S_body[1], C2S_body[0])/dtor
        azim_here[i] = (int)(round(np.arctan2(C2S_body[1], C2S_body[0])/dtor))
        elev_not_int_here[i] = 90. - np.arcsin(C2S_body_norm[2])/dtor
        elev_here[i] = 90 - (int)(round(np.arcsin(C2S_body_norm[2])/dtor))

        # CYG TO GP
        C2G_orbit = ecef_to_lvlh_nadir(rx[i, :], vx[i, :], C2G[i, :])
        C2G_body = np.matmul(T_orbit_to_sc, C2G_orbit)
        C2G_body_norm = C2G_body / np.linalg.norm(C2G_body)
        cygnss_r_dot_C2G_body =  np.dot(z_body_here, C2G_body_norm);
        elev_gps_from_cyg_body_here[i] = 90. -  np.arccos(cygnss_r_dot_C2G_body)/dtor; # 0 if GPS at horizon from CYGNSS, 90 if GPS right above CYGNSS
    date.append(date_here); azim_not_int.append(azim_not_int_here); azim.append(azim_here.astype(int)); elev_not_int.append(elev_not_int_here); elev.append(elev_here.astype(int)); elev_gps_from_cyg_body.append(elev_gps_from_cyg_body_here)
    
# date of the entire overpass. first_date is the oldest date among all csv files, last_date is the most recent date
date_entire = []; date_datetime_entire = []
nseconds = (int)((last_date - first_date).total_seconds()) + 1
nb_seconds_since_start = np.zeros([nseconds])
for isecond in range(nseconds):
    date_datetime_entire.append(first_date + timedelta(seconds = isecond))
    date_entire.append(datetime.strftime(date_datetime_entire[-1], "%Y-%m-%dT%H:%M:%S"))
    nb_seconds_since_start[isecond] = isecond
date_entire = np.array(date_entire)
date_datetime_entire = np.array(date_datetime_entire)

prn_selected = []; fom_selected = []; which_ant_selected = []; elev_selected = []; azim_selected = []; elev_not_int_selected = []; azim_not_int_selected = []
for isecond in range(nseconds):
    date_now = date_entire[isecond]
    date_now_date = datetime.strptime(date_now, "%Y-%m-%dT%H:%M:%S")
    elev_this_time = []; azim_this_time = []; prn_this_time = []
    elev_not_int_this_time = []; azim_not_int_this_time = [];
    elev_gps_from_cyg_body_this_time = [];
    # figure out the list of PRNs exisisting at this particular time
    for iprn in range(nprn):
        if date_now in date[iprn]:
            itime = date[iprn].index(date_now)
            elev_this_time.append(elev[iprn][itime])
            azim_this_time.append(azim[iprn][itime])
            elev_not_int_this_time.append(elev_not_int[iprn][itime])
            azim_not_int_this_time.append(azim_not_int[iprn][itime])
            elev_gps_from_cyg_body_this_time.append(elev_gps_from_cyg_body[iprn][itime])
            prn_this_time.append((int)(prn[iprn]))
    # in this list, select the PRNs with 4 highest FOM (with a few other small tricks that the onboard algo does)
    if isecond == 0:
        prn_selected_here, fom_selected_here, which_ant_selected_here, elev_selected_here, azim_selected_here, elev_not_int_selected_here, azim_not_int_selected_here = \
            select_highest_foms(elev_this_time, azim_this_time, elev_not_int_this_time, azim_not_int_this_time, elev_gps_from_cyg_body_this_time, prn_this_time, fom_map, [], [], 0) # since this is the first step, we don't reproduce the bug onbaord (that necessarily needs the previous step to be reproduced) to prn_previous_step is set to [].Same for which_ant_previous_step
    else:
        prn_selected_here, fom_selected_here, which_ant_selected_here, elev_selected_here, azim_selected_here, elev_not_int_selected_here, azim_not_int_selected_here = \
            select_highest_foms(elev_this_time, azim_this_time, elev_not_int_this_time, azim_not_int_this_time, elev_gps_from_cyg_body_this_time, prn_this_time, fom_map, prn_selected[-1], which_ant_selected[-1], reproduce_onboard_bug_antenna_selection)

    prn_selected.append(prn_selected_here); fom_selected.append(fom_selected_here); which_ant_selected.append(which_ant_selected_here); elev_selected.append(elev_selected_here); azim_selected.append(azim_selected_here); elev_not_int_selected.append(elev_not_int_selected_here); azim_not_int_selected.append(azim_not_int_selected_here);     

date_interest_start = '2018-10-31T18:43:10' #'2018-09-26T12:04:58'
date_interest_stop = '2018-10-31T18:55:20'#'2018-09-26T12:17:31'

date_interest_start_datetime =  datetime.strptime(date_interest_start, "%Y-%m-%dT%H:%M:%S")
date_interest_stop_datetime =  datetime.strptime(date_interest_stop, "%Y-%m-%dT%H:%M:%S")
itime_start = np.where(date_datetime_entire == date_interest_start_datetime)[0][0]
itime_stop = np.where(date_datetime_entire == date_interest_stop_datetime)[0][0]
inter_dur_sec = (date_interest_stop_datetime - date_interest_start_datetime).total_seconds()

# PLOT PRNS VS TIME
import matplotlib.colors as colors
import matplotlib.patches as mpatches
plot_selected_prn = 1 # if this is set to 1 then 2 PRNs selected for the overpass will also be reported at the bottom of the plot
color_gain = ['grey', 'blue', 'limegreen', 'red']
label_gain = ['0', '2-5', '6-10', '11-15']
handles_arr = []
for icat in range(len(label_gain)):
    handles_arr.append(mpatches.Patch(color=color_gain[icat], label=label_gain[icat]))

marker_ant = [10, 11]
dant_netcdf = [0, 0.15]
dant_spock = [0, -0.15]

# figure out the number of different prn (SpOCK and on-board) during that time intervals
prn_list = []
for itime in range(itime_start, itime_stop):
    for ispec in range(4):
        if ( prn_selected[itime][ispec] in prn_list ) == False:
            prn_list.append(prn_selected[itime][ispec])

prn_list = np.array(prn_list)
nprn = len(prn_list)
prn_list_sort = prn_list[np.argsort(prn_list)]

height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'PRN'
x_label = 'Real time'
ax_title =  date_entire[itime_start].replace('T', ' ') + ' to ' + date_entire[itime_stop][11:19] + ' UTC - sat-bop'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_title(ax_title, weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7,bottom=1,  left=1, right=1)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
gps_satbop_value = []
for itime in range(itime_start, itime_stop):
    gain_spock_now = fom_selected[itime]
    gain_sort_index = np.argsort(-gain_spock_now) # descending order
    for ispec in range(4):
        prn_spock = prn_selected[itime][ispec]
        if which_ant_selected[itime][ispec] == 2: # starboard
            iant_spock = 0
        elif which_ant_selected[itime][ispec] == 3: # starboard
            iant_spock = 1
        prn_spock_value = np.where(prn_list_sort == prn_spock)[0][0]


        xaxis = (nb_seconds_since_start[itime] - nb_seconds_since_start[itime_start])/60.
        if fom_selected[itime][ispec] == 0: # !!!!!!!!!!! if change gain limmits then need to change also label_gain
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[0], s = 20)   
        elif ((fom_selected[itime][ispec] >= 1) & (fom_selected[itime][ispec] <= 5)):
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[1], s = 20)   
        elif ((fom_selected[itime][ispec] >= 6) & (fom_selected[itime][ispec] <= 10)):
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[2], s = 20)   
        elif ((fom_selected[itime][ispec] >= 11) & (fom_selected[itime][ispec] <= 15)):
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[3], s = 20)   
        else:
            print "***! Error: the gain is not between 0 and 15. !***"; sys.exit()

tmax = (nb_seconds_since_start[itime_stop] - nb_seconds_since_start[itime_start])/60./2
tstar = tmax + 130/60.
ax.plot([tmax, tmax], [-0.6, -0.6 + 0.05], linewidth = 2, color = 'k')
ax.plot([tstar, tstar], [-0.6, -0.6 + 0.05], linewidth = 2, color = 'k')

dt_xlabel =  1. # min
xticks = np.arange(0, inter_dur_sec/60.+1, dt_xlabel)
date_list_str = []
date_list = [date_datetime_entire[itime_start] + timedelta(seconds=(x-xticks[0])*60.) for x in xticks]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[11:19] )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot, rotation=60, horizontalalignment = 'center')

ax.yaxis.set_ticks(np.arange(1, nprn+1))
ax.yaxis.set_ticklabels(prn_list_sort, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0)
ax.set_ylim([-0.6, nprn+0.5])
ax.text(tmax, ax.get_ylim()[0], ' Tmax', rotation = 90, fontsize = fontsize_plot, horizontalalignment  = 'center', verticalalignment = 'bottom')
ax.text(tstar, ax.get_ylim()[0], ' Tmax+130s', rotation = 90, fontsize = fontsize_plot, horizontalalignment  = 'center', verticalalignment = 'bottom')
legend = ax.legend( loc='center left',  bbox_to_anchor=(1, 0.5), fontsize = fontsize_plot, handles=handles_arr, ncol=1, frameon=False, title = 'PRN gain')
legend.get_title().set_fontsize(str(fontsize_plot)) 

fig_save_name = '/Users/cbv/satbop_' + satbop_output_dir[:-1] + '_itimeStart' + str(itime_start) + '_itimeStop' + str(itime_stop) + '_prn_select.pdf'
#'testfm0' + str(cygfm)+ '.pdf'#+str(itime_in) + '_score_3d_binom.pdf'#time_diagram_prn_spock_onboard_iday' + str(idate) + '_itimeStart' + str(itime_start) + '_itimeStop' + +str(itime_stop) + '.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
print fig_save_name
