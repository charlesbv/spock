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
satbop_output_dir = './' # directory where the csv output files are
satbop_pass = 1 # pass to look at (integer or string)
attitude = [0, 0, -90] # attitude of the FM: [pitch, roll, yaw] (in deg) 
order_rotation = [3,2,1] # rder_rotation order 1 means firt you do this rotation. for example: [2,1,3] -> first you roll then pitch then yaw

# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import glob
import ipdb
from datetime import datetime, timedelta
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
date_entire = []
nseconds = (int)((last_date - first_date).total_seconds()) + 1
for isecond in range(nseconds):
    date_entire.append(datetime.strftime(first_date + timedelta(seconds = isecond), "%Y-%m-%dT%H:%M:%S"))
date_entire = np.array(date_entire)

for isecond in range(nseconds):
    date_now = date_entire[isecond]
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
            prn_this_time.append(prn[iprn])
    # in this list, select the PRNs with 4 highest FOM (with a few other small tricks that the onboard algo does)
    select_highest_foms(elev_this_time, azim_this_time, elev_not_int_this_time, azim_not_int_this_time, elev_gps_from_cyg_body_this_time, prn_this_time, fom_map)
