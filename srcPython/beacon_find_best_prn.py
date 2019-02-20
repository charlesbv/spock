# For a given overpass of CYGNSS over the beacon station in Holloman, NM, this script selects the opitmum specular point, based on certain criteria (see in code what those criterias are)
# Inputs:
# - prefix_pass = prefix of pass to look at. Format: 'pass_N' where N is the overpass to manually select
# - csv_dir = directory of the csv files output by tds-bop


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
prefix_pass = 'pass_9'
csv_dir = 'outputCygnssMar'
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import glob
from beacon_read_csv import *
import numpy as np
if csv_dir[-1] != '/':
    csv_dir = csv_dir + '/'
deg2rad = np.pi/180
filename_list = glob.glob(csv_dir + prefix_pass + '*csv')
nfile = len(filename_list)
elev_tx_from_rx = []
max_elev_tx_from_rx = np.zeros([nfile])
for ifile in range(nfile):
    # Read the cvs file
    filename = filename_list[ifile]
    date, prn, target_lat, target_lon, target_alt, target_ecef_x, target_ecef_y, target_ecef_z, target_rx_sat_look_angle_az, target_rx_sat_look_angle_el, target_rx_sat_range, sp_lat, sp_lon, sp_ecef_pos_x, sp_ecef_pos_y, sp_ecef_pos_z, sp_gain, rx_sub_sat_lat, rx_sub_sat_lon, rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z, rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z, tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z, rx_power = beacon_read_csv(filename)
    # Compute the elevation of the GPS in the polar coord of CYGNSS
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
    elev_tx_from_rx.append( 90. - np.arccos(rx_dot_rx_to_tx_lvlh) / deg2rad )
    max_elev_tx_from_rx[ifile] = np.max(elev_tx_from_rx[-1])

