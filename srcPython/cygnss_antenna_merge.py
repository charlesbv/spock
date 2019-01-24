# This script reads in the port and starboard antenna patterns and merge them into one single antenna pattern.
# INPUTS:
# - starboard_antenna_filename: name of the starboard antenna patter file
# - port_antenna_filename: name of the port antenna patter file
# OUTPUT:
# - the merged antenna gain pattern is created as follow:
  # - for each elev/azim bin, the gain is the max gain between the gain of the starboard antenna and the gain of the port antenna for this bin
# ASSUMPTION:
# - the pattern files have the .agm extension and must be the same format as the onboard antenna patter files:
   # - 18 rows (elev) by 24 col (azim), ie same format as the onboard antenna pattern files
   # - elev from 0 to 90 (every 5 deg), azim from -180 to 180 (every 15 deg)

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
starboard_antenna_filename = '/Users/cbv/cspice/data/ant_1_starboard_ddmi_v1.agm'
port_antenna_filename = '/Users/cbv/cspice/data/ant_1_port_ddmi_v1.agm'
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import numpy as np

# Read starboard antenna gain
starboard_antenna_file = open(starboard_antenna_filename)
read_starboard = starboard_antenna_file.readlines()
nb_theta = len(read_starboard)
nb_phi = len(read_starboard[0].split(','))
theta_max = 90. 
dtheta = (int)( theta_max/nb_theta )
if dtheta != 5.:
    sys.exit('d_elev for starboard antenna should be 5 deg.')
theta_arr = np.arange(0, theta_max, dtheta)
phi_max = 180. 
dphi = (int)( phi_max * 2/nb_phi )
if dphi != 15.:
    sys.exit('d_azim for starboard antenna should be 15 deg.')
phi_arr = np.arange(-180, phi_max, dphi)
starboard_gain = np.zeros([nb_theta, nb_phi])
for itheta in range(nb_theta):
    for iphi in range(nb_phi):
        starboard_gain[itheta, iphi] = np.float( read_starboard[itheta].split(',')[iphi] )
starboard_antenna_file.close()

# Read starboard antenna gain
port_antenna_file = open(port_antenna_filename)
read_port = port_antenna_file.readlines()
nb_theta = len(read_port)
nb_phi = len(read_port[0].split(','))
theta_max = 90. 
dtheta = (int)( theta_max/nb_theta )
if dtheta != 5.:
    sys.exit('d_elev for port antenna should be 5 deg.')
theta_arr = np.arange(0, theta_max, dtheta)
phi_max = 180. 
dphi = (int)( phi_max * 2/nb_phi )
if dphi != 15.:
    sys.exit('d_azim for port antenna should be 15 deg.')
phi_arr = np.arange(-180, phi_max, dphi)
port_gain = np.zeros([nb_theta, nb_phi])
for itheta in range(nb_theta):
    for iphi in range(nb_phi):
        port_gain[itheta, iphi] = np.float( read_port[itheta].split(',')[iphi] )
port_antenna_file.close()

# Calculate the merged gain
merged_gain =  np.zeros([nb_theta, nb_phi])

starboard_antenna_filename_no_path = starboard_antenna_filename.split('/')[-1]
path_to_starboard_file = starboard_antenna_filename.replace(starboard_antenna_filename_no_path, '')
port_antenna_filename_no_path = port_antenna_filename.split('/')[-1]
merged_antenna_filename = path_to_starboard_file + 'merged_' + starboard_antenna_filename_no_path.replace('.agm', '_with_') + port_antenna_filename_no_path
merged_antenna_file = open(merged_antenna_filename, "w")
for itheta in range(nb_theta):
    row = ''
    for iphi in range(nb_phi):
        if starboard_gain[itheta, iphi] >= port_gain[itheta, iphi]:
            merged_gain[itheta, iphi] = starboard_gain[itheta, iphi]
        else:
            merged_gain[itheta, iphi] =	port_gain[itheta, iphi]
        if iphi != nb_phi-1:
            row = row + str((int)(merged_gain[itheta, iphi])) + ','
        else:
            row = row + str((int)(merged_gain[itheta, iphi]))
    print >> merged_antenna_file, row
merged_antenna_file.close()
