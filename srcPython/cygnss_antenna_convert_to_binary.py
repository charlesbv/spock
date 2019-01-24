# This script takes as input the antenna in txt or agm format and convert it
# into a binary file in the format readable by tds-bop.exe (beacon software)
# INPUTS:
# filename_gain_in (filename of the antenna gain file to convert to binary file)
# - el_start_deg (elevation start), az_start_deg (azimuth start)
# - el_stop_deg (elevation stop), az_stop_deg (azimuth stop)
# ASSUMPTIONS:
# 1- the file is n_el rows * n_az columns, no header
# 2- each element of a line (azimuth) is separated by a comma
# 3- if fthe format of the input file is different than from 1- and 2-, adpat
# the script accordingly

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
filename_gain_in = '/Users/cbv/cspice/data/ant_1_starboard_ddmi_v1.agm'
el_start_deg = 0; az_start_deg = -180
el_stop_deg = 90; az_stop_deg = 180
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import numpy as np
from array import array


file_gain_in = open(filename_gain_in)
read_f = file_gain_in.readlines()
n_el = len(read_f) # assumes no header
n_az = len(read_f[0].split(',')) # assumes azimuths are comma separated
el_start_deg = np.float(el_start_deg); az_start_deg = np.float(az_start_deg)
el_stop_deg = np.float(el_stop_deg); az_stop_deg = np.float(az_stop_deg)
el_inc_deg = (el_stop_deg - el_start_deg)/n_el
az_inc_deg = (az_stop_deg - az_start_deg)/n_az

gain = np.zeros([n_el, n_az])
filename_gain_out = filename_gain_in.split('.')[0] + '_test.bin'
file_gain_out = open(filename_gain_out, "wb")
# x1 x2 numAz numEl ('i')
x1 = 0; x2 = 0 # not used in cygnss_antenna_pattern.py
numAz = n_az; numEl = n_el # be consistent with cygnss_antenna_pattern.py
float_array = array('i', [x1, x2, numAz, numEl])
float_array.tofile(file_gain_out)
# az_start_deg el_start_deg az_inc_deg el_inc_deg ('d')
float_array = array('d', [az_start_deg, el_start_deg, az_inc_deg, el_inc_deg])
float_array.tofile(file_gain_out)
for iel in range(n_el):
    for iaz in range(n_az):
        gain[iel, iaz] = read_f[iel].split(',')[iaz]
for iaz in range(n_az):
    for iel in range(n_el):
        float_array = array('d', [gain[iel, iaz]])
        float_array.tofile(file_gain_out)

file_gain_out.close()
