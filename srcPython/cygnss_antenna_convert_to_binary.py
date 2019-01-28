# This script takes as input the antenna in txt or agm format and convert it
# into a binary file in the format readable by tds-bop.exe (beacon software)
# INPUTS:
# - filename_gain_in (filename of antenna gain file to convert to binary file)
# - res_map: resolution of the map (set to 'coarse' or 'fine')
# ASSUMPTIONS:
# - for the a description of the format of the input file, read the header of
# the script cygnss_antenna_pattern.py
# - the coarse on-board antennas (.agm) vary in elevation from 0 to 90 (step +5)
#  and in azimuth from -180 to 180 (step +15, -180 to 0 is port, 0 to 180 is
# starboard). However, it looks like tds-bop expects the coarse antennas (.bin)
# to vary in elevation from 90 to 0 (step -5) and in azimuth from 0 to 360
# (step +15, 0 to 180 is starboard, 180 to 360 is port). Therefore, for coarse
# files (i.e., res_map is 'coarse'), this script sets el_start_deg to 90,
# az_start_deg to 0, and el_inc_deg to -5 (actually el_inc_deg is kept at +5,
# see next point)
# - if the resolution res_map is set to "coarse", then, for an unknown
# reason, we need to set el_inc_deg to its opposite value. !!!!!! todo: we
# actually need to check that this statement is true

# old assumptions (not valid anymore):
# - when writing this function, I assumed that the format of the output file
# (i.e., the binary file) didn't have to follow a particular rule, as long as 
# the values of el_start_deg, az_start_deg, numEl, numAz, el_inc_deg, and
# az_inc_deg were correctly indicated in the header. !!!! todo: But for example,
# I'm not # sure either that tds-bop correctly handles azimuths from -180 to 0,
# and tds-bop might expect that elevations go from 90 to 0 (not 0 to 90)
# - however, if the resolution res_map is set to "coarse", then, for an unknown
# reason, we need to set el_inc_deg to its opposite value. !!!!!! todo: we
# actually need to check that this statement is true



# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
res_map = 'fine' # 'coarse' or 'fine'
filename_gain_in = '/Users/cbv/work/spockOut/beacon/ant_data/merged_ant_1_starboard_v6_with_ant_1_port_v6.txt'

# Coarse file names
# '/Users/cbv/cspice/data/ant_1_starboard_ddmi_v1.agm'
# '/Users/cbv/cspice/data/ant_1_port_ddmi_v1.agm'
# '/Users/cbv/cspice/data/merged_ant_1_starboard_ddmi_v1_with_ant_1_port_ddmi_v1.agm'

# Fine file names
# '/Users/cbv/work/spockOut/beacon/ant_data/ant_1_starboard_v6.txt'
# '/Users/cbv/work/spockOut/beacon/ant_data/ant_1_port_v6.txt'
# '/Users/cbv/work/spockOut/beacon/ant_data/merged_ant_1_starboard_v6_with_ant_1_port_v6.txt'


# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import numpy as np
from array import array
import ipdb

if ((res_map != 'coarse') & (res_map != 'fine')):
    sys.exit("!***!\n res_map must be set to 'coarse' or 'fine'\
            \n!***!")
                    
print 'Map resolution: ' + res_map

if res_map == 'coarse': # Coarse map parameters
    el_start_deg = 0; az_start_deg = -180
    numEl = 18; numAz = 24
    el_inc_deg = 5.; az_inc_deg = 15.
else: # Fine map parameters
    el_start_deg = 90; az_start_deg = 0
    numEl = 901; numAz = 3601
    el_inc_deg = -0.1; az_inc_deg = 0.1



file_gain_in = open(filename_gain_in)
read_f = file_gain_in.readlines()
el_start_deg = np.float(el_start_deg); az_start_deg = np.float(az_start_deg)
el_inc_deg = np.float(el_inc_deg); az_inc_deg = np.float(az_inc_deg)
numEl = (int)(numEl); numAz = (int)(numAz)

gain = np.zeros([numEl, numAz])
filename_gain_out = filename_gain_in.split('.')[0] + '_test_elevOpposite.bin'
file_gain_out = open(filename_gain_out, "wb")
# x1 x2 numAz numEl ('i')
x1 = 0; x2 = 0 # not used in cygnss_antenna_pattern.py
float_array = array('i', [x1, x2, numAz, numEl])
float_array.tofile(file_gain_out)
# az_start_deg el_start_deg az_inc_deg el_inc_deg ('d')
#
# if the resolution res_map is set to "coarse", then, for an unknown
# reason, we need to set el_inc_deg to its opposite value. !!!!!! todo: we
# actually need to check that this statement is true
if res_map == 'coarse':
    el_start_deg = 90. # although el_start_deg is 0 in .agm, el_start_deg
    # needs to be 90 in binary
    az_start_deg = 0 # although az_start_deg is -180 in .agm, az_start_deg
    # needs to be 0 in binary
float_array = array('d', [az_start_deg, el_start_deg, az_inc_deg, el_inc_deg])

float_array.tofile(file_gain_out)
if res_map == 'coarse': # in .agm coarse files, rows are elevations, columns are
    # azimuth
    gain_az_not_corrected = np.zeros([numEl, numAz])
    for iel in range(numEl):
        for iaz in range(numAz):
            # numEl-1-iel because elev in .agm goes from 0 to 90 but
            # in binary it elev needs to go from 90 to 0
            gain_az_not_corrected[numEl-1-iel, iaz] =read_f[iel].split(',')[iaz]
    # transformation: in .agm azim goes from -180 to 180 (-180 to 0 is port,
    # 180 to 0 is starboard) but in binary azim needs to go from 0 to 360
    # (0 to 180 is starboard, 180 to 360 is port)
    if numAz != 24:
        sys.exit('The number of azimuth bins should be 24 for the coarse file.')
    gain[:, 0:numAz/2] = gain_az_not_corrected[:, numAz/2:]
    gain[:, numAz/2:] = gain_az_not_corrected[:, 0:numAz/2]
    
else: # in .txt fine files, rows are azimuths, columns are elevations
    nheader = 2
    for iaz in range(numAz):
        print iaz, numAz-1
        for iel in range(numEl):
            gain[iel, iaz] = read_f[iaz+nheader].split(',')[iel]
    
for iaz in range(numAz):
    for iel in range(numEl-1,-1,-1):#!!!!!! should be for iel in range(numEl):
        float_array = array('d', [gain[iel, iaz]])
        float_array.tofile(file_gain_out)

file_gain_out.close()
