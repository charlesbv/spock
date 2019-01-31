# This script creates gaussian nadir pattern with no azimuth dependence
# and a set half power beam width (e.g. 30 deg). It creates the fine
# and coarse antennas, both in binary format following the format of
# tds-bop (see the header of the script cygnss_antenna_pattern.py for
# a complete description of the format


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
res_map = 'coarse' #resolution of the map (set to 'coarse' or 'fine')
fwhm_half = 30. # in degrees. Full Width at Half Maximum / 2 
max_gain = 15
# fwhm represents the elevation at which the
# gain is equal to half its max value (obtained at elevation = 0 deg)
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from struct import *
import sys
import ipdb
from array import array

def gaussian(x, sig, max_g):
    expo = np.exp(-x**2/(2*sig**2))
    f = expo*max_g
    return f


if ((res_map != 'coarse') & (res_map != 'fine')):
    sys.exit("!***!\n res_map must be set to 'coarse' or 'fine'\
            \n!***!")
                    
print 'Map resolution: ' + res_map

fwhm = fwhm_half * 2

sigma = fwhm / (2*np.sqrt(2*np.log(2))) # std deviation

if res_map == 'coarse': # Coarse map parameters
    el_start_deg = 0; az_start_deg = -180
    numEl = 18; numAz = 24
    el_inc_deg = 5.; az_inc_deg = 15.
else: # Fine map parameters
    el_start_deg = 90; az_start_deg = 0
    numEl = 901; numAz = 3601
    el_inc_deg = -0.1; az_inc_deg = 0.1


el_start_deg = np.float(el_start_deg); az_start_deg = np.float(az_start_deg)
el_inc_deg = np.float(el_inc_deg); az_inc_deg = np.float(az_inc_deg)
numEl = (int)(numEl); numAz = (int)(numAz)

gain = np.zeros([numEl, numAz])
filename_gain_out = 'gaussian_' + res_map + '.bin'
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
    el_inc_deg = -el_inc_deg # see header of cygnss_antenna_pattern.py
    el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
    el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)

    gain_az_not_corrected = np.zeros([numEl, numAz])
    for iel in range(numEl):
        gain_az_not_corrected[iel, :] = gaussian(el_deg[iel], sigma,
                                                         max_gain)
    # transformation: in .agm azim goes from -180 to 180 (-180 to 0 is port,
    # 180 to 0 is starboard) but in binary azim needs to go from 0 to 360
    # (0 to 180 is starboard, 180 to 360 is port)
    if numAz != 24:
        sys.exit('The number of azimuth bins should be 24 for the coarse file.')
    gain[:, 0:numAz/2] = gain_az_not_corrected[:, numAz/2:]
    gain[:, numAz/2:] = gain_az_not_corrected[:, 0:numAz/2]
    
else: # in .txt fine files, rows are azimuths, columns are elevations
    nheader = 2
    el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
    el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)

    for iel in range(numEl):
        gain[iel, :] = gaussian(el_deg[iel], sigma, max_gain)
    
for iaz in range(numAz):
    for iel in range(numEl-1,-1,-1):#!!!!!! should be for iel in range(numEl):
        float_array = array('d', [gain[iel, iaz]])
        float_array.tofile(file_gain_out)

file_gain_out.close()
