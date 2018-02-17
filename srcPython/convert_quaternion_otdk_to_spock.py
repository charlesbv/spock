# This scirpt converts a .a file created by ODTK to a .txt file with a format compatible with SpOCK. Attitude is specified as quaternions
# ASSUMPTIONS:
# - see section "PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT"
# - SpOCK requires the quaternions to be expressed in the SPICE format (see ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html) -> if the input file to convert is in the engineering format, then set the variable file_in_engineering_convention to 1 so that this scripts makes the conversion from engineering to SPICE format
import numpy as np
import sys
sys.path.append("../../kalman/spock_development_new_structure_kalman_dev/srcPython")
import os
from read_input_file import *
from read_output_file import *
from spock_main_input import *
from orbit_average import *
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
filename_odtk = "cyg_data/2F_20170828_214648_STKdefpred_v001.a" # ODTF file input to convert into SpOCK's format
file_in_engineering_convention = 1 # set to 1 if filename_odtk shows quaternion under the engineering format (see ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html). set to 0 if filename_odtk shows quaternion under the SPICE format (see ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/q2m_c.html)
filename_spock = "q_spice.txt" # name of the converted quaternion file -> input of SpOCK
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

file_odtk = open(filename_odtk) 
read_file_odtk = file_odtk.readlines()
iline = 0
found_epoch = 0
while (found_epoch == 0):
    if ('YYDDD' in read_file_odtk[iline]):
        epochyyddd = np.float(read_file_odtk[iline].split(":")[1])
        found_epoch = 1
    iline = iline + 1

found_hdr = 0
while (found_hdr == 0):
    if (len(read_file_odtk[iline].split()) > 0):
        if read_file_odtk[iline].split()[0] == "AttitudeTimeQuatAngVels":
            found_hdr = 1
    iline = iline + 1
nb_header = iline

n = len(read_file_odtk) - nb_header - 2
date_odtk = []
epochyyddd_no_decimal = (int)(epochyyddd)
epochyyddd_decimal_only = epochyyddd - epochyyddd_no_decimal
epochyyddd_date = datetime.strptime( str(epochyyddd_no_decimal), "%y%j" ) + timedelta( hours = epochyyddd_decimal_only*24 )
q_odtk = np.zeros([n, 4]) # quaternion read in odtk file 

file_spock = open(filename_spock, "w+")
print 
print >> file_spock, "#BEGINNINGOFHEADER"
print >> file_spock, "#ENDOFHEADER"

for iline in range(n):
    date_odtk_temp = np.float( read_file_odtk[iline+nb_header].split()[0] )
    date_odtk.append( datetime.strftime( epochyyddd_date + timedelta( seconds = date_odtk_temp ) , "%Y-%m-%dT%H:%M:%S%.%f")[:-3] )
    if file_in_engineering_convention == 1:
        q_odtk[iline, 0] = np.float( read_file_odtk[iline+nb_header].split()[1] )
        q_odtk[iline, 1] = np.float( read_file_odtk[iline+nb_header].split()[2] )
        q_odtk[iline, 2] = np.float( read_file_odtk[iline+nb_header].split()[3] )
        q_odtk[iline, 3] = np.float( read_file_odtk[iline+nb_header].split()[4] )
        # Conversion engineering to SPICE (SpOCK requires SPICE format)
        q_odtk[iline, 0] = q_odtk[iline, 3]
        q_odtk[iline, 1] = -q_odtk[iline, 0]
        q_odtk[iline, 2] = -q_odtk[iline, 1]
        q_odtk[iline, 3] = -q_odtk[iline, 2]

    else:
        q_odtk[iline, 0] = np.float( read_file_odtk[iline+nb_header].split()[1] )
        q_odtk[iline, 1] = np.float( read_file_odtk[iline+nb_header].split()[2] )
        q_odtk[iline, 2] = np.float( read_file_odtk[iline+nb_header].split()[3] )
        q_odtk[iline, 3] = np.float( read_file_odtk[iline+nb_header].split()[4] )
    print >> file_spock, date_odtk[-1], format(q_odtk[iline,0], '.14e'), format(q_odtk[iline,1], '.14e'), format(q_odtk[iline,2], '.14e'), format(q_odtk[iline,3], '.14e')
    
print >> file_spock, "#ENDOFFILE"
file_spock.close()

