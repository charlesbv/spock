# This script compares the specular point position and selection algorithm between two runs by SpOCK.
# It was build to compare the SP outputs using a merged antenna gain pattern file and the SP outputs using two antenna gain pattern files (port and starboard, like the onboard algorithm does).
# These tests were done in preparation of the CYGNSS beacon campaign in 2019.
# INPUTS:
# - filename1 = main input file of SpOCK of the 1st run
# - filename2 = main input file of SpOCK of the 2nd run
# ASSUMPTIONS
# - none 


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
filename1 = 'twoAntennas.txt'
filename2 = 'merged.txt'
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import ipdb
import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import numpy as np

from cygnss_read_spock_spec_bin import *
import read_input_file; reload(read_input_file); from read_input_file import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta

cygfm_to_spock_nb = [4,3,8,2,1,6,7,5]

cygfm = 3
isc =  cygfm_to_spock_nb[cygfm-1] - 1

# Read first SpOCK main input file
var_in, var_in_order = read_input_file(filename1)
dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, \
                                                            'dt_output')]; 
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_name_list')];
gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order,\
                                                                'gps_name')];
# Read specular point output of SpOCK from the first run
print 'Reading SP for first run...'
spec_spock_filename = output_file_path_list[isc] + "specular_" + \
    output_file_name_list[isc]
data_spec =cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'\
), gps_name_list_spock, dt_spock_output, 1) 
date_spock1 = data_spec[0]; lon_spock1 = data_spec[1]; lat_spock1 = data_spec[2]
gain_spock1 = data_spec[3]; gps_spock1 = data_spec[4]



# Read second SpOCK main input file
var_in, var_in_order = read_input_file(filename2)
dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, \
                                                            'dt_output')]; 
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_name_list')];
gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order,\
                                                                'gps_name')];
# Read specular point output of SpOCK from the second run
print 'Reading SP for second run...'
spec_spock_filename = output_file_path_list[isc] + "specular_" + \
    output_file_name_list[isc]
data_spec =cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'\
), gps_name_list_spock, dt_spock_output, 1) 
date_spock2 = data_spec[0]; lon_spock2 = data_spec[1]; lat_spock2 = data_spec[2]
gain_spock2 = data_spec[3]; gps_spock2 = data_spec[4]

