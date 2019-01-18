# This script compares the specular point output of SpOCK
# with the specular point output of tds-bop
# for the planning of the CYGNSS beacon campaign


import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import numpy as np

from cygnss_read_spock_spec_bin import *
from read_input_file import *
from find_in_read_input_order_variables import *

input_filename = 'spock.txt'
cygfm = 3

var_in, var_in_order = read_input_file(input_filename)
dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, \
                                                            'dt_output')]; 
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_name_list')];
gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order,\
                                                                'gps_name')];
cygfm_to_spock_nb = [4,3,8,2,1,6,7,5]
# ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
isc =  cygfm_to_spock_nb[cygfm-1] - 1
spec_spock_filename = output_file_path_list[isc] + "specular_" + \
    output_file_name_list[isc]

data_spec =cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'\
), gps_name_list_spock, dt_spock_output, 1) 
date_spock = data_spec[0]; lon_spock = data_spec[1]; lat_spock = data_spec[2]
gain_spock = data_spec[3]; gps_spock = data_spec[4]
normpower_spock = data_spec[5]; x_cyg_spock = data_spec[6]
y_cyg_spock = data_spec[7]
z_cyg_spock = data_spec[8]; x_gps_spock = data_spec[9]
y_gps_spock = data_spec[10]; z_gps_spock = data_spec[11]
x_spec_spock = data_spec[12]; y_spec_spock = data_spec[13]
z_spec_spock = data_spec[14]; nb_spec_spock = data_spec[15]
el_spec_spock = data_spec[16]; az_spec_spock = data_spec[17]
el_gps_from_cyg_spock = data_spec[18];  el_spec_not_int_spock = data_spec[19]
az_spec_not_int_spock = data_spec[20]

