# This script answers Chris's questions in his email on 08/12/2019

# ALGORITHM
from datetime import datetime, timedelta
import sys
import os
import ipdb
sys.path.append("/Users/cbv/work/spock/srcPython")
from read_input_file import *
from cygnss_read_spock_spec_bin import *
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt


## Read specular point locations
### Read SpOCK main input file to figure out stuff to then read the output
input_filename = 'lakes.txt'
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
gps_name = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
### Read SpOCK output files
nb_sc = 8 # !!!!!!!!!
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
lon_spec = []; lat_spec = []; gain_spec = []; gps_spec = []; date_spec = []
filename_spec_spock = []
nb_time_this_sc = []
print "Reading the SP positions..."
for isc in range(nb_sc):
    print isc, nb_sc-1, str(datetime.now())[0:19]
    which_sc = isc
    cyg = format(isc + 1, "02")
    filename_spec_spock.append( output_file_path_list[which_sc] + "specular_" + output_file_name_list[which_sc].replace(".txt",".bin") )

    data_spec = cygnss_read_spock_spec_bin(filename_spec_spock[-1], gps_name, dt, 0) 
    # date_spec.append(data_spec[0]); lon_spec.append(data_spec[1]); lat_spec.append(data_spec[2]); gain_spec.append(data_spec[3]); gps_spec.append(data_spec[4])
    
