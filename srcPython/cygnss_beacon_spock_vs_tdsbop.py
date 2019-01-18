# This script compares the specular point output of SpOCK
# with the specular point output of tds-bop
# for the planning of the CYGNSS beacon campaign
# Inputs:
# - input_filename: name of the SpOCK input filename
# - cygfm: which FM to look at
# ASSUMPTIONS:
# - before running this script, generate the ground station report:
# python report_coverage_ground_station_for_sift_parallel_sftp.py input_filename
import ipdb
import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import numpy as np

from cygnss_read_spock_spec_bin import *
import read_input_file; reload(read_input_file); from read_input_file import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta

input_filename = 'spock.txt'
cygfm = 3

# Read SpOCK main input file
var_in, var_in_order = read_input_file(input_filename)
dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, \
                                                            'dt_output')]; 
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order,\
                                                'output_file_name_list')];
gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order,\
                                                                'gps_name')];

# Read specular point output of SpOCK
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

# Select the times when the FM is in contact of the beacon station
## Read the ground station coverage report
gs_filename = output_file_path_list[isc] + 'coverage/' + 'report_all_by_' + \
              output_file_name_list[isc]
gs_file = open(gs_filename)
read_gs_file = gs_file.readlines()
ncontact = len(read_gs_file) - 2
date_coverage_spock = []
for icontact in range(ncontact):
    date_start = datetime.strptime(read_gs_file[2+icontact].split()[1] + 'T' +\
                                   read_gs_file[2+icontact].split()[2], "%Y-%m-%dT%H:%M:%S")
    date_stop = datetime.strptime(read_gs_file[2+icontact].split()[4] + 'T' +\
                                   read_gs_file[2+icontact].split()[5], "%Y-%m-%dT%H:%M:%S")
    date_coverage_spock.append([date_start, date_stop])
gs_file.close()

## For each contact, save the specular point information (positon, PRN, ...)
icontact = 0
ntime = len(date_spock)
date_spock_gs = []; lon_spock_gs = []; lat_spock_gs = []; gain_spock_gs = []
gps_spock_gs = []; x_cyg_spock_gs = []; y_cyg_spock_gs = []; z_cyg_spock_gs = []
x_gps_spock_gs = []; y_gps_spock_gs = []; z_gps_spock_gs = []
x_spec_spock_gs = []; y_spec_spock_gs = []; z_spec_spock_gs = []
el_spec_spock_gs = []; az_spec_spock_gs = []; el_spec_not_int_spock_gs = []
az_spec_not_int_spock_gs = []
itime = 0
while itime < ntime:
    if (date_spock[itime] >= date_coverage_spock[icontact][0]):
        date_spock_gs_sub = []; lon_spock_gs_sub = []; lat_spock_gs_sub = []; gain_spock_gs_sub = []
        gps_spock_gs_sub = []; x_cyg_spock_gs_sub = []; y_cyg_spock_gs_sub = []; z_cyg_spock_gs_sub = []
        x_gps_spock_gs_sub = []; y_gps_spock_gs_sub = []; z_gps_spock_gs_sub = []
        x_spec_spock_gs_sub = []; y_spec_spock_gs_sub = []; z_spec_spock_gs_sub = []
        el_spec_spock_gs_sub = []; az_spec_spock_gs_sub = []; el_spec_not_int_spock_gs_sub = []
        az_spec_not_int_spock_gs_sub = []
        while (date_spock[itime] <= date_coverage_spock[icontact][1]): # current date in
            # ground contact interval
            date_spock_gs_sub.append(date_spock[itime])
            lon_spock_gs_sub.append(lon_spock[itime])
            lat_spock_gs_sub.append(lat_spock[itime])
            gain_spock_gs_sub.append(gain_spock[itime])
            gps_spock_gs_sub.append(gps_spock[itime])
            x_cyg_spock_gs_sub.append(x_cyg_spock[itime])
            y_cyg_spock_gs_sub.append(y_cyg_spock[itime])
            z_cyg_spock_gs_sub.append(z_cyg_spock[itime])
            x_gps_spock_gs_sub.append(x_gps_spock[itime])
            y_gps_spock_gs_sub.append(y_gps_spock[itime])
            z_gps_spock_gs_sub.append(z_gps_spock[itime])
            x_spec_spock_gs_sub.append(x_spec_spock[itime])
            y_spec_spock_gs_sub.append(y_spec_spock[itime])
            z_spec_spock_gs_sub.append(z_spec_spock[itime])
            el_spec_spock_gs_sub.append(el_spec_spock[itime])
            az_spec_spock_gs_sub.append(az_spec_spock[itime])
            el_spec_not_int_spock_gs_sub.append(el_spec_not_int_spock[itime])
            az_spec_not_int_spock_gs_sub.append(az_spec_not_int_spock[itime])
            itime = itime + 1
        icontact = icontact + 1
        date_spock_gs.append(date_spock_gs_sub)
        lon_spock_gs.append(lon_spock_gs_sub)
        lat_spock_gs.append(lat_spock_gs_sub)
        gain_spock_gs.append(gain_spock_gs_sub)
        gps_spock_gs.append(gps_spock_gs_sub)
        x_cyg_spock_gs.append(x_cyg_spock_gs_sub)
        y_cyg_spock_gs.append(y_cyg_spock_gs_sub)
        z_cyg_spock_gs.append(z_cyg_spock_gs_sub)
        x_gps_spock_gs.append(x_gps_spock_gs_sub)
        y_gps_spock_gs.append(y_gps_spock_gs_sub)
        z_gps_spock_gs.append(z_gps_spock_gs_sub)
        x_spec_spock_gs.append(x_spec_spock_gs_sub)
        y_spec_spock_gs.append(y_spec_spock_gs_sub)
        z_spec_spock_gs.append(z_spec_spock_gs_sub)
        el_spec_spock_gs.append(el_spec_spock_gs_sub)
        az_spec_spock_gs.append(az_spec_spock_gs_sub)
        el_spec_not_int_spock_gs.append(el_spec_not_int_spock_gs_sub)
        az_spec_not_int_spock_gs.append(az_spec_not_int_spock_gs_sub)        
    else:        
        itime = itime + 1        
    
    if icontact >= ncontact:
        break
