# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

# cygnss_coverage_5years.py was written to answer Chris's email on Aug 19 2019

import sys
sys.path.append("/Users/cbv/work/spock/srcPython") 
import matplotlib
matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap, shiftgrid
from collections import *
from matplotlib import colors
import matplotlib.ticker as ticker
import ipdb
from read_input_file import *
from find_in_read_input_order_variables import *
from read_spec_spock import *
from cygnss_read_spock_spec_bin import *
import pickle


input_filename = sys.argv[1] #'spec_spock_2d_noIIF.txt'
var_in, var_in_order = read_input_file(input_filename)
dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
nsc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];
cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']


nspec = 4 
# Read SP positions
print 'Reading SP positions...'
for isc in range(nsc):
    print 'isc', isc
    spec_spock_filename = output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
    data_spec = cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'), gps_name_list_spock, dt_spock_output, 0) 
    date_spock = data_spec[0]; lon_spec_sc = data_spec[1]; lat_spec_sc = data_spec[2];   # !!!!!!! to take into account all gains: lon_spec_sc = data_spec[1]; lat_spec_sc = data_spec[2]; lon_spec_sc = data_spec[5]; lat_spec_sc = data_spec[6];  to take into account only SPs with non-0 gain
    gain_spec_sc = data_spec[3]; gps_name_sc = data_spec[4]; 
    if isc == 0:
        ntime = len(lon_spec_sc)
        lon_spec = np.zeros([ntime, nsc, nspec]) + 99999999*180/np.pi
        lat_spec = np.zeros([ntime, nsc, nspec]) + 99999999*180/np.pi
        gain_spec = np.zeros([ntime, nsc, nspec]) + 99999999*180/np.pi
    for itime in range(ntime):
        for ispec in range(len(lon_spec_sc[itime])):
            lon_spec[itime, isc, ispec] = lon_spec_sc[itime][ispec]
            lat_spec[itime, isc, ispec] = lat_spec_sc[itime][ispec]
            gain_spec[itime, isc, ispec] = gain_spec_sc[itime][ispec]
            

cell_width = 0.25#!!!!!!! should be 0.25 # in degrees
exclude_revisit_dt = cell_width*np.sqrt(2)*110/7.5 # if two visits are within exclude_revisit_dt seconds, then count the two visits as one visit
lat_max = 35
lat_min = -35
print_hour = 1. # print time every print_hour hour(s)
ncell_lon = (int)(360./cell_width) + 1
ncell_lat = (int)((lat_max - lat_min)/cell_width) + 1
ncell = ncell_lon*ncell_lat
cov_total = np.zeros([ncell_lon, ncell_lat])
time_visit = np.zeros([ncell_lon, ncell_lat, nsc*nspec]) -1 
date_ref = date_spock[0]
cov_time = np.zeros([ntime])
cell_covered = np.zeros([ncell_lon, ncell_lat])
ncell_covered_time = np.zeros([ntime])
itime = 0
#stop_recording_visit = 0
while ((itime < ntime)):
    for isc in range(nsc):
        for ispec in range(nspec):
            icell_lon = (int)(lon_spec[itime, isc, ispec] / cell_width)
            if ((lat_spec[itime, isc, ispec] >= lat_min) & (lat_spec[itime, isc, ispec] <= lat_max)): # theya re SPs at latitudes > 35 deg
                if gain_spec[itime, isc, ispec] >= 4: #!!!!!!! sure you want to keep the condition "if gain_spec[itime, isc, ispec]	> 0"
                    icell_lat = (int)((lat_spec[itime, isc, ispec] - lat_min) / cell_width)
                    #if stop_recording_visit == 0:
                    time_visit[icell_lon, icell_lat, isc*nspec + ispec] = (date_spock[itime] - date_ref).total_seconds()
                    cell_covered[icell_lon, icell_lat] = cell_covered[icell_lon, icell_lat] + 1
    ncell_covered_time[itime] = len(np.where(cell_covered != 0)[0])
    cov_time[itime] = ncell_covered_time[itime] * 100./ncell
    #if cov_time[itime]  >= 70.:
    if itime >= 30*3600 + 1:
        #stop_recording_visit = 1
        break
    if np.mod(itime, print_hour*3600) == 0:
        print itime, ntime, str(itime/3600) + 'h', format(cov_time[itime], ".0f") + '%'    
    itime = itime + 1
revisit_time = []

for icell_lon in range(ncell_lon):
    #print icell_lon, ncell_lon
    for icell_lat in range(ncell_lat):
        visiter = np.where(time_visit[icell_lon, icell_lat, :] != -1)[0]
        if len(visiter) > 0:
            time_visit_temp_not_sorted = time_visit[icell_lon, icell_lat, visiter]
            sort_time_visit = np.argsort(time_visit_temp_not_sorted)
            time_visit_temp = time_visit_temp_not_sorted[sort_time_visit]
            nvisit = len(time_visit_temp)
            time_visit_valid = []
            time_visit_valid.append(time_visit_temp[0])
            for ivisit in range(nvisit-1):
                if time_visit_temp[ivisit+1] - time_visit_temp[ivisit] > exclude_revisit_dt: #  if two visits are within exclude_revisit_dt seconds, then count the two visits as one visit 
                    time_visit_valid.append(time_visit_temp[ivisit+1])
                    revisit_time.append( time_visit_valid[-1] - time_visit_valid[-2] )
            if len(time_visit_valid) > 0: # this cell has been visited at least once
                cov_total[icell_lon, icell_lat] = cov_total[icell_lon, icell_lat] + 1
revisit_time = np.array(revisit_time)
ilon_cell_not_covered = np.where(cov_total == 0)[0]
ilat_cell_not_covered = np.where(cov_total == 0)[0]
ncell_not_covered = len(ilon_cell_not_covered)
ncell_covered = ncell - ncell_not_covered
perc_cov = ncell_covered * 100./ncell
mean_revisit_time = np.mean(revisit_time) / 3600. # in hours

print input_filename
print 'Percentage coverage: ' + format(perc_cov, ".1f")
print 'Mean revisit time: ' + format(mean_revisit_time, ".1f")
