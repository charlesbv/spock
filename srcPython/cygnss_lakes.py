# This script answers Chris's questions in his email on 08/12/2019

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
#input_filename = '082819.txt'
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

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

def cygnss_lakes(input_filename):
    ## Read specular point locations
    ### Read SpOCK main input file to figure out stuff to then read the output
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
    filename_spec_spock = []
    nb_time_this_sc = []
    min_lat_lon = [[32.992112, -107.335678], [13.900429, -89.585032], [13.713906, -89.107452], [37.157854, 99.701568]]
    max_lat_lon = [[32.890515, -107.241349], [13.831710, -89.515894], [13.629988, -88.986665], [36.548034, 100.745112]]
    name = ['Caballo', 'Coatepeque', 'LLopango', 'Qinghai']
    min_gain = 3

    nlake = len(name)
    min_lat_lon0_360 = []
    max_lat_lon0_360 = []
    for ilake in range(nlake):
        min_lat = min_lat_lon[ilake][0]
        min_lon = min_lat_lon[ilake][1]
        min_lon_corr = min_lon
        if min_lon < 0:
            min_lon_corr = 360 + min_lon
        min_lat_lon0_360.append([min_lat, min_lon_corr])
        max_lat = max_lat_lon[ilake][0]
        max_lon = max_lat_lon[ilake][1]
        max_lon_corr = max_lon
        if max_lon < 0:
            max_lon_corr = 360 + max_lon
        max_lat_lon0_360.append([max_lat, max_lon_corr])

    lon_spec = []; lat_spec = []; gain_spec = []
    visit_time = []
    visit_which_sp = [] # whcih sp
    visit_which_sc = [] # whcih sp
    for ilake in range(nlake):
        min_lat = min_lat_lon0_360[ilake][0]
        min_lon = min_lat_lon0_360[ilake][1]
        max_lat = max_lat_lon0_360[ilake][0]
        max_lon = max_lat_lon0_360[ilake][1]
        #print min_lat, min_lon, max_lat, max_lon
        if ilake == 0: # only read the sp position once (for the first lake)
            print "Reading the SP positions..."
            for isc in range(nb_sc):
                #print isc, nb_sc-1, str(datetime.now())[0:19]
                which_sc = isc
                cyg = format(isc + 1, "02")
                filename_spec_spock.append( output_file_path_list[which_sc] + "specular_" + output_file_name_list[which_sc].replace(".txt",".bin") )
                data_spec = cygnss_read_spock_spec_bin(filename_spec_spock[-1], gps_name, dt, 0)
                if isc == 0:
                    date_spec = np.array(data_spec[0])
                lon_spec.append(np.array(data_spec[1])); lat_spec.append(np.array(data_spec[2])); gain_spec.append(np.array(data_spec[3]));
        visit_time_lake = []
        visit_time_lake_temp = []
        visit_which_lake_temp = []
        visit_which_lake = []
        visit_sc_lake_temp = []
        visit_sc_lake = []

        for isc in range(nb_sc):
            index_sp_in = np.where((lon_spec[isc] >= min_lon) & (lon_spec[isc] <= max_lon) & (lat_spec[isc] >= min_lat) & (lat_spec[isc] <= max_lat) & (gain_spec[isc] > min_gain))
            time_sp_in = index_sp_in[0]
            which_sp_in = index_sp_in[1]
            date_spec_in = date_spec[time_sp_in]
            lon_spec_in = lon_spec[isc][time_sp_in, which_sp_in]
            lat_spec_in = lat_spec[isc][time_sp_in, which_sp_in]
            gain_spec_in = gain_spec[isc][time_sp_in, which_sp_in]
            visit_time_lake_temp = visit_time_lake_temp + list(time_sp_in) # visit_time_lake_temp is a concatenation of the visit times of all sc
            visit_which_lake_temp = visit_which_lake_temp + list(which_sp_in)
            visit_sc_lake_temp = visit_sc_lake_temp + list(np.zeros([len(which_sp_in)]) + isc) 

        visit_time_lake_temp_arr = np.array(visit_time_lake_temp)
        index_sort_visit_time_lake_temp_arr = np.argsort(visit_time_lake_temp_arr)
        visit_time_lake_temp_arr_sort = visit_time_lake_temp_arr[index_sort_visit_time_lake_temp_arr]
        visit_which_lake_temp_arr_sort = np.array(visit_which_lake_temp)[index_sort_visit_time_lake_temp_arr]
        visit_sc_lake_temp_arr_sort = np.array(visit_sc_lake_temp)[index_sort_visit_time_lake_temp_arr]    
        nvisit = len(visit_time_lake_temp_arr_sort)
        if nvisit > 0:
            visit_time_lake.append(visit_time_lake_temp_arr_sort[0])
            visit_which_lake.append(visit_which_lake_temp_arr_sort[0])
            visit_sc_lake.append(visit_sc_lake_temp_arr_sort[0])    
            max_dim_lake = np.max([max_lon - min_lon, max_lat - min_lat])
            min_time_between_two_visit = max_dim_lake*np.sqrt(2)*110/7.5
            for ivisit in range(1, nvisit):
                if visit_time_lake_temp_arr_sort[ivisit] - visit_time_lake[-1] > min_time_between_two_visit:
                    visit_time_lake.append(visit_time_lake_temp_arr_sort[ivisit])
                    visit_which_lake.append(visit_which_lake_temp_arr_sort[ivisit])
                    visit_sc_lake.append(visit_sc_lake_temp_arr_sort[ivisit])
        visit_time.append(visit_time_lake)
        visit_which_sp.append(visit_which_lake)
        visit_which_sc.append(visit_sc_lake)

    return  visit_time, visit_which_sp, visit_which_sc



