
# stopped at iday = 1673 (2004-10-24.txt crashed)
# This script downloads all TLEs for satellites de-1, viking, polar, and image during their lifetime an plot the argument of apogee (as read from the TLEs) as a function of tim
# ASSUMPTIONS:
# - in the SpOCK simulation, station names must include the character "NP" (if at the North Pole) or "SP" (if at the South Pole)
# - for the post processing of SPOCK in looking at the coverage of both the South and North, it only works if there are two sc (not 1, not 3, not 4 etc). The initial idea of this analysis is to look at one 1 sc sees the stations at the North AND the other sc sees the stations at the South. So it makes sense that we look at 2 sc and only 2 sc
# - see section "PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT"

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
download_tle = 0
polar_info = ['Polar', '23802', '1996-02-24', '2008-04-28', 'blue'] # important to put the dates when the sc was operational (and not just in space after being decomissionned)
image_info = ['IMAGE', '26113', '2000-03-25', '2005-12-18', 'red']
de1_info = ['DE-1', '12624', '1981-08-03', '1991-02-28', 'black']
viking_info = ['Viking', '16614', '1986-02-22', '1987-05-12', 'magenta']
arr_info = [polar_info, image_info, de1_info, viking_info]
pole_offset = 30.
#min_nb_stations = 6
fov = [15,8]#[7.5, 7.5, 7.5, 7.5]#[15, 8] # first umber for Polar, second for IMAGE !!!!!! important to put Polar first then IMAGE in the SpOCK simuation (#ORBIT section)
min_min_nb_stations = 6
max_min_nb_stations = 12
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

xaxis_start_date = '1981-08-03' # need ot have the day to be 01
xaxis_stop_date = '2008-05-01' # need ot have the day to be 01
xaxis_nticks = 8


import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
from spock_main_input_sgp4 import *
import os
from datetime import datetime, timedelta
import ipdb
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from convert_tle_date_to_date import *
from matplotlib.lines import Line2D
import pickle
from read_input_file import *
from find_in_read_input_order_variables import *
from report_coverage_ground_station_amarex import *
from read_output_file import *

def cov_both_pole(spock_input_filename): # !!!! before calling this function, need to call the function report_coverage_ground_station_amarex
    var_in, var_in_order = read_input_file(spock_input_filename)
    date_start_spock = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')];
    date_stop_spock = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')];
    dt_spock = var_in[find_in_read_input_order_variables(var_in_order, 'dt')];  # the coverage is output every dt, not every dt_output
    nb_steps_spock = (int)((date_stop_spock - date_start_spock).total_seconds())#/dt_spock) # need to calcualte it here because read_input_file gives the nb of time steps output considering dt_output, not dt. Here we're interested in dt
    cov_spock = pickle.load(open('pickle/' +spock_input_filename.replace('.txt', '.pickle')))
    nsc_spock = len(cov_spock)
    north_stations_spock = np.zeros([nsc_spock, nb_steps_spock])
    south_stations_spock = np.zeros([nsc_spock, nb_steps_spock])
    north = []; south = []
    for isc in range(nsc_spock): #[0,2]:#!!!!!!range(nsc_spock):
        ncov = len(cov_spock[isc])
        for icov in range(ncov):
            cov_start = (int)((datetime.strptime(cov_spock[isc][icov][0], "%Y-%m-%dT%H:%M:%S") - date_start_spock).total_seconds())
            cov_stop = (int)((datetime.strptime(cov_spock[isc][icov][1], "%Y-%m-%dT%H:%M:%S") - date_start_spock).total_seconds())
            min_aov = np.float(cov_spock[isc][icov][3])

            if 'NP' in cov_spock[isc][icov][2]: # North Pole station
                north_stations_spock[isc, cov_start: cov_stop+1] = north_stations_spock[isc, cov_start: cov_stop+1] + 1
            if 'SP' in cov_spock[isc][icov][2]: # South Pole station
                south_stations_spock[isc, cov_start: cov_stop+1] = south_stations_spock[isc, cov_start: cov_stop+1] + 1


        second_list_north = (np.where( (north_stations_spock[isc, :] >= min_nb_stations))[0]).tolist()
        north = north + list(set(second_list_north) - set(north)) # concatenate north and list( np.where( (north_stations_spock[isc, :] >= min_nb_stations)[0] )) removing duplicates
        second_list_south = (np.where( (south_stations_spock[isc, :] >= min_nb_stations))[0]).tolist()
        south = south + list(set(second_list_south) - set(south)) 
    north.sort()
    south.sort()

    north_and_south = list(set(north) & set(south))
    north_and_south.sort()
    north_south_two_sc = north_and_south
    # Below works only for two sc (I checked that it gives the same result at north_and_south)
    # north_south_two_sc = np.where( ((north_stations_spock[0, :] >= min_nb_stations) & # sc 1 sees north stations AND
    #                                 (south_stations_spock[1, :] >= min_nb_stations)) # sc 2 sees south stations
    #                                                                                  # OR 
    #                                | ((south_stations_spock[0, :] >= min_nb_stations) & # sc 1 sees south stations AND
    #                                   (north_stations_spock[1, :] >= min_nb_stations)) ) # sc 2 sees north stations
    # north_south_two_sc = north_south_two_sc[0]
    # end of Below works only for two sc

    nb_seconds_both_pole = np.float(len(north_south_two_sc))
    nb_hours_both_pole = nb_seconds_both_pole / 3600.

    return north_south_two_sc 


def cov_both_pole_latitude_limitation(spock_input_filename): # Same as cov_both_pole but only record times when latitude is above lat_min (defined in this function) (see Aaron's email on 08/28/2019)!!!! before calling this function, need to call the function report_coverage_ground_station_amarex
    var_in, var_in_order = read_input_file(spock_input_filename)
    date_start_spock = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')];
    date_stop_spock = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')];
    dt_spock = var_in[find_in_read_input_order_variables(var_in_order, 'dt')];  # the coverage is output every dt, not every dt_output
    dt_output_spock = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 

    nb_steps_spock = (int)((date_stop_spock - date_start_spock).total_seconds())#/dt_spock) # need to calcualte it here because read_input_file gives the nb of time steps output considering dt_output, not dt. Here we're interested in dt
    cov_spock = pickle.load(open('pickle/' +spock_input_filename.replace('.txt', '.pickle')))
    var_to_read = ['latitude']
    nsc_spock = len(cov_spock)
    north_stations_spock = np.zeros([nsc_spock, nb_steps_spock])
    south_stations_spock = np.zeros([nsc_spock, nb_steps_spock])
    north = []; south = []
    for isc in [0, 2]: #!!!!!range(nsc_spock):
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        lat = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
        nb_seconds_spock_amarex = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
        lat_min = 35. # only record times when lat is greater than lat_min (abosolute values)
        index_lat_ok = (np.where(np.abs(lat) >= lat_min)[0] * dt_output_spock).astype(int).tolist() # multiply by dt_output because we want these indices to represent seconds (not minutes if for instance dt_output = 60)
        index_lat_ok_temp_temp = (np.where(np.abs(lat) >= lat_min)[0]).astype(int)
        index_lat_ok_temp = nb_seconds_spock_amarex[index_lat_ok_temp_temp].astype(int) # multiply by dt_output because we want these indices to represent seconds (not minutes if for instance dt_output = 60)
        nlat_ok = len(index_lat_ok_temp)
        index_lat_ok = [] # index_lat_ok fills the gap between two cnosecutive time steps dt_output (ie it add all seconds between two dt_output time steps)
        for ilat in range(1, nlat_ok):
            if (index_lat_ok_temp[ilat] - index_lat_ok_temp[ilat-1]) <= dt_output_spock:
                index_lat_ok = index_lat_ok + range(index_lat_ok_temp[ilat-1], index_lat_ok_temp[ilat])
        ncov = len(cov_spock[isc])
        for icov in range(ncov):
            cov_start = (int)((datetime.strptime(cov_spock[isc][icov][0], "%Y-%m-%dT%H:%M:%S") - date_start_spock).total_seconds())
            cov_stop = (int)((datetime.strptime(cov_spock[isc][icov][1], "%Y-%m-%dT%H:%M:%S") - date_start_spock).total_seconds())
            min_aov = np.float(cov_spock[isc][icov][3])

            if 'NP' in cov_spock[isc][icov][2]: # North Pole station
                north_stations_spock[isc, cov_start: cov_stop+1] = north_stations_spock[isc, cov_start: cov_stop+1] + 1
            if 'SP' in cov_spock[isc][icov][2]: # South Pole station
                south_stations_spock[isc, cov_start: cov_stop+1] = south_stations_spock[isc, cov_start: cov_stop+1] + 1


        second_list_north = (np.where( (north_stations_spock[isc, :] >= min_nb_stations))[0]).tolist()
        north = north + list(set(second_list_north) - set(north)) # concatenate north and list( np.where( (north_stations_spock[isc, :] >= min_nb_stations)[0] )) removing duplicates
        second_list_south = (np.where( (south_stations_spock[isc, :] >= min_nb_stations))[0]).tolist()
        south_lat_limit_temp = list(set(second_list_south) & set(index_lat_ok)) # select times when the sc lat is > lat_min (doesn't matter if the interesection of index_lat_ok is with south or north
        south = south + list(set(south_lat_limit_temp) - set(south))
    north.sort()
    south.sort()

    north_and_south = list(set(north) & set(south))
    north_and_south.sort()
    north_south_two_sc_lat_limit = north_and_south
    nb_seconds_both_pole = np.float(len(north_south_two_sc_lat_limit))
    nb_hours_both_pole = nb_seconds_both_pole / 3600.

    return north_south_two_sc_lat_limit 


def cov_pole_one_sc(spock_input_filename, isc): # !!!! before calling this function, need to call the function report_coverage_ground_station_amarex. isc is the sc # in Section #ORBBIT of spock_input_filename that you're interested in
    var_in, var_in_order = read_input_file(spock_input_filename)
    date_start_spock = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')];
    date_stop_spock = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')];
    dt_spock = var_in[find_in_read_input_order_variables(var_in_order, 'dt')];  # the coverage is output every dt, not every dt_output
    nb_steps_spock = (int)((date_stop_spock - date_start_spock).total_seconds())#/dt_spock) # need to calcualte it here because read_input_file gives the nb of time steps output considering dt_output, not dt. Here we're interested in dt
    cov_spock = pickle.load(open('pickle/' +spock_input_filename.replace('.txt', '.pickle')))
    nsc_spock = len(cov_spock)
    north_stations_spock = np.zeros([nsc_spock, nb_steps_spock])
    south_stations_spock = np.zeros([nsc_spock, nb_steps_spock])

    ncov = len(cov_spock[isc])
    for icov in range(ncov):
        cov_start = (int)((datetime.strptime(cov_spock[isc][icov][0], "%Y-%m-%dT%H:%M:%S") - date_start_spock).total_seconds())
        cov_stop = (int)((datetime.strptime(cov_spock[isc][icov][1], "%Y-%m-%dT%H:%M:%S") - date_start_spock).total_seconds())
        min_aov = np.float(cov_spock[isc][icov][3])

        if 'NP' in cov_spock[isc][icov][2]: # North Pole station
            north_stations_spock[isc, cov_start: cov_stop+1] = north_stations_spock[isc, cov_start: cov_stop+1] + 1
        if 'SP' in cov_spock[isc][icov][2]: # South Pole station
            south_stations_spock[isc, cov_start: cov_stop+1] = south_stations_spock[isc, cov_start: cov_stop+1] + 1
    north_south = np.where( (north_stations_spock[isc, :] >= min_nb_stations) | # sc sees north stations OR
                                   (south_stations_spock[isc, :] >= min_nb_stations)) # sc sees south stations AND
    north_south = north_south[0]
    nb_seconds_pole = np.float(len(north_south))
    nb_hours_pole = nb_seconds_pole / 3600.

    return north_south

#bla = cov_both_pole('032800.txt')
fov = np.array(fov)

# # amarex
# amarex_filename = 'amarex4.txt'
# #report_coverage_ground_station_amarex(amarex_filename, fov)
# print 'Done writing the report of ground stations'
# nb_events_amarex = np.zeros([max_min_nb_stations - min_min_nb_stations+1, 13])
# var_in, var_in_order = read_input_file(amarex_filename)
# date_start_spock_amarex = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')];
# date_stop_spock_amarex = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')];
# nday_amarex = (date_stop_spock_amarex - date_start_spock_amarex).days
# cov_spock_both_pole_amarex = np.zeros([nday_amarex])
# imin_station = 0
# for min_nb_stations in range(min_min_nb_stations, max_min_nb_stations+1):
#     print 'min_nb_stations', min_nb_stations
#     cov_spock_both_pole_amarex_temp = np.array(cov_both_pole(amarex_filename)) # !!!!!!!! looks like the latitude limitation doens't make a difference with 4 sc. In theory here it should be np.array(cov_both_pole_latitude_limitation(amarex_filename)) 
#     # While cov_both_pole is used for a simulation of one day, cov_both_pole_latitude_limitation not necesarily. Howver, we want daily statistics so we need to subdivide cov_spock_both_pole_amarex into daily sub arrays
#     index_start = 0
#     for iday in range(nday_amarex):
#         index_stop_temp = np.where(cov_spock_both_pole_amarex_temp <= 86400*(iday+1))[0] + 1
#         if len(index_stop_temp) > 0:
#             index_stop = index_stop_temp[-1]
#             cov_spock_both_pole_amarex_daily = cov_spock_both_pole_amarex_temp[index_start:index_stop]
#             index_start = index_stop
#             cov_spock_both_pole_amarex[iday] = np.float(len(cov_spock_both_pole_amarex_daily))
#     print cov_spock_both_pole_amarex/3600.
#     print ('      ' + str((int)(min_nb_stations)) + '       ')
#     for ihour_obs in range(0,(int)(np.max(cov_spock_both_pole_amarex)/3600)+1):
#         nb_events_amarex[imin_station, ihour_obs] = len(np.where(cov_spock_both_pole_amarex >= ihour_obs*3600)[0])
#         print  (str(nb_events_amarex[imin_station, ihour_obs] ).replace('0','X').zfill(6).replace('0',' ').replace('X','0')),
#     imin_station = imin_station + 1
# raise Exception
# all below is for polar, image and the other non-amarex sc
nsc = len(arr_info)
arg_apogee = []
arg_apogee_trans = []
epoch_start = []
epoch_stop = []
tle_epoch = []
handles_arr = []
for isc in range(nsc):
    epoch_start.append(datetime.strptime(arr_info[isc][2], "%Y-%m-%d"))
    epoch_stop.append(datetime.strptime(arr_info[isc][3], "%Y-%m-%d"))
    handles_arr.append(Line2D([0], [0], color =  arr_info[isc][-1], lw=4, label = arr_info[isc][0]))

xaxis_start_date_date = datetime.strptime(xaxis_start_date, "%Y-%m-%d")
date_ref = xaxis_start_date_date#np.min(epoch_start)
nb_seconds_since_ref = []
read_tle_file = []
for isc in range(nsc):
    if download_tle == 1:
        os.system('python download_tle.py ' + arr_info[isc][1] + ' ' + arr_info[isc][2] + ' ' + arr_info[isc][3] )
    tle_filename = arr_info[isc][1] + '_' + arr_info[isc][2] + '_' + arr_info[isc][3] + '.txt'
    tle_file = open(tle_filename)
    read_tle_file.append(tle_file.readlines())
    ntle = len(read_tle_file[-1]) / 2
    arg_apogee_sc = []
    arg_apogee_trans_sc = []
    tle_epoch_sc = []
    nb_seconds_since_ref_sc = []
    for itle in range(ntle):
        iline = itle * 2
        arg_apogee_temp = np.mod(np.float(read_tle_file[-1][iline + 1][34:42]) + 180., 360)
        arg_apogee_sc.append(arg_apogee_temp)
        if ((arg_apogee_temp >= 0) & (arg_apogee_temp < 90)):
            arg_apogee_trans_temp = arg_apogee_temp
        elif ((arg_apogee_temp >= 90) & (arg_apogee_temp < 270)):
            arg_apogee_trans_temp = 180 - arg_apogee_temp
        elif ((arg_apogee_temp >= 270) & (arg_apogee_temp < 360)):
            arg_apogee_trans_temp = arg_apogee_temp - 360
        arg_apogee_trans_sc.append(arg_apogee_trans_temp)
        tle_epoch_sc_raw = read_tle_file[-1][iline][18:32]
        tle_epoch_sc_temp = convert_tle_date_to_date(tle_epoch_sc_raw)
        tle_epoch_sc.append(tle_epoch_sc_temp)
        nb_seconds_since_ref_sc.append((tle_epoch_sc_temp - date_ref).total_seconds())
    arg_apogee.append(arg_apogee_sc)
    arg_apogee_trans.append(arg_apogee_trans_sc)
    tle_epoch.append(np.array(tle_epoch_sc))
    nb_seconds_since_ref.append(nb_seconds_since_ref_sc)


# FIGURING OUT WITH SPOCK WHEN POLAR AND IMAGE SEE THE NORTH AND SOUTH POLES
isc_polar = 0 # !!!!!! must be the sc # in Section #ORBIT of spock_input_filename and be the sc # in arr_info
isc_image = 1 # !!!!!! must be the sc # in Section #ORBIT of spock_input_filename and be the sc # in arr_info
tle_start = tle_epoch[isc_image][0]# only select the times when the two sc were operational at the same time. In this case, these times coorespond to the times when IMAGE was operational
tle_stop = tle_epoch[isc_image][-1]
# for each day between tle_start and tle_stop, propagate Polar and IMAGE for this entire day, using the most recent TLEs for these two sc
date_start_spock_ana = datetime.strftime(tle_start + timedelta(days = 1), "%Y-%m-%d")[0:10] + 'T00:00:00'
date_stop_spock_ana = datetime.strftime(tle_stop, "%Y-%m-%d")[0:10] + 'T00:00:00'
date_start_spock_ana_date = datetime.strptime(date_start_spock_ana, "%Y-%m-%dT%H:%M:%S")
date_stop_spock_ana_date = datetime.strptime(date_stop_spock_ana, "%Y-%m-%dT%H:%M:%S")
nb_day_spock_ana = (date_stop_spock_ana_date - date_start_spock_ana_date).days + 1
date_spock_ana = np.array([date_start_spock_ana_date + timedelta(days=i) for i in np.arange(0, nb_day_spock_ana, 1)])
nb_day_spock_ana = len(date_spock_ana)
imin_station = 0
nb_events = np.zeros([max_min_nb_stations - min_min_nb_stations + 1, 13])
date_events = ['2000-08-04 01:15',
'2001-05-28 03:45',
'2001-08-13 22:41',
'2002-11-05 11:21',
'2001-07-02 04:51',
'2002-10-23 11:30',
'2001-08-03 23:40',
'2001-11-15 17:40',
'2002-11-03 04:26'
]
nevent_interest  = len(date_events)

for min_nb_stations in range(min_min_nb_stations, max_min_nb_stations+1):
    cov_nb_seconds_since_start = []
    cov_spock_both_pole = []
    cov_spock_pole_one_sc = []
    ov_spock_both_pole = []
    for iday in range(nb_day_spock_ana):
        date_start_spock_date = date_spock_ana[iday]
        date_stop_spock_date = date_spock_ana[iday] + timedelta(days = 1)
        date_start_spock = datetime.strftime(date_start_spock_date, "%Y-%m-%dT%H:%M:%S")
        date_stop_spock = datetime.strftime(date_stop_spock_date, "%Y-%m-%dT%H:%M:%S")
        cov_nb_seconds_since_start.append((date_start_spock_date - date_spock_ana[0]).total_seconds())
        spock_input_filename = date_start_spock[0:10] + '.txt'
    #     # Create TLE file for Polar and IMAGE ##!!!!!!! put Polar first, then IMAGE. Otherwise need to change the fov order
    #     ## Polar
    #     itle_polar = np.where(tle_epoch[isc_polar] < date_start_spock_date)[0][-1]
    #     tle_filename_polar = date_start_spock[0:10] + '_' + arr_info[isc_polar][0].lower() + '.txt'
    #     tle_file_polar = open(tle_filename_polar, "w")
    #     print >> tle_file_polar, read_tle_file[isc_polar][itle_polar*2].replace('\r','').replace('\n','')
    #     print >> tle_file_polar, read_tle_file[isc_polar][itle_polar*2+1].replace('\r','').replace('\n','')
    #     tle_file_polar.close()
    #     ## Image
    #     itle_image = np.where(tle_epoch[isc_image] < date_start_spock_date)[0][-1]
    #     tle_filename_image = date_start_spock[0:10] + '_' + arr_info[isc_image][0].lower() + '.txt'
    #     tle_file_image = open(tle_filename_image, "w")
    #     print >> tle_file_image, read_tle_file[isc_image][itle_image*2].replace('\r','').replace('\n','')
    #     print >> tle_file_image, read_tle_file[isc_image][itle_image*2+1].replace('\r','').replace('\n','')
    #     tle_file_image.close()
    #     # Create SpOCK main input file
    #     spock_main_input_sgp4(
    #      spock_input_filename,
    #     # for TIME section
    #     date_start_spock,
    #     date_stop_spock,
    #     60.,
    #     # for SPACECRAFT section
    #     2,
    #         '0',
    #     29,
    #         "cygnss_geometry_2016_acco09.txt",
    #     # for ORBIT section
    #     tle_filename_polar + '\n' + tle_filename_image,
    #     # for FORCES section
    #     8,
    #     'drag sun_gravity moon_gravity',
    #     'static',
    #     # for OUTPUT section
    #     'out/out',
    #     60.,
    #     # for ATTITUDE section
    #     "nadir",#"(0;-22;0) (0;0;0)",#"nadir", #!!!!!
    #     # for GROUNDS_STATIONS section
    #     "stations.txt",#"my_ground_stations.txt"
    #      # for SPICE section
    #      '/Users/cbv/cspice/data',
    #      # for DENSITY_MOD section
    # 1
    #     )
    #     # if iday > -1:
    #     # # Run SpOCK
    #     #     #os.system('mpirun -np 2 spock_amarex ' + spock_input_filename)
    #     #     #Run report_coverage_ground_station_amarex to create the pickle of the statistic coverge
    #     #     report_coverage_ground_station_amarex(spock_input_filename, fov) # !!!!!! important to put Polar first then IMAGE in th
        # e SpOCK simuation (#ORBIT section)
        #Calculate the number of seconds during whcih Polar sees the North pole and IMAGE sees the South pole or vice versa
        cov_spock_both_pole.append( np.float(len(cov_both_pole(spock_input_filename))) )
        cov_spock_pole_one_sc.append( np.float(len(cov_pole_one_sc(spock_input_filename, isc_image))) )
        print iday, nb_day_spock_ana, str(date_start_spock_date)[:10], cov_spock_both_pole[-1] #cov_spock_pole_one_sc[-1]#cov_spock_both_pole[-1]
        # if iday == 0:
        #     raise Exception
        
    cov_spock_both_pole_arr = np.array(cov_spock_both_pole)
    cov_spock_pole_one_sc_arr = np.array(cov_spock_pole_one_sc)

    print ('      ' + str((int)(min_nb_stations)) + '       ')
    for ihour_obs in range(0,(int)(np.max(cov_spock_both_pole_arr)/3600) + 1):
        nb_events[imin_station, ihour_obs] = len(np.where(cov_spock_both_pole_arr >= ihour_obs*3600)[0])
        print  (str(len(np.where(cov_spock_both_pole_arr >= ihour_obs*3600)[0])).replace('0','X').zfill(6).replace('0',' ').replace('X','0')),
    imin_station = imin_station + 1
    # for ievent in range(nevent_interest):
    # #ievent = 0
    #     spock_input_filename = date_events[ievent][:10]+ '.txt'#'2000-08-04.txt'
    #     north_south_two_sc_out = cov_both_pole(spock_input_filename)
    #     #################
    #     event_time = datetime.strptime(date_events[ievent], "%Y-%m-%d %H:%M")
    #     date_start_spock_date = datetime.strptime(date_events[ievent][:10], "%Y-%m-%d")
    #     nb_steps_spock = 86400
    #     date_cov_arr =  np.array([date_start_spock_date + timedelta(seconds=i) for i in np.arange(0, nb_steps_spock, 1)])
    #     print date_events[ievent], (event_time in date_cov_arr[north_south_two_sc_out])

raise Exception

xaxis_start_plot_str = '2000-03-01' # need ot have the day to be 01
xaxis_stop_plot_str = '2006-01-01' # need ot have the day to be 01
xaxis_nticks_one_sc = 8

xaxis_start_plot = datetime.strptime(xaxis_start_plot_str, "%Y-%m-%d")

height_fig = 15.  # the width is calculated as height_fig * 4/3. 
fontsize_plot = 25
ratio_fig_size = 4./3
fig_title = 'Percentage of hours per day that IMAGE sees at least ' + str((int)(min_nb_stations)) + ' stations of a pole'
y_label = '% hours/day'
x_label = 'Real time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                  
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
ax = fig.add_subplot(gs[0, 0])
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           

# len(cov_spock_pole_one_sc_arr) is the number of days -> xaxis is in seconds
start_xaxis_one_sc = (date_spock_ana[0] - xaxis_start_plot).total_seconds()
xaxis_one_sc = np.arange(start_xaxis_one_sc, start_xaxis_one_sc + len(cov_spock_pole_one_sc_arr) * 24*3600, 24*3600)
ax.plot(xaxis_one_sc, cov_spock_pole_one_sc_arr/(24*3600.) * 100, linewidth = 2, color = 'black')

# x axis label is in real time
xaxis_stop_plot = datetime.strptime(xaxis_stop_plot_str, "%Y-%m-%d")
delta_xaxis_start_stop = (xaxis_stop_plot - xaxis_start_plot).total_seconds() / 3600./ 24 / 30
nb_month_per_tick = (int)(delta_xaxis_start_stop / xaxis_nticks_one_sc) + 1


xticks = np.zeros([xaxis_nticks_one_sc+1])
date_list_str = []
xprevious = xaxis_start_plot
date_list_str.append(datetime.strftime(xaxis_start_plot, "%b") + '\n' + str(xaxis_start_plot)[:4])
print nb_month_per_tick
for itick in range(xaxis_nticks_one_sc-1):
    #ipdb.set_trace()
    xtick_temp = datetime.strftime(xprevious, "%Y-%m")
    xtick_temp_year = (int)(xtick_temp[:4])
    xtick_temp_month = (int)(xtick_temp[5:7])
    if xtick_temp_month + nb_month_per_tick > 12:
        xtick_temp_year_new = xtick_temp_year + (int)((xtick_temp_month + nb_month_per_tick)/12.)
        xtick_temp_month_new = np.mod(xtick_temp_month + nb_month_per_tick, 12)
        if xtick_temp_month_new == 0: # month 24, 36, 48, etc is Dec and the year should not be incremented by 2, 3, 4, etc but only 1, 2, 3, etc
            xtick_temp_month_new = 12
            xtick_temp_year_new = xtick_temp_year_new - 1
    else:
        xtick_temp_month_new = xtick_temp_month + nb_month_per_tick
        xtick_temp_year_new = xtick_temp_year
    xtick_temp_str = str(xtick_temp_year_new) + '-' + str(xtick_temp_month_new).zfill(2) + '-01'
    xtick_temp_date = datetime.strptime(xtick_temp_str, "%Y-%m-%d")
    xticks[itick+1] = (xtick_temp_date - xaxis_start_plot).total_seconds() # xticks[0] is xaxis_start_plot
    date_list_str.append(datetime.strftime(xtick_temp_date, "%b") + '\n' + xtick_temp_str[:4])
    print str(xprevious)[0:10], str(xtick_temp_date)[0:10]
    xprevious = xtick_temp_date
xticks[-1] = (xaxis_stop_plot  - xaxis_start_plot).total_seconds()
date_list_str.append(datetime.strftime(xaxis_stop_plot, "%b") + '\n' + str(xaxis_stop_plot)[:4])
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])


fig.set_figheight(height_fig)
fig.set_figwidth(height_fig*ratio_fig_size)
fig_save_name = 'cov_image_' + str((int)(min_nb_stations)) + '_stations_day.pdf'

fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')

raise Exception
height_fig = 15.  # the width is calculated as height_fig * 4/3. 
fontsize_plot = 25
ratio_fig_size = 4./3
fig_title = 'Histogram of # hours/day during which Polar and IMAGE see both poles together (one pole each)'
y_label = 'Percentage (%)'
x_label = '% hours'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.958,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                                                                                                                  
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
ax = fig.add_subplot(gs[0, 0])
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                                                                                                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
range_min = 0
range_max = 5
nbins = 10
hist_along_data = np.histogram(cov_spock_both_pole_arr/3600., range = [range_min, range_max], bins = nbins)
bin_array_temp = hist_along_data[1]
bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
binsize_actual = bin_array[1] - bin_array[0]
hist_along = hist_along_data[0] * 100. / len(cov_spock_both_pole_arr)
ax.bar(bin_array, hist_along, binsize_actual)
ax.set_ylim([0, np.max(hist_along)*1.1])
ax.set_xlim([np.min(bin_array_temp), np.max(bin_array_temp)])

fig.set_figheight(height_fig)
fig.set_figwidth(height_fig*ratio_fig_size)
fig_save_name = 'dist_nb_hours_both_poles.pdf'

fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')

raise Exception



### Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3


fig_title = 'Argument of apogee as a function of time'
y_label = 'Argument of apogee ' + u'(\N{DEGREE SIGN})' #'Real time'
x_label = 'Real time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in bold

for isc in range(nsc):
    ntle = len(arg_apogee_trans[isc])
    arg_apogee_trans_sc_arr = np.array(arg_apogee_trans[isc])
    in_pole_offset = np.where((np.abs(90 - arg_apogee_trans_sc_arr) <= pole_offset) | (np.abs(-90 - arg_apogee_trans_sc_arr) <= pole_offset))[0]
    out_pole_offset = np.where((np.abs(90 - arg_apogee_trans_sc_arr) > pole_offset) & (np.abs(-90 - arg_apogee_trans_sc_arr) > pole_offset))[0]
    nb_seconds_since_ref_sc_arr = np.array(nb_seconds_since_ref[isc])
    ax.scatter(nb_seconds_since_ref_sc_arr[in_pole_offset], arg_apogee_trans_sc_arr[in_pole_offset], s = 5, color = arr_info[isc][-1], label = arr_info[isc][0])
    ax.scatter(nb_seconds_since_ref_sc_arr[out_pole_offset], arg_apogee_trans_sc_arr[out_pole_offset], s = 5, color = arr_info[isc][-1], label = arr_info[isc][0])

    
# x axis label is in real time
xaxis_stop_date_date = datetime.strptime(xaxis_stop_date, "%Y-%m-%d")
delta_xaxis_start_stop = (xaxis_stop_date_date - xaxis_start_date_date).total_seconds() / 3600./ 24 / 30
nb_month_per_tick = (int)(delta_xaxis_start_stop / xaxis_nticks) + 1


xticks = np.zeros([xaxis_nticks+1])
date_list_str = []
xprevious = xaxis_start_date_date
date_list_str.append(datetime.strftime(xaxis_start_date_date, "%b") + '\n' + str(xaxis_start_date_date)[:4])
print nb_month_per_tick
for itick in range(xaxis_nticks-1):
    #ipdb.set_trace()
    xtick_temp = datetime.strftime(xprevious, "%Y-%m")
    xtick_temp_year = (int)(xtick_temp[:4])
    xtick_temp_month = (int)(xtick_temp[5:7])
    if xtick_temp_month + nb_month_per_tick > 12:
        xtick_temp_year_new = xtick_temp_year + (int)((xtick_temp_month + nb_month_per_tick)/12.)
        xtick_temp_month_new = np.mod(xtick_temp_month + nb_month_per_tick, 12)
        if xtick_temp_month_new == 0: # month 24, 36, 48, etc is Dec and the year should not be incremented by 2, 3, 4, etc but only 1, 2, 3, etc
            xtick_temp_month_new = 12
            xtick_temp_year_new = xtick_temp_year_new - 1
    else:
        xtick_temp_month_new = xtick_temp_month + nb_month_per_tick
        xtick_temp_year_new = xtick_temp_year
    xtick_temp_str = str(xtick_temp_year_new) + '-' + str(xtick_temp_month_new).zfill(2) + '-01'
    xtick_temp_date = datetime.strptime(xtick_temp_str, "%Y-%m-%d")
    xticks[itick+1] = (xtick_temp_date - xaxis_start_date_date).total_seconds() # xticks[0] is xaxis_start_date_date
    date_list_str.append(datetime.strftime(xtick_temp_date, "%b") + '\n' + xtick_temp_str[:4])
    print str(xprevious)[0:10], str(xtick_temp_date)[0:10]
    xprevious = xtick_temp_date
xticks[-1] = (xaxis_stop_date_date  - xaxis_start_date_date).total_seconds()
date_list_str.append(datetime.strftime(xaxis_stop_date_date, "%b") + '\n' + str(xaxis_stop_date_date)[:4])
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Satellite", fontsize = fontsize_plot,  handles=handles_arr)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'arg_apog_transform.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




fig_title = 'Argument of apogee as a function of time'
y_label = 'Argument of apogee ' + u'(\N{DEGREE SIGN})' #'Real time'
x_label = 'Real time'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in bold

for isc in range(nsc):
    ax.scatter(nb_seconds_since_ref[isc], arg_apogee[isc], s = 5, color = arr_info[isc][-1], label = arr_info[isc][0])

    
# x axis label is in real time
xaxis_stop_date_date = datetime.strptime(xaxis_stop_date, "%Y-%m-%d")
delta_xaxis_start_stop = (xaxis_stop_date_date - xaxis_start_date_date).total_seconds() / 3600./ 24 / 30
nb_month_per_tick = (int)(delta_xaxis_start_stop / xaxis_nticks) + 1


xticks = np.zeros([xaxis_nticks+1])
date_list_str = []
xprevious = xaxis_start_date_date
date_list_str.append(datetime.strftime(xaxis_start_date_date, "%b") + '\n' + str(xaxis_start_date_date)[:4])
print nb_month_per_tick
for itick in range(xaxis_nticks-1):
    #ipdb.set_trace()
    xtick_temp = datetime.strftime(xprevious, "%Y-%m")
    xtick_temp_year = (int)(xtick_temp[:4])
    xtick_temp_month = (int)(xtick_temp[5:7])
    if xtick_temp_month + nb_month_per_tick > 12:
        xtick_temp_year_new = xtick_temp_year + (int)((xtick_temp_month + nb_month_per_tick)/12.)
        xtick_temp_month_new = np.mod(xtick_temp_month + nb_month_per_tick, 12)
        if xtick_temp_month_new == 0: # month 24, 36, 48, etc is Dec and the year should not be incremented by 2, 3, 4, etc but only 1, 2, 3, etc
            xtick_temp_month_new = 12
            xtick_temp_year_new = xtick_temp_year_new - 1
    else:
        xtick_temp_month_new = xtick_temp_month + nb_month_per_tick
        xtick_temp_year_new = xtick_temp_year
    xtick_temp_str = str(xtick_temp_year_new) + '-' + str(xtick_temp_month_new).zfill(2) + '-01'
    xtick_temp_date = datetime.strptime(xtick_temp_str, "%Y-%m-%d")
    xticks[itick+1] = (xtick_temp_date - xaxis_start_date_date).total_seconds() # xticks[0] is xaxis_start_date_date
    date_list_str.append(datetime.strftime(xtick_temp_date, "%b") + '\n' + xtick_temp_str[:4])
    print str(xprevious)[0:10], str(xtick_temp_date)[0:10]
    xprevious = xtick_temp_date
xticks[-1] = (xaxis_stop_date_date  - xaxis_start_date_date).total_seconds()
date_list_str.append(datetime.strftime(xaxis_stop_date_date, "%b") + '\n' + str(xaxis_stop_date_date)[:4])
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Satellite", fontsize = fontsize_plot,  handles=handles_arr)
legend.get_title().set_fontsize(str(fontsize_plot))


fig_save_name = 'arg_apog.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


