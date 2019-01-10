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

# This file is a copy of gps_read_spock_spec_bin.py on June 6 2018. It was written for the analysis with Scottt and CHris on sampling with all GNSS constellations (GPS, Galielo, Glonass, etc). The main difference with gps_read_spock_spec_bin.py is that less output is written if debug = 0  (lonl lat, and heading of sat are not printed by find_specular_points_gnssAll.c is debug = 0)
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
import numpy as np
from struct import *
from datetime import datetime, timedelta
from read_input_file import *


def gnss_read_spock_spec_bin(filename_spec, gps_name, dt_output, lonGrid, latGrid, debug):

    minLonGrid = lonGrid[0]
    maxLonGrid = lonGrid[1]
    dlon = lonGrid[2] # in degree
    minLatGrid = latGrid[0]
    maxLatGrid = latGrid[1]
    dlat = latGrid[2] # in degree
    nlon = (int)((maxLonGrid - minLonGrid)/dlon)
    nlat = (int)((maxLatGrid - minLatGrid)/dlat)
    grid_spec = np.zeros([nlat, nlon])
    grid_gps = np.zeros([nlat, nlon])

    #filename_spec = "/Users/cbv/cygnss/sift_temp/spock_out/test2_heading/test2_heading4/specular_test2_heading4.bin"    

    spec_file = open(filename_spec, "rb")

    lon_sat = [] ; lat_sat = [] ; heading_sat = [] ; date_now_all = []; ecef_x_spec = [] ; ecef_y_spec = [] ; ecef_z_spec = [] ; lon_spec = [] ; lat_spec = [] ; gain_spec = []; ecef_x_gps = []; ecef_y_gps = []; ecef_z_gps = []; ecef_vx_gps = []; ecef_vy_gps = []; ecef_vz_gps = []; ecef_x_sat = [] ; ecef_y_sat = [] ; ecef_z_sat = [] ; ecef_vx_sat = []; ecef_vy_sat = []; ecef_vz_sat = []; normpower_spec = []; elev_spec = [] ; elev_gps_from_cyg = []; elev_spec_not_int = []; azim_spec = [] ; azim_spec_not_int = []; prn = []

    iGps_temp = spec_file.read(4)
    count = 0
    while iGps_temp != "":

        iGps = unpack('i', iGps_temp)[0] 
        iPt = unpack('i', spec_file.read(4))[0]  # 0
        iPtInner  = iPt#( iPt % (int) (dt_output) );

        if count == 0:
            iPtInner_previous = iPtInner
            lon_sat_this_time = [] ; lat_sat_this_time = [] ; heading_sat_this_time = [] ;  ecef_x_spec_this_time = [] ; ecef_y_spec_this_time = [] ; ecef_z_spec_this_time = [] ; lon_spec_this_time = [] ; lat_spec_this_time = [] ; gain_spec_this_time = []; ecef_x_gps_this_time = []; ecef_y_gps_this_time = []; ecef_z_gps_this_time = []; ecef_vx_gps_this_time = []; ecef_vy_gps_this_time = []; ecef_vz_gps_this_time = []; ecef_x_sat_this_time = [] ; ecef_y_sat_this_time = [] ; ecef_z_sat_this_time = [] ; ecef_vx_sat_this_time = []; ecef_vy_sat_this_time = []; ecef_vz_sat_this_time = []; normpower_spec_this_time = []; elev_spec_this_time = [] ; elev_gps_from_cyg_this_time = []; elev_spec_not_int_this_time = []; azim_spec_this_time = [] ; azim_spec_not_int_this_time = []; prn_this_time = []

        count = 1

        if iPtInner != iPtInner_previous: # new time step
            lon_sat.append(lon_sat_this_time) ; lat_sat.append(lat_sat_this_time) ; heading_sat.append(heading_sat_this_time) ; date_now_all.append(date_now); ecef_x_spec.append(ecef_x_spec_this_time) ; ecef_y_spec.append(ecef_y_spec_this_time) ; ecef_z_spec.append(ecef_z_spec_this_time) ; lon_spec.append(lon_spec_this_time) ; lat_spec.append(lat_spec_this_time) ; gain_spec.append(gain_spec_this_time); ecef_x_gps.append(ecef_x_gps_this_time); ecef_y_gps.append(ecef_y_gps_this_time); ecef_z_gps.append(ecef_z_gps_this_time); ecef_vx_gps.append(ecef_vx_gps_this_time); ecef_vy_gps.append(ecef_vy_gps_this_time); ecef_vz_gps.append(ecef_vz_gps_this_time); ecef_x_sat.append(ecef_x_sat_this_time) ; ecef_y_sat.append(ecef_y_sat_this_time) ; ecef_z_sat.append(ecef_z_sat_this_time) ; ecef_vx_sat.append(ecef_vx_sat_this_time); ecef_vy_sat.append(ecef_vy_sat_this_time); ecef_vz_sat.append(ecef_vz_sat_this_time); normpower_spec.append(normpower_spec_this_time); elev_spec.append(elev_spec_this_time) ; elev_gps_from_cyg.append(elev_gps_from_cyg_this_time); elev_spec_not_int.append(elev_spec_not_int_this_time); azim_spec.append(azim_spec_this_time) ; azim_spec_not_int.append(azim_spec_not_int_this_time); prn.append(prn_this_time)

            lon_sat_this_time = [] ; lat_sat_this_time = [] ; heading_sat_this_time = [] ; ecef_x_spec_this_time = [] ; ecef_y_spec_this_time = [] ; ecef_z_spec_this_time = [] ; lon_spec_this_time = [] ; lat_spec_this_time = [] ; gain_spec_this_time = []; ecef_x_gps_this_time = []; ecef_y_gps_this_time = []; ecef_z_gps_this_time = []; ecef_vx_gps_this_time = []; ecef_vy_gps_this_time = []; ecef_vz_gps_this_time = []; ecef_x_sat_this_time = [] ; ecef_y_sat_this_time = [] ; ecef_z_sat_this_time = [] ; ecef_vx_sat_this_time = []; ecef_vy_sat_this_time = []; ecef_vz_sat_this_time = []; normpower_spec_this_time = []; elev_spec_this_time = [] ; elev_gps_from_cyg_this_time = []; elev_spec_not_int_this_time = []; azim_spec_this_time = [] ; azim_spec_not_int_this_time = []; prn_this_time =[]
            if ((iPt % (1440*30)) == 0):
                print iPt, iPt / (1440. * 30)
#        if (iPtInner  == 0):
        time_ymdhmsm  = np.zeros([7])
        for sss in range(7):
            time_ymdhmsm[sss]  = unpack('f', spec_file.read(4))[0]
        date_now_str_temp  = str((int)(time_ymdhmsm[0] )) + '-' + str((int)(time_ymdhmsm[1])) + '-' + str((int)(time_ymdhmsm[2])) + 'T' + str((int)(time_ymdhmsm[3])) + ':' + str((int)(time_ymdhmsm[4])) + ':' + str((int)(time_ymdhmsm[5]))
        date_now  = datetime.strptime(date_now_str_temp, "%Y-%m-%dT%H:%M:%S")
        date_now_save  = date_now
#         else:
#             date_now  = date_now_save + timedelta(seconds  = iPtInner)


        prn_this_time.append(gps_name[iGps])        
        if debug == 1:
            lon_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) 
            lat_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) 
            heading_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) 
            


            ecef_x_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # -5341.445801
            ecef_y_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # -1642.567505
            ecef_z_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # -3074.054443
        lon_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # 197.093384
        lat_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # -28.814657
        ilon = (int)((lon_spec_this_time[-1] - minLonGrid)/ dlon)
        ilat = (int)((lat_spec_this_time[-1] - minLatGrid)/ dlat)                
        if ( ilon == 360 ): # very rarely lon_spec_this_time[-1] = 360 instead of 0, for some reason
            ilon = 0
        if ( lon_spec_this_time[-1] >= 360.01 ):
            print lon_spec_this_time[-1], ilon, lat_spec_this_time[-1], ilat
        gain_spec_this_time.append(unpack('B', spec_file.read(1))[0] ) # 8
        if debug == 1:
            ecef_x_gps_this_time.append(unpack('f', spec_file.read(4))[0] ) # -22841.478516
            ecef_y_gps_this_time.append(unpack('f', spec_file.read(4))[0] ) # -11299.916992
            ecef_z_gps_this_time.append(unpack('f', spec_file.read(4))[0] ) # 7354.112793
            ecef_vx_gps_this_time.append(unpack('f', spec_file.read(4))[0] ) # 
            ecef_vy_gps_this_time.append(unpack('f', spec_file.read(4))[0] ) # 
            ecef_vz_gps_this_time.append(unpack('f', spec_file.read(4))[0] ) # 
            ecef_x_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) # -5451.068848
            ecef_y_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) # -1508.591064
            ecef_z_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) # -3941.083740
            ecef_vx_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) # 
            ecef_vy_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) # 
            ecef_vz_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) # 
            normpower_spec_this_time.append(unpack('B', spec_file.read(1))[0] ) # 8
            elev_spec_this_time.append(unpack('B', spec_file.read(1))[0] ) # 51
            azim_spec_temp  =  unpack('h', spec_file.read(2))[0]   # 281


            azim_spec_not_int_temp  = unpack('f', spec_file.read(4))[0]    # 280.973572

            if (azim_spec_temp < 0):
                azim_spec_temp  = 360 + azim_spec_temp;
            if (azim_spec_not_int_temp < 0):
                azim_spec_not_int_temp  = 360 + azim_spec_not_int_temp;

            azim_spec_this_time.append(azim_spec_temp)
            azim_spec_not_int_this_time.append(azim_spec_not_int_temp)


        elev_gps_from_cyg_this_time.append(unpack('f', spec_file.read(4))[0] ) # 24.299160
        elev_spec_not_int_this_time.append(unpack('f', spec_file.read(4))[0] ) # 51.279995
        if ((elev_spec_not_int_this_time[-1] <= 80.) & (gain_spec_this_time[-1] >= 2)):
            grid_spec[ilat, ilon] = grid_spec[ilat, ilon] + 1
        if elev_gps_from_cyg_this_time[-1] >= 10:
            grid_gps[ilat, ilon] = grid_gps[ilat, ilon] + 1



        iGps_temp = spec_file.read(4)
        iPtInner_previous = iPtInner
        #print  date_now_all, lon_spec, lat_spec, gain_spec, 'PRN_' + gps_name[iGps], normpower_spec, ecef_x_sat, ecef_y_sat, ecef_z_sat, ecef_x_gps, ecef_y_gps, ecef_z_gps, ecef_x_spec, ecef_y_spec, ecef_z_spec, elev_spec, azim_spec, elev_gps_from_cyg, elev_spec_not_int, azim_spec_not_int

    lon_sat.append(lon_sat_this_time) ; lat_sat.append(lat_sat_this_time) ; heading_sat.append(heading_sat_this_time) ; date_now_all.append(date_now); ecef_x_spec.append(ecef_x_spec_this_time) ; ecef_y_spec.append(ecef_y_spec_this_time) ; ecef_z_spec.append(ecef_z_spec_this_time) ; lon_spec.append(lon_spec_this_time) ; lat_spec.append(lat_spec_this_time) ; gain_spec.append(gain_spec_this_time); ecef_x_gps.append(ecef_x_gps_this_time); ecef_y_gps.append(ecef_y_gps_this_time); ecef_z_gps.append(ecef_z_gps_this_time); ecef_vx_gps.append(ecef_vx_gps_this_time); ecef_vy_gps.append(ecef_vy_gps_this_time); ecef_vz_gps.append(ecef_vz_gps_this_time); ecef_x_sat.append(ecef_x_sat_this_time) ; ecef_y_sat.append(ecef_y_sat_this_time) ; ecef_z_sat.append(ecef_z_sat_this_time) ; ecef_vx_sat.append(ecef_vx_sat_this_time); ecef_vy_sat.append(ecef_vy_sat_this_time); ecef_vz_sat.append(ecef_vz_sat_this_time); normpower_spec.append(normpower_spec_this_time); elev_spec.append(elev_spec_this_time) ; elev_gps_from_cyg.append(elev_gps_from_cyg_this_time); elev_spec_not_int.append(elev_spec_not_int_this_time); azim_spec.append(azim_spec_this_time) ; azim_spec_not_int.append(azim_spec_not_int_this_time); prn.append(prn_this_time)# for the last time step
    nb_spec_pts = len(lon_sat_this_time)
    spec_file.close()

    data = []
    if debug == 1:
        data.append(date_now_all); data.append(lon_spec); data.append(lat_spec); data.append(gain_spec); data.append(prn); data.append(normpower_spec); data.append(ecef_x_sat); data.append(ecef_y_sat); data.append(ecef_z_sat); data.append(ecef_x_gps); data.append(ecef_y_gps); data.append(ecef_z_gps); data.append( ecef_x_spec); data.append(ecef_y_spec); data.append(ecef_z_spec); data.append(nb_spec_pts); data.append(elev_spec); data.append(azim_spec); data.append(elev_gps_from_cyg); data.append(elev_spec_not_int); data.append(azim_spec_not_int); data.append(ecef_vx_sat); data.append(ecef_vy_sat); data.append(ecef_vz_sat); data.append(grid_spec); data.append(grid_gps)
    else:
        data.append(date_now_all); data.append(lon_spec); data.append(lat_spec); data.append(gain_spec); data.append(prn); data.append(grid_spec); data.append(grid_gps)
    return data
