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

# copy of cygnss_read_netcdf_and_convert_to_eci_and_oe.py on 04-22-2019.
# In addition to what is done cygnss_read_netcdf_and_convert_to_eci_and_oe.py, this script writes the eci r/v postion in the same format as the high drag eng pvt eci files (eg spock_FM5_20171216_eng_pvt_query-13527_eci.txt)


# This script is a copy of cygnss_read_netcdf.py on 2019/02/22. The only difference is that it also returns the ECI r/v and the osculating elements of the CYGNSS (from the ECEF of the netcdf file)
# This script reads a given netcdf file filename and return information on the satellite and SPs
# inputs: 
# - filename: name of netcdf file to read
# outputs:
# - lat_cyg, lon_cyg: latitude and longitude of cyg
# - lat_sp, lon_sp: latitude and longitude of sp
# - fom_sp, prn_sp: figure of merit (0-15) and prn of sp


# Assumptions:
# - see section PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
import ipdb
import matplotlib
from datetime import datetime, timedelta
import numpy as np
import os
import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cart2kep import *
from ecef2eci import *
from os import listdir
from read_input_file import *
from read_output_file import *
from cygnss_read_spock_spec import *
from cygnss_read_spock_spec_bin import *
from netCDF4 import Dataset
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
#from ecef2eci import *
#from eci_to_lvlh import *
import pickle
from cygnss_name_to_norad_id import *
import os.path

# load_spice_override: if 1 then ignored; if 0 then the spice kernels (to convert ecef to eci and orbital elements) won't be loaded, even though it's the first iteration of reading filename
def cygnss_read_netcdf_to_eci_observation(filename_list):
    # filename = '/Users/cbv/cygnss/netcdf/2018/013/cyg03.ddmi.s20180113-000000-e20180113-235959.l1.power-brcs.a21.d21.nc'
    # load_spice_override = 1
    nfile = len(filename_list)    
    for ifile in range(nfile):
        filename = filename_list[ifile]
        print ifile, nfile-1
        if ifile == 0:
            filename_eci  = filename.replace('.nc', '_eci.txt')
            file_eci = open(filename_eci, 'w')
            print >> file_eci,'#Date position(km/s) velocity(km/s)'
            print >> file_eci,'#START'
        time_gain_0 = []
        x_spec = []
        y_spec = []
        z_spec = []
        lat_spec = [] 
        quality_flags = [] 
        lon_spec = [] 

        x_cyg = []
        y_cyg = []
        z_cyg = []

        lat_cyg = []
        lon_cyg = []

        r_eci_cyg = []
        v_eci_cyg = []

        sma = []
        eccentricity = []
        inclination = []
        long_an = []
        w = []
        f = []
        period = []
        phase_angle = []


        pitch_cyg = []
        roll_cyg = []
        yaw_cyg = []

        x_gps = []
        y_gps = []
        z_gps = []

        vx_cyg = []
        vy_cyg = []
        vz_cyg = []


        gain = []
        az_spec = []
        el_spec = []
        az_orbit_spec = []
        el_orbit_spec = []

        fom = []
        index = []
        gps = []
        rx_to_sp_range = []
        tx_to_sp_range = []
        rcg = []
        date_flight = []
        date_flight_rounded = []
        date_flight_rounded_date = []
        index_in_spock_date_same = [] 
        index_in_spock_not_interpolated_date_same = []
        nb_seconds_since_initial_epoch_spock = []

        fh = Dataset(filename, mode='r')
        # nc_attrs = fh.ncattrs()
        # nc_dims = [dim for dim in fh.dimensions]  # list of nc dimensions
        # nc_vars = [var for var in fh.variables]  # list of nc variables

        x_spec_temp = fh.variables['sp_pos_x'][:] # X component of the specular point position in the ECEF coordinate system, in meters, at ddm_timestamp_utc, as calculated on the ground.
        y_spec_temp = fh.variables['sp_pos_y'][:]
        z_spec_temp = fh.variables['sp_pos_z'][:]

        lat_spec_temp = fh.variables['sp_lat'][:]
        lon_spec_temp = fh.variables['sp_lon'][:]

        x_cyg_temp = fh.variables['sc_pos_x'][:]
        y_cyg_temp = fh.variables['sc_pos_y'][:]
        z_cyg_temp= fh.variables['sc_pos_z'][:]

        lat_cyg_temp = fh.variables['sc_lat'][:]
        lon_cyg_temp = fh.variables['sc_lon'][:]

        pitch_cyg_temp = fh.variables['sc_pitch'][:] # Spacecraft pitch angle relative to the orbit frame, in radians at ddm_timestamp_utc
        roll_cyg_temp = fh.variables['sc_roll'][:]
        yaw_cyg_temp = fh.variables['sc_yaw'][:]

        x_gps_temp = fh.variables['tx_pos_x'][:]
        y_gps_temp = fh.variables['tx_pos_y'][:]
        z_gps_temp= fh.variables['tx_pos_z'][:]

        quality_flags_temp = fh.variables['quality_flags'][:]


        gain_temp = fh.variables['sp_rx_gain'][:] # The receive antenna gain in the direction of the specular point, in dBi, at ddm_timestamp_utc
        az_spec_temp = fh.variables['sp_az_body'][:] # Let line A be the line that extends from the spacecraft to the specular point, at ddm_timestamp_utc. Let line B be the projection of line A onto the spacecraft body frame XY plane. sp_az_body is the angle between the spacecraft body frame +X axis and line B, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems Definitions.
        el_spec_temp = fh.variables['sp_theta_body'][:] # The angle between the spacecraft body frame +Z axis and the line extending from the spacecraft to the specular point, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems Definitions.

        az_orbit_spec_temp = fh.variables['sp_az_orbit'][:] # Let line A be the line that extends from the spacecraft to the specular point at ddm_timestamp_utc. Let line B be the projection of line A onto the orbit frame XY plane. sp_az_orbit is the angle between the orbit frame +X axis (the velocity vector) and line B, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems D
        el_orbit_spec_temp = fh.variables['sp_theta_orbit'][:] # The angle between the orbit frame +Z axis and the line extending from the spacecraft to the specular point, in degrees, at ddm_timestamp_utc. See UM Doc. 148-0336, CYGNSS Science Data Processing Coordinate Systems Definitions.

        fom_temp = fh.variables['prn_fig_of_merit'][:] # The RCG Figure of Merit (FOM) for the DDM. Ranges from 0 through 15.The DDMI selects the four strongest xb.specular points (SP) for DDM production. It ranks the strength of SPs using an antenna RCG map. The map converts the position of the SP in antenna azimuth and declination angles to an RCG FOM. 0 represents the least FOM value. 15 represents the greatest FOM value.

        gps_temp = fh.variables['prn_code'][:] # The PRN code of the GPS signal associated with the DDM. Ranges from 0 to 32. 0 = reflectometry channel idle. 1 through 32 = GPS PRN codes.
        vx_cyg_temp = fh.variables['sc_vel_x'][:]
        vy_cyg_temp = fh.variables['sc_vel_y'][:]
        vz_cyg_temp= fh.variables['sc_vel_z'][:]

        rx_to_sp_range_temp = np.double(fh.variables['rx_to_sp_range'][:])
        tx_to_sp_range_temp = np.double(fh.variables['tx_to_sp_range'][:])

        time_flight = fh.variables['ddm_timestamp_utc'][:]
        time_coverage_start = fh.getncattr(fh.ncattrs()[fh.ncattrs().index('time_coverage_start')])
        time_coverage_start_datetime = datetime.strptime(time_coverage_start[:-4], "%Y-%m-%dT%H:%M:%S.%f") 
        #fh.close()
        nb_time_flight_temp = len(x_cyg_temp)
        date_flight_t = []
        date_flight_rounded_temp = []
        time_remove_list = []
        itime = -1
        prog = 0


        while itime < nb_time_flight_temp-1:
            if itime*100./(nb_time_flight_temp-1) > prog:
                #print str(prog)+'%'
                prog = prog + 10
            itime = itime + 1

            time_remove = 0
            date_flight_temp_date = time_coverage_start_datetime + timedelta(microseconds = round(time_flight[itime]*10**6))

            date_flight_temp = datetime.strftime(date_flight_temp_date, "%Y-%m-%dT%H:%M:%S.%f" )
            if ( date_flight_temp.split('.')[1][0] == '9' ): # round to next second
                date_flight_temp_date = datetime.strptime(date_flight_temp, "%Y-%m-%dT%H:%M:%S.%f")
                date_flight_date = date_flight_temp_date + timedelta(seconds = 1)
                date_flight_date_rounded_temp = datetime.strftime(date_flight_date, "%Y-%m-%dT%H:%M:%S.%f").split('.')[0]
                date_flight_date_rounded = datetime.strptime(date_flight_date_rounded_temp, "%Y-%m-%dT%H:%M:%S")
                date_flight_str_rounded = date_flight_date_rounded_temp
            elif ( date_flight_temp.split('.')[1][0] == '0' ): # round to next second
                date_flight_date_rounded = datetime.strptime(date_flight_temp.split('.')[0], "%Y-%m-%dT%H:%M:%S" )
                date_flight_str_rounded = date_flight_temp.split('.')[0]
            else: #if time can't be rounded by less than 100 ms
                time_remove = 1


            if ( time_remove == 1 ): # remove time if can't be rounded by ess than 100 ms 
                time_remove_list.append(itime)
            else:
                if type(x_spec_temp) == ma.core.MaskedArray:
                    x_spec.append(x_spec_temp.data[itime]/1000.)
                else:
                    x_spec.append(x_spec_temp[itime]/1000.)
                if type(y_spec_temp) == ma.core.MaskedArray:
                    y_spec.append(y_spec_temp.data[itime]/1000.)
                else:
                    y_spec.append(y_spec_temp[itime]/1000.)
                if type(z_spec_temp) == ma.core.MaskedArray:
                    z_spec.append(z_spec_temp.data[itime]/1000.)
                else:
                    z_spec.append(z_spec_temp[itime]/1000.)
                if type(gain_temp) == ma.core.MaskedArray:
                    gain.append(gain_temp.data[itime])
                else:
                    gain.append(gain_temp[itime])
                if type(az_spec_temp) == ma.core.MaskedArray:
                    az_spec.append(az_spec_temp.data[itime])
                else:
                    az_spec.append(az_spec_temp[itime])

                if type(el_spec_temp) == ma.core.MaskedArray:
                    el_spec.append(el_spec_temp.data[itime])
                else:
                    el_spec.append(el_spec_temp[itime])

                if type(az_orbit_spec_temp) == ma.core.MaskedArray:
                    az_orbit_spec.append(az_orbit_spec_temp.data[itime])
                else:
                    az_orbit_spec.append(az_orbit_spec_temp[itime])

                if type(el_orbit_spec_temp) == ma.core.MaskedArray:
                    el_orbit_spec.append(el_orbit_spec_temp.data[itime])
                else:
                    el_orbit_spec.append(el_orbit_spec_temp[itime])

                if type(fom_temp) == ma.core.MaskedArray:
                    fom.append(fom_temp.data[itime])
                else:
                    fom.append(fom_temp[itime])

                if type(gps_temp) == ma.core.MaskedArray:
                    gps.append(gps_temp.data[itime])
                else:
                    gps.append(gps_temp[itime])

                if type(rx_to_sp_range_temp) == ma.core.MaskedArray:
                    rx_to_sp_range.append(rx_to_sp_range_temp.data[itime])
                else:
                    rx_to_sp_range.append(rx_to_sp_range_temp[itime])
                if type(tx_to_sp_range_temp) == ma.core.MaskedArray:
                    tx_to_sp_range.append(tx_to_sp_range_temp.data[itime])
                else:
                    tx_to_sp_range.append(tx_to_sp_range_temp[itime])


                if type(lon_cyg_temp) == ma.core.MaskedArray:
                    lon_cyg.append(lon_cyg_temp.data[itime])
                else:
                    lon_cyg.append(lon_cyg_temp[itime])

                if type(lat_cyg_temp) == ma.core.MaskedArray:
                    lat_cyg.append(lat_cyg_temp.data[itime])
                else:
                    lat_cyg.append(lat_cyg_temp[itime])

                if type(lat_spec_temp) == ma.core.MaskedArray:
                    lat_spec.append(lat_spec_temp.data[itime])
                else:
                    lat_spec.append(lat_spec_temp[itime])
                if type(quality_flags_temp) == ma.core.MaskedArray:
                    quality_flags.append(quality_flags_temp.data[itime])
                else:
                    quality_flags.append(quality_flags_temp[itime])

                if type(lon_spec_temp) == ma.core.MaskedArray:
                    lon_spec.append(lon_spec_temp.data[itime])
                else:
                    lon_spec.append(lon_spec_temp[itime])

                if type(x_cyg_temp) == ma.core.MaskedArray:
                    x_cyg.append(x_cyg_temp.data[itime])
                else:
                    x_cyg.append(x_cyg_temp[itime])

                if type(y_cyg_temp) == ma.core.MaskedArray:
                    y_cyg.append(y_cyg_temp.data[itime])
                else:
                    y_cyg.append(y_cyg_temp[itime])

                if type(z_cyg_temp) == ma.core.MaskedArray:
                    z_cyg.append(z_cyg_temp.data[itime])
                else:
                    z_cyg.append(z_cyg_temp[itime])

                if type(vx_cyg_temp) == ma.core.MaskedArray:
                    vx_cyg.append(vx_cyg_temp.data[itime])
                else:
                    vx_cyg.append(vx_cyg_temp[itime])

                if type(vy_cyg_temp) == ma.core.MaskedArray:
                    vy_cyg.append(vy_cyg_temp.data[itime])
                else:
                    vy_cyg.append(vy_cyg_temp[itime])

                if type(vz_cyg_temp) == ma.core.MaskedArray:
                    vz_cyg.append(vz_cyg_temp.data[itime])
                else:
                    vz_cyg.append(vz_cyg_temp[itime])
                r_ecef = [x_cyg[-1]/1000., y_cyg[-1]/1000., z_cyg[-1]/1000.]
                v_ecef = [vx_cyg[-1]/1000., vy_cyg[-1]/1000., vz_cyg[-1]/1000.]
                if ((len(x_cyg) == 1) & (ifile == 0)):
                    load_spice = 1
                else:
                    load_spice = 0
                # if (load_spice_override == 0):
                #     load_spice = 0
                r_eci_temp, v_eci_temp = ecef2eci(r_ecef, v_ecef, date_flight_str_rounded, load_spice)
                r_eci_cyg.append(r_eci_temp*1000.)
                v_eci_cyg.append(v_eci_temp*1000.)
                sma_temp, eccentricity_temp, inclination_temp, long_an_temp,\
                    w_temp, f_temp, period_temp, phase_angle_temp\
                    = cart2kep(r_eci_temp, v_eci_temp, date_flight_str_rounded, load_spice)
                sma.append(sma_temp); eccentricity.append(eccentricity_temp);
                inclination.append(inclination_temp); long_an.append(long_an_temp);
                w.append(w_temp); f.append(f_temp); period.append(period_temp);
                phase_angle.append(phase_angle_temp)
                date_flight_rounded.append(date_flight_str_rounded)
                date_flight_rounded_date.append(date_flight_date_rounded)
                print >> file_eci, date_flight_rounded[-1], r_eci_temp[0], r_eci_temp[1], r_eci_temp[2], v_eci_temp[0], v_eci_temp[1], v_eci_temp[2]
    file_eci.close()
    ipdb.set_trace(); plt.ion(); fig, ax = plt.subplots(); ax.plot(np.abs(pitch_cyg_temp) * 180. / np.pi); fig, ax = plt.subplots(); ax.plot(np.abs(roll_cyg_temp) * 180. / np.pi); fig, ax = plt.subplots(); ax.plot(np.abs(yaw_cyg_temp) * 180. / np.pi);plt.show(); plt.show()
    return date_flight_rounded, lon_cyg, lat_cyg, lon_spec, lat_spec, fom, gps,\
        x_cyg, y_cyg, z_cyg, vx_cyg, vy_cyg, vz_cyg,date_flight_rounded_date,\
        r_eci_cyg, v_eci_cyg, sma, eccentricity, inclination, long_an, w, \
        f, period, phase_angle

# 
# fig, ax = plt.subplots(); ax.plot(np.abs(roll_cyg_temp) * 180. / np.pi)
# fig, ax = plt.subplots(); ax.plot(np.abs(yaw_cyg_temp) * 180. / np.pi)
