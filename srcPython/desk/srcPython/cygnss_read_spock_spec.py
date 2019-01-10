# This function reads the specular point position files output by SpOCK. It returns the time/lon/lat/gain of the specular points as well as the GPS corresponding to the specular point. Also the NormPower of the spec, xyz of CYGNSS, xyz of GPS, xyz spec (ECEF, km)
# The arguments of the function are
# - filename_spec: the name of the specular file to read (path included)
# lon, lat, gain, gps are lists. For each list, each element of the list correspond to one time so it has nb_spec elements (where nb_spec is the number of specular points at that time (4 usually, if not always))

import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
import os
from read_input_file import *
from datetime import datetime, timedelta


def cygnss_read_spock_spec(filename_spec):
    #filename_spec = "./spock_out/spock_spec_start_2017-06-03T00_34_00_end_2017-06-03T00_50_00/spock_spec_start_2017-06-03T00_34_00_end_2017-06-03T00_50_004/specular_spock_spec_start_2017-06-03T00_34_00_end_2017-06-03T00_50_004.txt"

    # READ THE POSITIONS OF THE SPECULAR POINTS, THE GPS AND THE CYGNSS SATELLITES

    interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)

    count_less_than_4_spec = 0
    time_less_than_4_spec = []
    nb_spec_this_time = []
    file_specular = open(filename_spec)
    read_file_specular  = file_specular.readlines()
    # Nb of lines in the spec file header
    nb_lines_header_output_file_spec = 0
    while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
        nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
    nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1

    # determine number od spec points (usually 4 but this number can be changed in find_specular_points.c
    nb_spec_pts = 0
    date_temp = read_file_specular[nb_lines_header_output_file_spec].split()[0]
    date_temp_date = datetime.strptime(date_temp, "%Y-%m-%dT%H:%M:%S")
    date_temp_date_now = date_temp_date
    while (date_temp_date_now == date_temp_date):
        nb_spec_pts = nb_spec_pts + 1
        date_temp = read_file_specular[nb_lines_header_output_file_spec + nb_spec_pts].split()[0]
        date_temp_date_now = datetime.strptime(date_temp, "%Y-%m-%dT%H:%M:%S")



    ispec_save = 0
    j = -1
    lon = []; lat = []; gain = []; gps = []; date = []
    normpower = []; x_cyg = []; y_cyg = []; z_cyg = []; x_gps = []; y_gps = []; z_gps = []
    x_spec = []; y_spec = []; z_spec = []; elev_spec = []; azim_spec = []; elev_gps_from_cyg = []
    elev_spec_not_int = []; azim_spec_not_int = []; 
    order = []
    while (ispec_save < len(read_file_specular)-1-nb_spec_pts):
        lon_this_time = []
        lat_this_time = []
        gain_this_time = []
        gps_this_time = []
        normpower_this_time = []
        x_cyg_this_time = []
        y_cyg_this_time = []
        z_cyg_this_time = []
        x_gps_this_time = []
        y_gps_this_time = []
        z_gps_this_time = []
        x_spec_this_time = []
        y_spec_this_time = []
        z_spec_this_time = []
        elev_spec_this_time = []
        azim_spec_this_time = []
        elev_spec_not_int_this_time = []
        azim_spec_not_int_this_time = []
        elev_gps_from_cyg_this_time = []

        j = j + 1
        time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0]
        time_spec_sublist_temp_ini = datetime.strptime(time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
        date.append(time_spec_sublist_temp_ini)
        lon_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1]) )
        lat_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2]) )
        gain_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[3]) )
        gps_this_time.append( (int)(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[4].split('PRN_')[1]) )
        normpower_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[5]) )
        x_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[6]) )
        y_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[7]) )
        z_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[8]) )
        x_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[9]) )
        y_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[10]) )
        z_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[11]) )
        x_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[12]) )
        y_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[13]) )
        z_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[14]) )
        elev_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[15]) )
        azim_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[16]) )
        elev_gps_from_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[17]))
        elev_spec_not_int_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[18]) )
        azim_spec_not_int_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[19]) )



        ispec = 1
        while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == time_spec_sublist_temp_ini):
            lon_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1]) )
            lat_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2]) )
            gain_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3]) )
            gps_this_time.append( (int)(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[4].split('PRN_')[1]) )
            normpower_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[5]) )
            x_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[6]) )
            y_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[7]) )
            z_cyg_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[8]) )
            x_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[9]) )
            y_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[10]) )
            z_gps_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[11]) )
            x_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[12]) )
            y_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[13]) )
            z_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[14]) )
            elev_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[15]) )
            azim_spec_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[16]) )
            elev_gps_from_cyg_this_time.append(np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[17]))
            elev_spec_not_int_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[18]) )
            azim_spec_not_int_this_time.append( np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[19]) )


            ispec = ispec + 1
            if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular)):
                break
            if (len(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()) == 0): # sometimes the spec file ends with two blank lines
                break
        if ispec < 4:
            count_less_than_4_spec = count_less_than_4_spec + 1
            time_less_than_4_spec.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec-1].split()[0])
            nb_spec_this_time.append(ispec)
        ispec_save = ispec + ispec_save
        lon.append(lon_this_time)
        lat.append(lat_this_time)
        gain.append(gain_this_time)
        gps.append(gps_this_time)
        normpower.append(normpower_this_time)
        x_cyg.append(x_cyg_this_time)
        y_cyg.append(y_cyg_this_time)
        z_cyg.append(z_cyg_this_time)
        x_gps.append(x_gps_this_time)
        y_gps.append(y_gps_this_time)
        z_gps.append(z_gps_this_time)
        x_spec.append(x_spec_this_time)
        y_spec.append(y_spec_this_time)
        z_spec.append(z_spec_this_time)
        elev_spec.append(elev_spec_this_time)
        azim_spec.append(azim_spec_this_time)
        elev_gps_from_cyg.append(elev_gps_from_cyg_this_time)
        elev_spec_not_int.append(elev_spec_not_int_this_time)
        azim_spec_not_int.append(azim_spec_not_int_this_time)
        order_this_time = np.zeros([nb_spec_pts])
        index_sort_normpower = np.argsort(normpower_this_time)
        order_now = 0
        ispec_order_now = 0
        order_this_time[index_sort_normpower[ispec_order_now]] = 0
        for ispec_order in range(1,nb_spec_pts):
            if normpower_this_time[index_sort_normpower[ispec_order]] == normpower_this_time[index_sort_normpower[ispec_order-1]] :
                order_this_time[index_sort_normpower[ispec_order]] = order_now
            else:
                order_now = order_now + 1
                order_this_time[index_sort_normpower[ispec_order]] = order_now
        order.append(order_this_time)
    file_specular.close() 
    return date, lon, lat, gain, gps, normpower, x_cyg, y_cyg, z_cyg, x_gps, y_gps, z_gps,  x_spec, y_spec, z_spec, nb_spec_pts, elev_spec, azim_spec, elev_gps_from_cyg, elev_spec_not_int, azim_spec_not_int



