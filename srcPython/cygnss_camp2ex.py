# This script reads the lon/lat/alt and sc solar zenith angle (sza) from SpOCK and slecte the times when the sc is within 1000 km of the Clar air base in the Philippine (gs_lon, gs_lat). For those times, it outputs the lon/lat/alt and solar zenith angle
# Inputs:
# - input_filename: SpOCK main input filename
# Asssumptions:
# - see sections PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
input_filename = 'camp2exTest.txt'
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
from read_input_file import *
from read_output_file import *

gs_lon = 120.538067
gs_lat = 15.192336


# Read the position of the sc and the sza
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 

dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')];
date_start = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')];
date_stop = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')];
date_start_str = str(date_start).replace('-','')[:8]
date_stop_str = str(date_stop).replace('-','')[:8]
nb_sc = 8 
label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
label_arr_conversion = [3, 2, 7, 1, 0, 5, 6, 4]

for isc_temp in range(nb_sc): # !!!!!! nb_sc
    isc = label_arr_conversion[isc_temp]
    print label_arr[isc]
    var_to_read = ['longitude', 'latitude', 'altitude', 'solar_zenith']
    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
    longitude = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]
    latitude = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
    altitude = var_out[find_in_read_input_order_variables(var_out_order, 'altitude')]
    solar_zenith = var_out[find_in_read_input_order_variables(var_out_order, 'solar_zenith')]

    if isc_temp == 0:
        date = var_out[find_in_read_input_order_variables(var_out_order, 'date_round_sec')]
        nb_seconds_since_start = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
        nstep = len(date)

    # Compute the distance from the sc subsatellite point to the "ground station" (Clark air base)
    lon_sc = longitude * np.pi / 180
    lat_sc = latitude * np.pi / 180
    gs_lon_rad = gs_lon * np.pi / 180
    gs_lat_rad = gs_lat * np.pi / 180

    # approximate radius of earth in km at 15 deg latitude (source: https://rechneronline.de/earth-radius/)
    earth_radius = 6376.

    dlon = gs_lon_rad - lon_sc
    dlat = gs_lat_rad - lat_sc

    a = np.sin(dlat / 2)**2 + np.cos(lat_sc) * np.cos(gs_lat_rad) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    distance = earth_radius * c


    # Write the results in an output file: wheenver the subsatellite point is less than 1000 km from the air base, print the lon/lat/alt/sza
    output_filename = 'out/' + input_filename.replace('.txt', '_out_' + date_start_str + '_to_' + date_stop_str + '_' + label_arr[isc] + '.txt')
    output_file = open(output_filename, "w")
    print >> output_file, "# File auto-generated by the Spacecraft Orbital Characterization Kit (SpOCK)"
    print >> output_file, "# Position and solar zenith angle (SZA) of " + label_arr[isc] 
    print >> output_file, "# Lat/lon in deg\
    \n# Alt in km\
    \n# SZA in deg: 0 if the Sun is at the zenith of the sc\
    \n#             90 if the Sun is at the horizon of sc (aft)\
    \n#             180 if the center of the Earth is between the sc and the Sun\
    \n#             270 if the Sun is at the horizon of sc (fore)\
    \n#\
    \n#   TIME (UTC)       Lat     Lon     Alt    SZA"
    for istep in range(nstep):
        if distance[istep] < 1000.:
            print >> output_file, date[istep] + ' ' + format(latitude[istep], ".3f") + ' ' + format( longitude[istep], ".3f") + ' ' + format( altitude[istep], ".3f") + ' ' + format( solar_zenith[istep], ".3f")

    output_file.close()

