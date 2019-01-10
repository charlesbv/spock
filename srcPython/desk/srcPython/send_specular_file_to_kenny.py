import sys
import os
from get_prop_dir import *
from read_input_file import *


main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2] 

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6]; satellite_to_plot = input_variables[7]
nb_satellites = input_variables[4]; filename_gps_tle = input_variables[16]

name_mission = 'other/kenny/august/'

for isat in range(nb_satellites):
    os.system("rsync -av " + satellite_to_plot_path[isat] + "interpolated_position_LLA_ECEF_" + satellite_to_plot[isat] + " " + satellite_to_plot_path[isat] + "specular_" + satellite_to_plot[isat] + " " + main_input_file_name + " " + filename_gps_tle +  " srbwks2014-0008.engin.umich.edu:./" + name_mission)
