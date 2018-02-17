import numpy as np
from datetime import datetime, timedelta
from get_prop_dir import *

def out(run_dir, filename, which_sat, list_variable): # this is a new version of formerly called read_output_file: run_dir is the name of the run directory (no path needed); filename is the name of the propagator input file (no path needed); which_sat is the satellite to compute among all satellites propagated in the propagtor input file (which_sat starts at 1)
# filename = '/home/cbv/PropSim/output/run_cygnss1_one_year_in_a_row/cygnss2_day1/cygnss2_day1.txt'
# list_variable = ['position',
#                  'velocity',
#                  'longitude',
#                  'latitude',
#                  'altitude',
#                  'raan',
#                  'true_anomaly',
#                  'arg_perigee',
#                  'right_asc',
#                  'local_time']

    # OUTPUT FILE NAMES
    ## Open propgaator input file to find the name of the propagator output files
    filename = get_prop_dir(1) + run_dir + '/' + 'input/main_input/' + filename
    file_input = open(filename, "r")
    read_file_input = file_input.readlines()
    ### NUMBER OF SATELLITES
    nb_satellites =  (int)(read_file_input[6].split()[0]) 
    output_file_name_list = []
    output_file_path_list = []
    for i in range(nb_satellites):
        name_with_extension = (read_file_input[8].split(',')[i]).replace("\n","").replace(" ","") 
        name_no_extension = (name_with_extension.split('.')[0]).replace("\n","").replace(" ","") 
        output_file_path_list.append(get_prop_dir(1) + run_dir + "/output/run_" + ((read_file_input[8].split(', ')[0]).split('.')[0]).replace("\n","") + "/" + name_no_extension + "/")
        output_file_name_list.append(name_with_extension)
    filename_output = output_file_path_list[which_sat - 1] + output_file_name_list[which_sat - 1]

    file_to_read = open(filename_output, "r")   
    read_file_to_read = file_to_read.readlines()
    nb_lines_header = 10
    n = len(read_file_to_read) - nb_lines_header
    nb_variables = len(list_variable)
    date = []
    position = np.zeros([n, 3])
    velocity = np.zeros([n, 3])
    longitude = np.zeros([n])
    latitude = np.zeros([n])
    altitude = np.zeros([n])
    sma = np.zeros([n])
    inclination = np.zeros([n])
    eccentricity = np.zeros([n])
    true_anomaly = np.zeros([n])
    raan = np.zeros([n])
    arg_perigee = np.zeros([n])
    right_asc = np.zeros([n])
    local_time = np.zeros([n])
    density = np.zeros([n])
    calculate_position = 0
    calculate_velocity = 0
    calculate_longitude = 0
    calculate_latitude = 0
    calculate_altitude = 0
    calculate_sma = 0
    calculate_inclination = 0
    calculate_eccentricity = 0
    calculate_true_anomaly = 0
    calculate_raan = 0
    calculate_arg_perigee = 0
    calculate_right_asc = 0
    calculate_local_time = 0
    calculate_density = 0
    for j in range(nb_variables):
        if ( list_variable[j] == 'position' ):
            calculate_position = 1
        if ( list_variable[j] == 'velocity' ):
            calculate_velocity = 1
        if ( list_variable[j] == 'longitude' ):
            calculate_longitude = 1
        if ( list_variable[j] == 'latitude' ):
            calculate_latitude = 1
        if ( list_variable[j] == 'altitude' ):
            calculate_altitude = 1
        if ( list_variable[j] == 'sma' ):
            calculate_sma = 1
        if ( list_variable[j] == 'inclination' ):
            calculate_inclination = 1
        if ( list_variable[j] == 'eccentricity' ):
            calculate_eccentricity = 1
        if ( list_variable[j] == 'true_anomaly' ):
            calculate_true_anomaly = 1
        if ( list_variable[j] == 'raan' ):
            calculate_raan = 1
        if ( list_variable[j] == 'arg_perigee' ):
            calculate_arg_perigee = 1
        if ( list_variable[j] == 'right_asc' ):
            calculate_right_asc = 1
        if ( list_variable[j] == 'local_time' ):
            calculate_local_time = 1
        if ( list_variable[j] == 'density' ):
            calculate_density = 1


    for i in range(n):
        date.append(read_file_to_read[i+nb_lines_header].split()[0] + " " + read_file_to_read[i+nb_lines_header].split()[1])
        # if i == 1:
        #     dt_here = ( datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S") - datetime.strptime(date[-2], "%Y/%m/%d %H:%M:%S") ).days * 24 * 3600 + ( datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S") - datetime.strptime(date[-2], "%Y/%m/%d %H:%M:%S") ).seconds; 
        # if i > 0:
        #     if ( datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S") - datetime.strptime(date[-2], "%Y/%m/%d %H:%M:%S") ).days * 24 * 3600 + ( datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S") - datetime.strptime(date[-2], "%Y/%m/%d %H:%M:%S") ).seconds != dt_here:
        #         print date[-1]
        if (calculate_position == 1):
            position[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[2])
            position[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[3])
            position[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[4])
        if (calculate_velocity == 1):
            velocity[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[5])
            velocity[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[6])
            velocity[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[7])
        if (calculate_longitude == 1):
            longitude[i] = np.float(read_file_to_read[i+nb_lines_header].split()[8])
        if (calculate_latitude == 1):
            latitude[i] = np.float(read_file_to_read[i+nb_lines_header].split()[9])
        if (calculate_altitude == 1):
            altitude[i] = np.float(read_file_to_read[i+nb_lines_header].split()[10])
        if (calculate_sma == 1):
            sma[i] = np.float(read_file_to_read[i+nb_lines_header].split()[11])
        if (calculate_inclination == 1):
            inclination[i] = np.float(read_file_to_read[i+nb_lines_header].split()[12])
        if (calculate_eccentricity == 1):
            eccentricity[i] = np.float(read_file_to_read[i+nb_lines_header].split()[13])
        if (calculate_true_anomaly == 1):
            true_anomaly[i] = np.float(read_file_to_read[i+nb_lines_header].split()[14])
        if (calculate_raan == 1):
            raan[i] = np.float(read_file_to_read[i+nb_lines_header].split()[15])
        if (calculate_arg_perigee == 1):
            arg_perigee[i] = np.float(read_file_to_read[i+nb_lines_header].split()[16])
        if (calculate_right_asc == 1):
            right_asc[i] = np.float(read_file_to_read[i+nb_lines_header].split()[17])
        if (calculate_local_time == 1):
            local_time[i] = np.float(read_file_to_read[i+nb_lines_header].split()[18])

    variables = []
    order_variables = []
    variables.append(date)
    order_variables.append("date | " + str(len(order_variables)))
    if (calculate_position == 1):
        variables.append(position)
        order_variables.append("position | " + str(len(order_variables)))
    if (calculate_velocity == 1):
        variables.append(velocity)
        order_variables.append("velocity | " + str(len(order_variables)))
    if (calculate_longitude == 1):
        variables.append(longitude)
        order_variables.append("longitude | " + str(len(order_variables)))
    if (calculate_latitude == 1):
        variables.append(latitude)
        order_variables.append("latitude | " + str(len(order_variables)))
    if (calculate_altitude == 1):
        variables.append(altitude)
        order_variables.append("altitude | " + str(len(order_variables)))
    if (calculate_sma == 1):
        variables.append(sma)
        order_variables.append("sma | " + str(len(order_variables)))
    if (calculate_inclination == 1):
        variables.append(inclination)
        order_variables.append("inclination | " + str(len(order_variables)))
    if (calculate_eccentricity == 1):
        variables.append(eccentricity)
        order_variables.append("eccentricity | " + str(len(order_variables)))
    if (calculate_true_anomaly == 1):
        variables.append(true_anomaly)
        order_variables.append("true_anomaly | " + str(len(order_variables)))
    if (calculate_raan == 1):
        variables.append(raan)
        order_variables.append("raan | " + str(len(order_variables)))
    if (calculate_arg_perigee == 1):
        variables.append(arg_perigee)
        order_variables.append("arg_perigee | " + str(len(order_variables)))
    if (calculate_right_asc == 1):
        variables.append(right_asc)
        order_variables.append("right_asc | " + str(len(order_variables)))
    if (calculate_local_time == 1):
        variables.append(local_time)
        order_variables.append("local_time | " + str(len(order_variables)))

    if ( calculate_density == 1 ):
        density_filename = ""
        for j in range(1,len(filename_output.split('/')[0:-1])):
            density_filename = density_filename + "/" + filename_output.split('/')[0:-1][j]                
        density_filename = density_filename + "/density_" + filename_output.split('/')[-1]
        density_file = open(density_filename, "r")
        density_file_read = density_file.readlines()
        for i in range(1,n):
            density[i] = np.float(density_file_read[i-1].split()[1]) / 10**9
        density[0] = density[1] # since the ouptut file misses the first time step of the similution, we just make this density equal to the one at the second time step
        density_file.close()
        variables.append(density)
        order_variables.append("density | " + str(len(order_variables)))

    return variables, order_variables
