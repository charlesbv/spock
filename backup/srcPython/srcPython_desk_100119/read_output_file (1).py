import numpy as np
from datetime import datetime, timedelta

def read_output_file(filename, list_variable):
#     filename =  'spock/minxss_1_0_eci/minxss_1_0_eci1/minxss_1_0_eci1.txt'
#     list_variable = ['power']
# list_variable = ['position',
#                  'velocity',
#                  'longitude',
#                  'latitude',
#                  'altitude',
#                  'raan',
#                  'true_anomaly',
#                  'argument_perigee',
#                  'right_asc',
#                  'local_time']


    file_to_read = open(filename, "r")   
    read_file_to_read = file_to_read.readlines()
    nb_lines_header = 10
    n = len(read_file_to_read) - nb_lines_header
    nb_variables = len(list_variable)
    if 'given_output' in list_variable:
        file_given_output =  open('/'.join(filename.split('/')[:-1]) + "/given_output_" + filename.split('/')[-1])
        read_file_given_output = file_given_output.readlines() # !!!!! first date of given output is not necesarrily the same as for the other variables 
    if 'power' in list_variable:
        file_power =  open('/'.join(filename.split('/')[:-1]) + "/power_" + filename.split('/')[-1])
        read_file_power = file_power.readlines() # !!!!! first date of given output is not necesarrily the same as for the other variables 
    date = []
    position = np.zeros([n, 3])
    radius = np.zeros([n])
    speed = np.zeros([n])
    velocity = np.zeros([n, 3])
    acceleration = np.zeros([n, 3])
    acceleration_lvlh = np.zeros([n, 3])
    acceleration_lvlh_gravity = np.zeros([n, 3])
    acceleration_lvlh_drag = np.zeros([n, 3])
    acceleration_eci_drag = np.zeros([n, 3])
    longitude = np.zeros([n])
    latitude = np.zeros([n])
    altitude = np.zeros([n])
    sma = np.zeros([n])
    inclination = np.zeros([n])
    eccentricity = np.zeros([n])
    true_anomaly = np.zeros([n])
    raan = np.zeros([n])
    given_output = [] # !!!!! first date of given output is not necesarrily the same as for the other variables 
    power = np.zeros([n]) # !!!!! first date of given output is not necesarrily the same as for the other variables 
    argument_perigee = np.zeros([n])
    right_asc = np.zeros([n])
    local_time = np.zeros([n])
    density = np.zeros([n])
    temperature = np.zeros([n])
    cd = np.zeros([n])
    tot_area_drag = np.zeros([n])
    radius_perigee = np.zeros([n])
    radius_apogee = np.zeros([n])
    calculate_position = 0
    calculate_radius = 0
    calculate_speed = 0
    calculate_velocity = 0
    calculate_acceleration = 0
    calculate_acceleration_lvlh = 0
    calculate_acceleration_lvlh_gravity = 0
    calculate_acceleration_lvlh_drag = 0
    calculate_acceleration_eci_drag = 0
    calculate_longitude = 0
    calculate_latitude = 0
    calculate_altitude = 0
    calculate_sma = 0
    calculate_inclination = 0
    calculate_eccentricity = 0
    calculate_true_anomaly = 0
    calculate_raan = 0
    calculate_given_output = 0
    calculate_power = 0
    calculate_argument_perigee = 0
    calculate_right_asc = 0
    calculate_local_time = 0
    calculate_density = 0
    calculate_temperature = 0
    calculate_cd = 0
    calculate_tot_area_drag = 0
    calculate_radius_perigee = 0
    calculate_radius_apogee = 0
    for j in range(nb_variables):
        if ( list_variable[j] == 'position' ):
            calculate_position = 1
        if ( list_variable[j] == 'radius' ):
            calculate_radius = 1
            calculate_position = 1
        if ( list_variable[j] == 'velocity' ):
            calculate_velocity = 1
        if ( list_variable[j] == 'speed' ):
            calculate_speed = 1
            calculate_velocity = 1
        if ( list_variable[j] == 'acceleration' ):
            calculate_acceleration = 1
        if ( list_variable[j] == 'acceleration_lvlh' ):
            calculate_acceleration_lvlh = 1
        if ( list_variable[j] == 'acceleration_lvlh_gravity' ):
            calculate_acceleration_lvlh_gravity = 1
        if ( list_variable[j] == 'acceleration_lvlh_drag' ):
            calculate_acceleration_lvlh_drag = 1
        if ( list_variable[j] == 'acceleration_eci_drag' ):
            calculate_acceleration_eci_drag = 1
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
        if ( list_variable[j] == 'given_output' ):
            calculate_given_output = 1
        if ( list_variable[j] == 'power' ):
            calculate_power = 1
        if ( list_variable[j] == 'argument_perigee' ):
            calculate_argument_perigee = 1
        if ( list_variable[j] == 'right_asc' ):
            calculate_right_asc = 1
        if ( list_variable[j] == 'local_time' ):
            calculate_local_time = 1
        if ( list_variable[j] == 'density' ):
            calculate_density = 1
        if ( list_variable[j] == 'temperature' ):
            calculate_temperature = 1
        if ( list_variable[j] == 'cd' ):
            calculate_cd = 1
        if ( list_variable[j] == 'tot_area_drag' ):
            calculate_tot_area_drag = 1
        if ( list_variable[j] == 'radius_perigee' ):
            calculate_radius_perigee = 1
        if ( list_variable[j] == 'radius_apogee' ):
            calculate_radius_apogee = 1


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
        if (calculate_radius == 1):
            radius[i] = np.sqrt( position[i,0]*position[i,0] + position[i,1]*position[i,1] + position[i,2]*position[i,2] )
        if (calculate_velocity == 1):
            velocity[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[5])
            velocity[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[6])
            velocity[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[7])
        if (calculate_speed == 1):
            speed[i] = np.sqrt( velocity[i,0]*velocity[i,0] + velocity[i,1]*velocity[i,1] + velocity[i,2]*velocity[i,2] )

        if (calculate_longitude == 1):
            longitude[i] = np.float(read_file_to_read[i+nb_lines_header].split()[8])
        if (calculate_latitude == 1):
            latitude[i] = np.float(read_file_to_read[i+nb_lines_header].split()[9])
        if (calculate_altitude == 1):
            altitude[i] = np.float(read_file_to_read[i+nb_lines_header].split()[10])
        if ( (calculate_radius_perigee == 1) | (calculate_radius_apogee == 1) ):
            calculate_sma = 1
            calculate_eccentricity = 1
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
        if (calculate_power == 1):
            power[i] = np.sum(np.array(read_file_power[i+13].split()[2:-2]).astype(np.float)) # !!!!! first date of given output is not necesarrily the same as for the other variables
        if (calculate_given_output == 1):
            given_output.append(read_file_given_output[i].split()[-1]) # !!!!! first date of given output is not necesarrily the same as for the other variables
        if (calculate_argument_perigee == 1):
            argument_perigee[i] = np.float(read_file_to_read[i+nb_lines_header].split()[16])
        if (calculate_right_asc == 1):
            right_asc[i] = np.float(read_file_to_read[i+nb_lines_header].split()[17])
        if (calculate_local_time == 1):
            local_time[i] = np.float(read_file_to_read[i+nb_lines_header].split()[18])
        if (calculate_acceleration == 1):
            acceleration[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[19])
            acceleration[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[20])
            acceleration[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[21])
        if (calculate_acceleration_lvlh == 1):
            acceleration_lvlh[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[22])
            acceleration_lvlh[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[23])
            acceleration_lvlh[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[24])
        if (calculate_acceleration_lvlh_gravity == 1):
            acceleration_lvlh_gravity[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[25])
            acceleration_lvlh_gravity[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[26])
            acceleration_lvlh_gravity[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[27])
        if (calculate_acceleration_lvlh_drag == 1):
            acceleration_lvlh_drag[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[28])
            acceleration_lvlh_drag[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[29])
            acceleration_lvlh_drag[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[30])
        if (calculate_acceleration_eci_drag == 1):
            acceleration_eci_drag[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[31])
            acceleration_eci_drag[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[32])
            acceleration_eci_drag[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[33])

        if (calculate_radius_perigee == 1):
            radius_perigee[i] = sma[i]*(1 - eccentricity[i])
        if (calculate_radius_apogee == 1):
            radius_apogee[i] = sma[i]*(1 + eccentricity[i])
    variables = []
    order_variables = []
    variables.append(date)
    order_variables.append("date | " + str(len(order_variables)))
    if (calculate_position == 1):
        variables.append(position)
        order_variables.append("position | " + str(len(order_variables)))
    if (calculate_radius == 1):
        variables.append(radius)
        order_variables.append("radius | " + str(len(order_variables)))
    if (calculate_velocity == 1):
        variables.append(velocity)
        order_variables.append("velocity | " + str(len(order_variables)))
    if (calculate_speed == 1):
        variables.append(speed)
        order_variables.append("speed  | " + str(len(order_variables)))

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
    if (calculate_given_output == 1): # !!!!! first date of given output is not necesarrily the same as for the other variables 
        variables.append(given_output)
        order_variables.append("given_output | " + str(len(order_variables)))
    if (calculate_power == 1): # !!!!! first date of given output is not necesarrily the same as for the other variables 
        variables.append(power)
        order_variables.append("power | " + str(len(order_variables)))

    if (calculate_argument_perigee == 1):
        variables.append(argument_perigee)
        order_variables.append("argument_perigee | " + str(len(order_variables)))
    if (calculate_right_asc == 1):
        variables.append(right_asc)
        order_variables.append("right_asc | " + str(len(order_variables)))
    if (calculate_local_time == 1):
        variables.append(local_time)
        order_variables.append("local_time | " + str(len(order_variables)))
    if (calculate_acceleration == 1):
        variables.append(acceleration)
        order_variables.append("acceleration | " + str(len(order_variables)))
    if (calculate_acceleration_lvlh == 1):
        variables.append(acceleration_lvlh)
        order_variables.append("acceleration_lvlh | " + str(len(order_variables)))
    if (calculate_acceleration_lvlh_gravity == 1):
        variables.append(acceleration_lvlh_gravity)
        order_variables.append("acceleration_lvlh_gravity | " + str(len(order_variables)))
    if (calculate_acceleration_lvlh_drag == 1):
        variables.append(acceleration_lvlh_drag)
        order_variables.append("acceleration_lvlh_drag | " + str(len(order_variables)))
    if (calculate_acceleration_eci_drag == 1):
        variables.append(acceleration_eci_drag)
        order_variables.append("acceleration_eci_drag | " + str(len(order_variables)))


    if ( (calculate_density == 1 ) | (calculate_temperature == 1) | (calculate_cd == 1) | (calculate_tot_area_drag == 1) ):
        #        density_filename = ""
        #         for j in range(1,len(filename.split('/')[0:-1])):
        #             density_filename = density_filename + "/" + filename.split('/')[0:-1][j]         
        #         density_filename = density_filename + "/density_" + filename.split('/')[-1]

        density_filename = '/'.join(filename.split('/')[:-1]) + "/density_" + filename.split('/')[-1]
        density_file = open(density_filename, "r")
        density_file_read = density_file.readlines()
        for i in range(1,n):
            if calculate_density == 1:
                density[i] = np.float(density_file_read[i-1].split()[1]) / 10**9 # in kg/m^3
            if calculate_temperature == 1:
                temperature[i] = np.float(density_file_read[i-1].split()[2]) # in K
            if calculate_cd == 1:
                cd[i] = np.float(density_file_read[i-1].split()[3]) # in K
            if calculate_tot_area_drag == 1:
                tot_area_drag[i] = np.float(density_file_read[i-1].split()[4]) # in K

        if calculate_density == 1:
            density[0] = density[1] # since the ouptut file misses the first time step of the similution, we just make this density equal to the one at the second time step
        if calculate_temperature ==1:
            temperature[0] = temperature[1] # since the ouptut file misses the first time step of the similution, we just make this temperature equal to the one at the second time step
        if calculate_cd ==1:
            cd[0] = cd[1] # since the ouptut file misses the first time step of the similution, we just make this cd equal to the one at the second time step
        if calculate_tot_area_drag ==1:
            tot_area_drag[0] = tot_area_drag[1] # since the ouptut file misses the first time step of the similution, we just make this tot_area_drag equal to the one at the second time step

        density_file.close()
        if calculate_density == 1:
            variables.append(density)
            order_variables.append("density | " + str(len(order_variables)))
        if calculate_temperature == 1:
            variables.append(temperature)
            order_variables.append("temperature | " + str(len(order_variables)))
        if calculate_cd == 1:
            variables.append(cd)
            order_variables.append("cd | " + str(len(order_variables)))
        if calculate_tot_area_drag == 1:
            variables.append(tot_area_drag)
            order_variables.append("tot_area_drag | " + str(len(order_variables)))

            
    if (calculate_radius_perigee == 1):
        variables.append(radius_perigee)
        order_variables.append("radius_perigee | " + str(len(order_variables)))

    if (calculate_radius_apogee == 1):
        variables.append(radius_apogee)
        order_variables.append("radius_apogee | " + str(len(order_variables)))


    return variables, order_variables
