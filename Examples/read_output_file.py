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

import numpy as np
from datetime import datetime, timedelta
import ipdb
def read_output_file(filename, list_variable):


    file_to_read = open(filename, "r")   
    read_file_to_read = file_to_read.readlines()
    nb_lines_header = 10
    n = len(read_file_to_read) - nb_lines_header
    nb_variables = len(list_variable)
    if (('ap' in list_variable) | ('f107'  in list_variable)):
        file_given_output =  open('/'.join(filename.split('/')[:-1]) + "/given_output_" + filename.split('/')[-1])
        read_file_given_output = file_given_output.readlines() # !!!!! first date of given output is not necesarrily the same as for the other variables 
    if 'power' in list_variable:
        file_power =  open('/'.join(filename.split('/')[:-1]) + "/power_" + filename.split('/')[-1])
        read_file_power = file_power.readlines() # !!!!! first date of given output is not necesarrily the same as for the other variables 
    date = []
    date_datetime = []
    date_round_sec = [] 
    date_datetime_round_sec = []
    orb_number = []
    previous_orb = 0
    position = np.zeros([n, 3])
    position_ecef = np.zeros([n, 3])
    velocity_ecef = np.zeros([n, 3])
    radius = np.zeros([n])
    speed = np.zeros([n])
    velocity = np.zeros([n, 3])
    acceleration = np.zeros([n, 3])
    acceleration_lvlh = np.zeros([n, 3])
    acceleration_lvlh_gravity = np.zeros([n, 3])
    acceleration_lvlh_drag = np.zeros([n, 3])
    acceleration_lvlh_drag_mag = np.zeros([n])
    acceleration_eci_drag = np.zeros([n, 3])
    longitude = np.zeros([n])
    nb_seconds_since_start = np.zeros([n])
    latitude = np.zeros([n])
    altitude = np.zeros([n])
    sma = np.zeros([n])
    sma_ave = np.zeros([n])
    inclination = np.zeros([n])
    eccentricity = np.zeros([n])
    solar_zenith = np.zeros([n])
    true_anomaly = np.zeros([n])
    phase_angle = np.zeros([n]) # true_anomaly + argument_perigee
    phase_rate = np.zeros([n]) # ( phase_angle[i+1] - phase_angle[i] ) / dt (0 for i = 0)
    raan = np.zeros([n])
    beta = np.zeros([n]) # beta angle (orbit plane, Sun)
    ap = [] # !!!!! first date of given output is not necesarrily the same as for the other variables 
    f107 = [] # !!!!! first date of given output is not necesarrily the same as for the other variables 
    date_given_output = []
    power = np.zeros([n]) # !!!!! first date of given output is not necesarrily the same as for the other variables 
    argument_perigee = np.zeros([n])
    argument_perigee_ave = np.zeros([n])
    right_asc = np.zeros([n])
    local_time = np.zeros([n])
    density = np.zeros([n])
    temperature = np.zeros([n])
    cd = np.zeros([n])
    tot_area_drag = np.zeros([n])
    radius_perigee = np.zeros([n])
    radius_apogee = np.zeros([n])
    calculate_position = 0
    calculate_position_ecef = 0
    calculate_velocity_ecef = 0
    calculate_position_tle = 0
    calculate_velocity_tle = 0

    calculate_radius = 0
    calculate_speed = 0
    calculate_velocity = 0
    calculate_acceleration = 0
    calculate_acceleration_lvlh = 0
    calculate_acceleration_lvlh_gravity = 0
    calculate_acceleration_lvlh_drag = 0
    calculate_acceleration_lvlh_drag_mag = 0
    calculate_acceleration_eci_drag = 0
    calculate_longitude = 0
    calculate_nb_seconds_since_start = 0
    calculate_latitude = 0
    calculate_altitude = 0
    calculate_sma = 0
    calculate_sma_ave = 0
    calculate_inclination = 0
    calculate_eccentricity = 0
    calculate_solar_zenith = 0
    calculate_true_anomaly = 0
    calculate_phase_angle = 0
    calculate_phase_rate = 0
    calculate_raan = 0
    calculate_beta = 0
    calculate_ap = 0
    calculate_f107 = 0
    calculate_power = 0
    calculate_argument_perigee = 0
    calculate_argument_perigee_ave = 0
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
        if ( list_variable[j] == 'position_ecef' ):
            calculate_position_ecef = 1
        if ( list_variable[j] == 'velocity_ecef' ):
            calculate_velocity_ecef = 1
        if ( list_variable[j] == 'position_tle' ):
            calculate_position_tle = 1
        if ( list_variable[j] == 'velocity_tle' ):
            calculate_velocity_tle = 1

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
        if ( list_variable[j] == 'acceleration_lvlh_drag_mag' ):
            calculate_acceleration_lvlh_drag_mag = 1
            calculate_acceleration_lvlh_drag = 1
        if ( list_variable[j] == 'acceleration_eci_drag' ):
            calculate_acceleration_eci_drag = 1
        if ( list_variable[j] == 'longitude' ):
            calculate_longitude = 1
        if ( list_variable[j] == 'phase_rate' ):
            calculate_phase_rate = 1
        if ( ( list_variable[j] == 'nb_seconds_since_start' ) | (calculate_phase_rate == 1)):
            calculate_nb_seconds_since_start = 1

        if ( list_variable[j] == 'latitude' ):
            calculate_latitude = 1
        if ( list_variable[j] == 'altitude' ):
            calculate_altitude = 1
        if ( list_variable[j] == 'sma' ):
            calculate_sma = 1
            calculate_sma_ave = 1

        if ( list_variable[j] == 'inclination' ):
            calculate_inclination = 1
        if ( list_variable[j] == 'eccentricity' ):
            calculate_eccentricity = 1
        if ( list_variable[j] == 'solar_zenith' ):
            calculate_solar_zenith = 1

        if (list_variable[j] == 'phase_angle' ):
            calculate_phase_angle = 1
        if ( (list_variable[j] == 'true_anomaly' ) | (calculate_phase_angle == 1) ):
            calculate_true_anomaly = 1
        if ( list_variable[j] == 'raan' ):
            calculate_raan = 1
        if ( list_variable[j] == 'beta' ):
            calculate_beta = 1

        if ( list_variable[j] == 'ap' ):
            calculate_ap = 1
        if ( list_variable[j] == 'f107' ):
            calculate_f107 = 1
        if ( list_variable[j] == 'power' ):
            calculate_power = 1
        if ( ( list_variable[j] == 'argument_perigee' ) | (calculate_phase_angle == 1) ):
            calculate_argument_perigee = 1
            calculate_argument_perigee_ave = 1
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
    if ((calculate_position_ecef) | (calculate_velocity_ecef)):
        file_ecef =  open('/'.join(filename.split('/')[:-1]) + "/ECEF_" + filename.split('/')[-1])
        read_file_ecef = file_ecef.readlines() # !!!!! first date of given output is not necesarrily the same as for the other variables 
        nb_lines_header_ecef = 11
    if ((calculate_position_tle) | (calculate_velocity_tle)):
        file_tle =  open('/'.join(filename.split('/')[:-1]) + "/TLE_" + filename.split('/')[-1])
        read_file_tle = file_tle.readlines() # !!!!! first date of given output is not necesarrily the same as for the other variables 
        nb_lines_header_tle = 12
        n_tle = len(read_file_tle) - nb_lines_header_tle
        position_tle = np.zeros([n_tle, 3]) # if TLe were used ot initilize position, this includes the postinon form the tle epoch to the constellation epoch. this is reported by SpOCK in the output file TLE_
        velocity_tle = np.zeros([n_tle, 3])
        date_tle = []
    if ((calculate_ap == 1) | (calculate_f107 == 1)):
        n_given_output = len(read_file_given_output)
        for i in range(n_given_output):
            date_given_output.append(read_file_given_output[i].split()[0])
            if (calculate_ap == 1):
                ap.append(read_file_given_output[i].split()[3]) # !!!!! first date of given output is not necesarrily the same as for the other variables
            if (calculate_f107 == 1):
                f107.append(read_file_given_output[i].split()[1]) # !!!!! first date of given output is not necesarrily the same as for the other variables

    if ((calculate_position_tle) | (calculate_velocity_tle)):
        for i in range(n_tle):
            date_tle.append(read_file_tle[i+nb_lines_header_tle].split()[0] + " " + read_file_tle[i+nb_lines_header_tle].split()[1])
            if (calculate_position_tle == 1):
                position_tle[i,0] = np.float(read_file_tle[i+nb_lines_header_tle].split()[2])
                position_tle[i,1] = np.float(read_file_tle[i+nb_lines_header_tle].split()[3])
                position_tle[i,2] = np.float(read_file_tle[i+nb_lines_header_tle].split()[4])
            if (calculate_velocity_tle == 1):
                velocity_tle[i,0] = np.float(read_file_tle[i+nb_lines_header_tle].split()[5])
                velocity_tle[i,1] = np.float(read_file_tle[i+nb_lines_header_tle].split()[6])
                velocity_tle[i,2] = np.float(read_file_tle[i+nb_lines_header_tle].split()[7])

    
    for i in range(n):
        date.append(read_file_to_read[i+nb_lines_header].split()[0] + " " + read_file_to_read[i+nb_lines_header].split()[1])
        date_datetime.append(datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S.%f"))

        # round to nearest second
        date_here = date[-1]
        if  ( (int)(date_here.split('.')[1][0]) >= 5 ):
            date_here = date_here.split('.')[0]
            date_here = datetime.strptime(date_here, "%Y/%m/%d %H:%M:%S") + timedelta(seconds = 1)
            date_here = datetime.strftime(date_here, "%Y/%m/%d %H:%M:%S")
            date_here = date_here.split('.')[0]
        else:
            date_here = date_here.split('.')[0]
        date_round_sec.append(date_here)
        date_datetime_round_sec.append(datetime.strptime(date_here, "%Y/%m/%d %H:%M:%S"))

        nb_seconds_since_start[i] = (datetime.strptime(date[-1], "%Y/%m/%d %H:%M:%S.%f") - datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")).total_seconds()

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
        if (calculate_position_ecef == 1):
            position_ecef[i,0] = np.float(read_file_ecef[i+nb_lines_header_ecef].split()[2])
            position_ecef[i,1] = np.float(read_file_ecef[i+nb_lines_header_ecef].split()[3])
            position_ecef[i,2] = np.float(read_file_ecef[i+nb_lines_header_ecef].split()[4])
        if (calculate_velocity_ecef == 1):
            velocity_ecef[i,0] = np.float(read_file_ecef[i+nb_lines_header_ecef].split()[5])
            velocity_ecef[i,1] = np.float(read_file_ecef[i+nb_lines_header_ecef].split()[6])
            velocity_ecef[i,2] = np.float(read_file_ecef[i+nb_lines_header_ecef].split()[7])

        if (calculate_argument_perigee == 1):
            argument_perigee[i] = np.float(read_file_to_read[i+nb_lines_header].split()[16])
        if (calculate_phase_angle == 1):
            phase_angle[i] = np.float(np.mod(true_anomaly[i] + argument_perigee[i], 360))
            #phase_angle[i] = np.float(read_file_to_read[i+nb_lines_header].split()[38])
        if (calculate_phase_rate == 1):
            if i >= 1:
                phase_rate[i] = ( ( true_anomaly[i] + argument_perigee[i] ) - ( true_anomaly[i-1] + argument_perigee[i-1] ) ) / ( nb_seconds_since_start[i] - nb_seconds_since_start[i-1] )
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
        if (calculate_acceleration_lvlh_drag_mag == 1):
            acceleration_lvlh_drag_mag[i] = np.linalg.norm(acceleration_lvlh_drag[i,:])
        if (calculate_acceleration_eci_drag == 1):
            acceleration_eci_drag[i,0] = np.float(read_file_to_read[i+nb_lines_header].split()[31])
            acceleration_eci_drag[i,1] = np.float(read_file_to_read[i+nb_lines_header].split()[32])
            acceleration_eci_drag[i,2] = np.float(read_file_to_read[i+nb_lines_header].split()[33])
        if (calculate_beta == 1):
            beta[i] = np.float(read_file_to_read[i+nb_lines_header].split()[37])
        if (calculate_radius_perigee == 1):
            radius_perigee[i] = sma[i]*(1 - eccentricity[i])
        if (calculate_radius_apogee == 1):
            radius_apogee[i] = sma[i]*(1 + eccentricity[i])
        if (('ORB' in read_file_to_read[i+nb_lines_header]) & (i != 0)): # new orbit starts (ignore first line of output file)
            if (calculate_argument_perigee_ave == 1):
                argument_perigee_ave[previous_orb:i+1] = np.float(read_file_to_read[i+nb_lines_header].split()[39]) # 38 is phase angle
            if (calculate_sma_ave == 1):
                sma_ave[previous_orb:i+1] = np.float(read_file_to_read[i+nb_lines_header].split()[40])
            orb_number_temp = (int)(read_file_to_read[i+nb_lines_header].split('ORB')[1].split()[0])
            orb_number.append([i, orb_number_temp]) #[index when new orbit, orbit number]
            previous_orb = i
            # to finish up: from the last orbit until the last time step of the simulation
            argument_perigee_ave[previous_orb:] = np.float(read_file_to_read[i+nb_lines_header].split()[39]) # same value as previous orbit average value
            sma_ave[previous_orb:] = np.float(read_file_to_read[i+nb_lines_header].split()[40]) # same value as previous orbit average value
        if (calculate_solar_zenith == 1):
            solar_zenith[i] = np.float(read_file_to_read[i+nb_lines_header].split()[42])

    variables = []
    order_variables = []
    variables.append(date)
    order_variables.append("date | " + str(len(order_variables)))
    variables.append(date_datetime)
    order_variables.append("date_datetime | " + str(len(order_variables)))

    variables.append(date_round_sec)
    order_variables.append("date_round_sec | " + str(len(order_variables)))
    variables.append(date_datetime_round_sec)
    order_variables.append("date_datetime_round_sec | " + str(len(order_variables)))

    variables.append(orb_number)
    order_variables.append("orb_number | " + str(len(order_variables)))

    if (calculate_position == 1):
        variables.append(position)
        order_variables.append("position | " + str(len(order_variables)))
    if (calculate_position_ecef == 1):
        variables.append(position_ecef)
        order_variables.append("position_ecef | " + str(len(order_variables)))
    if (calculate_velocity_ecef == 1):
        variables.append(velocity_ecef)
        order_variables.append("velocity_ecef | " + str(len(order_variables)))

    if ( (calculate_position_tle == 1) | (calculate_velocity_tle == 1)):
        variables.append(date_tle)
        order_variables.append("date_tle | " + str(len(order_variables)))
    if (calculate_position_tle == 1):
        variables.append(position_tle)
        order_variables.append("position_tle | " + str(len(order_variables)))
    if (calculate_velocity_tle == 1):
        variables.append(velocity_tle)
        order_variables.append("velocity_tle | " + str(len(order_variables)))

    if (calculate_radius == 1):
        variables.append(radius)
        order_variables.append("radius | " + str(len(order_variables)))
    if (calculate_velocity == 1):
        variables.append(velocity)
        order_variables.append("velocity | " + str(len(order_variables)))
    if (calculate_speed == 1):
        variables.append(speed)
        order_variables.append("speed  | " + str(len(order_variables)))
    variables.append(nb_seconds_since_start)
    order_variables.append("nb_seconds_since_start | " + str(len(order_variables)))

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
    if (calculate_sma_ave == 1):
        variables.append(sma_ave)
        order_variables.append("sma_ave | " + str(len(order_variables)))

    if (calculate_inclination == 1):
        variables.append(inclination)
        order_variables.append("inclination | " + str(len(order_variables)))
    if (calculate_eccentricity == 1):
        variables.append(eccentricity)
        order_variables.append("eccentricity | " + str(len(order_variables)))
    if (calculate_solar_zenith == 1):
        variables.append(solar_zenith)
        order_variables.append("solar_zenith | " + str(len(order_variables)))

    if (calculate_true_anomaly == 1):
        variables.append(true_anomaly)
        order_variables.append("true_anomaly | " + str(len(order_variables)))
    if (calculate_phase_angle == 1):
        variables.append(phase_angle)
        order_variables.append("phase_angle | " + str(len(order_variables)))
    if (calculate_phase_rate == 1):
        variables.append(phase_rate)
        order_variables.append("phase_rate | " + str(len(order_variables)))

    if (calculate_raan == 1):
        variables.append(raan)
        order_variables.append("raan | " + str(len(order_variables)))
    if (calculate_beta == 1):
        variables.append(beta)
        order_variables.append("beta | " + str(len(order_variables)))
    if (calculate_ap == 1): # !!!!! first date of given output is not necesarrily the same as for the other variables 
        variables.append(ap)
        order_variables.append("ap | " + str(len(order_variables)))
    if (calculate_f107 == 1): # !!!!! first date of given output is not necesarrily the same as for the other variables 
        variables.append(f107)
        order_variables.append("f107 | " + str(len(order_variables)))
    if ((calculate_f107 == 1) | (calculate_ap == 1)): # !!!!! first date of given output is not necesarrily the same as for the other variables 
        variables.append(date_given_output)
        order_variables.append("date_given_output | " + str(len(order_variables)))

    if (calculate_power == 1): # !!!!! first date of given output is not necesarrily the same as for the other variables 
        variables.append(power)
        order_variables.append("power | " + str(len(order_variables)))

    if (calculate_argument_perigee == 1):
        variables.append(argument_perigee)
        order_variables.append("argument_perigee | " + str(len(order_variables)))

    if (calculate_argument_perigee_ave == 1):
        variables.append(argument_perigee_ave)
        order_variables.append("argument_perigee_ave | " + str(len(order_variables)))
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
    if (calculate_acceleration_lvlh_drag_mag == 1):
        variables.append(acceleration_lvlh_drag_mag)
        order_variables.append("acceleration_lvlh_drag_mag | " + str(len(order_variables)))

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
        for i in range(n):
            if calculate_density == 1:
                density[i] = np.float(density_file_read[i+nb_lines_header].split()[2]) / 10**9 # in kg/m^3
            if calculate_temperature == 1:
                temperature[i] = np.float(density_file_read[i+nb_lines_header].split()[3]) # in K
            if calculate_cd == 1:
                cd[i] = np.float(density_file_read[i+nb_lines_header].split()[4]) # in K
            if calculate_tot_area_drag == 1:
                tot_area_drag[i] = np.float(density_file_read[i+nb_lines_header].split()[5]) * 1e10  # in cm^2 (output in km^2 by SpOCK)

#         if calculate_density == 1:
#             density[0] = density[1] # since the ouptut file misses the first time step of the similution, we just make this density equal to the one at the second time step
#         if calculate_temperature ==1:
#             temperature[0] = temperature[1] # since the ouptut file misses the first time step of the similution, we just make this temperature equal to the one at the second time step
#         if calculate_cd ==1:
#             cd[0] = cd[1] # since the ouptut file misses the first time step of the similution, we just make this cd equal to the one at the second time step
#         if calculate_tot_area_drag ==1:
#             tot_area_drag[0] = tot_area_drag[1] # since the ouptut file misses the first time step of the similution, we just make this tot_area_drag equal to the one at the second time step

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

