import glob
import numpy as np
from datetime import datetime, timedelta
import time
import os
import fileinput
def read_input_file(filename):
    #filename = "minxss_1_0_omniweb.txt"

    filename_no_path = filename.split('/')[-1]
    filename_no_path_no_extension = filename_no_path.replace(".txt","")
    file_input = open(filename, "r")
    read_file_input = file_input.readlines()

     # for i in range(len(filename.split('/')) - 3):
     #     if (i > 0):
     #         run_dir = run_dir + "/" + filename.split('/')[i]
    run_dir = '/'.join(filename.split('/')[:-1]) + "/"

    if run_dir == "/":
        run_dir = "./"

    # Find all sections
    name_section = '#TIME'
    found_section = 0
    line_start_section_time = 0
    while ( ( line_start_section_time < len(read_file_input) ) & ( ( name_section in read_file_input[line_start_section_time] ) == 0 ) ):
        line_start_section_time = line_start_section_time + 1

    name_section = '#SPACECRAFT'
    found_section = 0
    line_start_section_spacecraft = 0
    while ( ( line_start_section_spacecraft < len(read_file_input) ) & ( ( name_section in read_file_input[line_start_section_spacecraft] ) == 0 ) ):
        line_start_section_spacecraft = line_start_section_spacecraft + 1

    name_section = '#ORBIT'
    found_section = 0
    line_start_section_orbit = 0
    while ( ( line_start_section_orbit < len(read_file_input) ) & ( ( name_section in read_file_input[line_start_section_orbit] ) == 0 ) ):
        line_start_section_orbit = line_start_section_orbit + 1

    name_section = '#FORCES'
    found_section = 0
    line_start_section_forces = 0
    while ( ( line_start_section_forces < len(read_file_input) ) & ( ( name_section in read_file_input[line_start_section_forces] ) == 0 ) ):
        line_start_section_forces = line_start_section_forces + 1

    name_section = '#KALMAN'
    found_section = 0
    line_start_section_kalman = 0
    while ( ( line_start_section_kalman < len(read_file_input) ) & ( ( name_section in read_file_input[line_start_section_kalman] ) == 0 ) ):
        line_start_section_kalman = line_start_section_kalman + 1
        if ( line_start_section_kalman == len(read_file_input)):
             break

    name_section = '#OUTPUT'
    found_section = 0
    line_start_section_output = 0
    while ( ( line_start_section_output < len(read_file_input) ) & ( ( name_section in read_file_input[line_start_section_output] ) == 0 ) ):
        line_start_section_output = line_start_section_output + 1

    line_out = read_file_input[1+line_start_section_output].split()[0].replace("\n", "").replace("\r","")
    name_section = '#STORMS'
    found_section = 0
    line_start_section_storms = 0
    while ( ( line_start_section_storms < len(read_file_input) -1 ) & ( ( name_section in read_file_input[line_start_section_storms] ) == 0 ) ):
        line_start_section_storms = line_start_section_storms + 1

    name_section = '##ENSEMBLES_COE'
    found_section = 0
    line_start_section_ensembles_coe = 0
    while ( ( line_start_section_ensembles_coe < len(read_file_input) -1 ) & ( ( name_section in read_file_input[line_start_section_ensembles_coe] ) == 0 ) ):
        line_start_section_ensembles_coe = line_start_section_ensembles_coe + 1

    name_section = '#ATTITUDE'
    found_section = 0
    line_start_section_attitude = 0
    while ( ( line_start_section_attitude < len(read_file_input) -1 ) & ( ( name_section in read_file_input[line_start_section_attitude] ) == 0 ) ):
        line_start_section_attitude = line_start_section_attitude + 1

    name_section = '##ENSEMBLES_ATTITUDE'
    found_section = 0
    line_start_section_ensembles_attitude = 0
    while ( ( line_start_section_ensembles_attitude < len(read_file_input) -1 ) & ( ( name_section in read_file_input[line_start_section_ensembles_attitude] ) == 0 ) ):
        line_start_section_ensembles_attitude = line_start_section_ensembles_attitude + 1

    name_section = '#GROUND_STATIONS'
    found_section = 0
    line_start_section_ground_stations = 0
    while ( ( line_start_section_ground_stations < len(read_file_input) -1 ) & ( ( name_section in read_file_input[line_start_section_ground_stations] ) == 0 ) ):
        line_start_section_ground_stations = line_start_section_ground_stations + 1


    name_section = '#OUTPUT_ENSEMBLES'
    found_section = 0
    line_start_section_output_ensembles = 0
    while ( ( line_start_section_output_ensembles < len(read_file_input) -1 ) & ( ( name_section in read_file_input[line_start_section_output_ensembles] ) == 0 ) ):
        line_start_section_output_ensembles = line_start_section_output_ensembles + 1


    # NUMBER OF SATELLITES
    nb_sc =  (int)(read_file_input[1 + line_start_section_spacecraft].split()[0]) 


    # OUTPUT FILE NAMES
    output_file_name_list = []
    output_file_path_list = []
    collision_filename = ''
    if '/' in line_out:
        path_run_dir = '/'.join(line_out.split('/')[:-1]) + "/"
    else:
        path_run_dir = './'
    if path_run_dir[0] == '~':
            path_run_dir=  os.path.expanduser('~') + path_run_dir[1:]
    common_name = line_out.split('/')[-1]
    if common_name != "now":
        if common_name != "out":
            common_name = line_out.split('/')[-1]
        else:
            common_name = filename_no_path_no_extension

        for i in range(nb_sc):
            output_file_path_list.append( path_run_dir + common_name + '/' + common_name + str(i+1) + "/")
            output_file_name_list.append(common_name + str(i+1) + ".txt")
            collision_filename =  path_run_dir  + common_name + "_collision.txt"

    else:
        common_name = max(glob.iglob(path_run_dir + '*/'), key=os.path.getctime).split('/')[-2]
        for i in range(nb_sc):
            output_file_path_list.append( path_run_dir + common_name + '/' + common_name + str(i+1) + "/")
            output_file_name_list.append(common_name + str(i+1) + ".txt")
            collision_filename =  path_run_dir  + common_name + "_collision.txt"





    # TIME
    ## START DATE
    start_date_now_n = 0
    while True:
        try:
            date_start = read_file_input[1+line_start_section_time].split()[0]
            date_start = datetime.strptime(date_start, "%Y-%m-%dT%H:%M")
            break
        except ValueError:
            try: 
                date_start = read_file_input[1+line_start_section_time].split()[0]
                date_start = datetime.strptime(date_start, "%Y-%m-%dT%H:%M:%S")
                break
            except ValueError:
                try: 
                    date_start = read_file_input[1+line_start_section_time].split()[0]
                    date_start = datetime.strptime(date_start, "%Y-%m-%dT%H:%M:%S.%f")
                    break
                except ValueError:
                    try: 
                        date_start = read_file_input[1+line_start_section_time].split()[0] + " " + read_file_input[1+line_start_section_time].split()[1]
                        date_start = datetime.strptime(date_start, "%Y-%m-%d %H:%M")
                        break
                    except ValueError:
                        try: 
                            date_start = read_file_input[1+line_start_section_time].split()[0] + " " + read_file_input[1+line_start_section_time].split()[1]
                            date_start = datetime.strptime(date_start, "%Y-%m-%d %H:%M:%S")
                            break
                        except ValueError:
                            try: 
                                if ("now" != read_file_input[1+line_start_section_time].split()[0]):
                                    raise ValueError
                                n_from_now = np.float(read_file_input[1+line_start_section_time].split()[1])
                                start_date_now_n = 1
                                break
                            except ValueError:
                                print "The start date shold have one of these formats: %Y-%m-%d %H:%M:%S or %Y-%m-%d %H:%M or %Y-%m-%dT%H:%M:%S or %Y-%m-%dT%H:%M or 'now n' (n being a number). The program will stop."
                                raise Exception


    if (start_date_now_n == 1): # if the user inputs the start date with 'now n' (n being the number of (decimal) hours to run the propagation for)
    #     if line_out.replace("\n", "") != "now":
    #         if line_out.replace("\n", "") != "out":
    #             one_of_the_filename_output = "spock/" + line_out.replace("\n","") + "/" + line_out.replace("\n","") + str(1) + "/" + line_out.replace("\n","")  + str(1) + ".txt"
    #         else:
    #             one_of_the_filename_output = "spock/" + filename_no_path_no_extension + "/" + filename_no_path_no_extension + str(1) + "/" + filename_no_path
    #     else:
    #         newest_run = max(glob.iglob('spock/*/'), key=os.path.getctime)
    #         one_of_the_filename_output =  newest_run + newest_run.split('/')[-2] +  "1/" + newest_run.split('/')[-2] + "1.txt"

    #    one_of_the_file_output = open(run_dir + one_of_the_filename_output, "r")
        one_of_the_file_output = open(output_file_path_list[0] + output_file_name_list[0], "r")
        one_of_the_file_output_read = one_of_the_file_output.readlines()
        date_start = one_of_the_file_output_read[10].split()[0] + " " + one_of_the_file_output_read[10].split()[1]
        date_start = datetime.strptime(date_start, "%Y/%m/%d %H:%M:%S")
        one_of_the_file_output.close()

        date_stop = date_start + timedelta(hours = n_from_now)
    ## END DATE
    else:
        while True:
            try:
                date_stop = read_file_input[2+line_start_section_time].split()[0]
                date_stop = datetime.strptime(date_stop, "%Y-%m-%dT%H:%M")
                break
            except ValueError:
                try: 
                    date_stop = read_file_input[2+line_start_section_time].split()[0]
                    date_stop = datetime.strptime(date_stop, "%Y-%m-%dT%H:%M:%S")
                    break
                except ValueError:
                    try: 
                        date_stop = read_file_input[2+line_start_section_time].split()[0]
                        date_stop = datetime.strptime(date_stop, "%Y-%m-%dT%H:%M:%S.%f")
                        break
                    except ValueError:
                        try: 
                            date_stop = read_file_input[2+line_start_section_time].split()[0] + " " + read_file_input[2+line_start_section_time].split()[1]
                            date_stop = datetime.strptime(date_stop, "%Y-%m-%d %H:%M")
                            break
                        except ValueError:
                            try: 
                                date_stop = read_file_input[2+line_start_section_time].split()[0] + " " + read_file_input[2+line_start_section_time].split()[1]
                                date_stop = datetime.strptime(date_stop, "%Y-%m-%d %H:%M:%S")
                                break
                            except ValueError:
                                print "The end date should have one of these formats: %Y-%m-%d %H:%M:%S or %Y-%m-%d %H:%M or %Y-%m-%dT%H:%M:%S or %Y-%m-%dT%H:%M. The program will stop."
                                raise Exception

    ## TIME STEP AND NUMBER OF STEPS IN THE SIMULATION
    delta_date = date_stop - date_start
    nb_seconds_simu = delta_date.days * 24 * 3600 + delta_date.seconds + delta_date.microseconds/10**6.
    dt = np.float(read_file_input[3+line_start_section_time].split()[0])
    dt_output = np.float(read_file_input[2+line_start_section_output].split()[0])
    if dt_output < dt: # the output time step is set to the integration time step if the user put a smaller output time step by mistake
        dt_output = dt
    nb_steps = (int)(nb_seconds_simu/dt_output) + 1
    if ( np.mod(nb_seconds_simu, dt_output) != 0 ):
        nb_steps = nb_steps + 1


    # NUMBER OF GPS
    gps_name = []
    nb_gps = 0
    filename_gps_tle = ''
    if (read_file_input[2+line_start_section_spacecraft].split()[0] != '0'):
        filename_gps_tle = read_file_input[2+line_start_section_spacecraft].split()[0]
        file_gps_tle = open(run_dir + filename_gps_tle, "r")
        read_file_gps_tle = file_gps_tle.readlines()
        for i in range(len(read_file_gps_tle)):
            if read_file_gps_tle[i].split()[0] == '1':
                gps_name.append( read_file_gps_tle[i-1].split()[1] )
                nb_gps = nb_gps + 1



    #ORBIT
    if (read_file_input[line_start_section_orbit + 1].split()[0].replace("\n","") == "collision"):
        filename_collision = read_file_input[line_start_section_orbit + 2].split()[0]
        file_collision = open(run_dir + filename_collision)
        read_file_collision = file_collision.readlines()
        find_nb_ens = 0
        while (read_file_collision[find_nb_ens].split(' ')[0].replace("\n","") != "#NB_ENSEMBLES_COLLISION" ):
            find_nb_ens = find_nb_ens + 1
        nb_ensembles_coe = (int)( read_file_collision[find_nb_ens + 1].split()[0] )

        file_collision.close()

    #STORMS
    storm_name = []
    if (line_start_section_storms == len(read_file_input) - 1):
        nb_storms = 0
    else:
        nb_storms = (int)(len( read_file_input[1+line_start_section_storms].split(',')))
    for i in range(nb_storms):
        storm_name.append( (read_file_input[1+line_start_section_storms].split(',')[i]).replace("\n","").replace(" ","") )

    #FORCES
    nb_ensembles_density = 0
    compute_drag = 0
    if ( "drag" in read_file_input[line_start_section_forces+2] ):
        compute_drag = 1
        if ( read_file_input[line_start_section_forces+3].split()[0] == "dynamic" ):
            if ( read_file_input[line_start_section_forces+4].split()[0] == "swpc" ):
                if ( len(read_file_input[line_start_section_forces+4].split()) >= 2 ):
                    nb_ensembles_density = (int) (read_file_input[line_start_section_forces+4].split()[1])

    # MASS
    mass =  np.float(read_file_input[3+line_start_section_spacecraft].split()[0]) 

    # SOLAR CELL EFFICIENCY
    solar_cell_eff = np.float(read_file_input[4+line_start_section_spacecraft].split()[0]) 

    # GEOMETRY FILE NAME
    geometry_filename_temp = read_file_input[5+line_start_section_spacecraft].split()[0]
    if geometry_filename_temp != "ballistic_coefficient":
        geometry_filename = ""
        # for i in range(len(filename.split('/')) - 2):
        #     if (i > 0):
        #         geometry_filename = geometry_filename + "/" + filename.split('/')[i]
        #geometry_filename = geometry_filename + '/geometry/' + geometry_filename_temp
        geometry_filename = geometry_filename_temp
        geometry_file = open(run_dir + geometry_filename)
        read_geometry_file = geometry_file.readlines()
        n_header_geo = 0
        while (read_geometry_file[n_header_geo].split()[0] != '#ENDOFHEADER'):
            n_header_geo = n_header_geo + 1
        n_header_geo = n_header_geo + 1
        ## Figure out if there is a section #NB_ENSEMBLES_CD
        section_ensemble_cd = 0
        nb_surfaces = 0
        for iline in range(len(read_geometry_file) - n_header_geo):
            if '#NB_ENSEMBLES_CD' in read_geometry_file[iline + n_header_geo]:
                section_ensemble_cd = 1
                nb_ensembles_cd = (int)(read_geometry_file[iline + n_header_geo + 1].split()[0] )
                nb_surfaces = -1
        if nb_surfaces == 0:
            nb_ensembles_cd = 0
        cdd = []
        area = []
        for iline in range(len(read_geometry_file) - n_header_geo):
            if len(read_geometry_file[iline + n_header_geo].split()) > 0:
                if read_geometry_file[iline + n_header_geo][0] == '#':
                    cd_temp = np.float( read_geometry_file[iline + n_header_geo+4].split()[0] )
                    area_temp = np.float( read_geometry_file[iline + n_header_geo+2].split()[0] ) / (100000**2) # cm^2 to km^2
                    cdd.append(cd_temp)
                    area.append(area_temp) # area in km^2
                    nb_surfaces = nb_surfaces + 1
        geometry_file.close()

    else:
        nb_surfaces = 1
        nb_ensembles_cd = 0
    # # TYPE OF ORBIT INITIALIZATION
    # type_orbit_ini = read_file_input[14].split()[0]

    # # ORBIT INITIALIZATION FOR EACH SPACECRAFT
    # ## TO DO

    # # NUMBER OF ENSEMBLES ON THE ORBITAL ELEMENTS
    if ( read_file_input[line_start_section_orbit + 1].split()[0].replace("\n","") != "collision" ):
        if ( ( line_start_section_ensembles_coe == len(read_file_input) - 1 ) ):
            nb_ensembles_coe = 0
        else:
            nb_ensembles_coe = (int)(read_file_input[1+line_start_section_ensembles_coe].split()[0])

    # # STANDARD DEVIATION FOR THE ENSEMBLES ON THE ORBITAL ELEMENTS FOR EACH MAIN SPACECRAFT
    # ## TO DO

    # # ATTITUDE TYPE
    # attitude_type = read_file_input[22].split()[0]

    # NUMBER OF ENSEMBLES ON THE ATTITUDE
    if (line_start_section_ensembles_attitude == len(read_file_input) - 1):
        nb_ensembles_attitude = 0
    else:
        nb_ensembles_attitude = (int)(read_file_input[1+line_start_section_ensembles_attitude].split()[0])

    # # ENSEMBLES TO OUTPUT
    ensembles_to_output = []
    if line_start_section_output_ensembles < len(read_file_input) - 1:
        for i in range(len(read_file_input[line_start_section_output_ensembles+1].replace(" ","").split(','))):
            if (i < len(read_file_input[line_start_section_output_ensembles+1].replace(" ","").split(',')) - 1):
                ensembles_to_output.append(read_file_input[line_start_section_output_ensembles+1].replace(" ","").split(',')[i].replace("\n",""))
            else:
                ensembles_to_output.append(read_file_input[line_start_section_output_ensembles+1].replace(" ","").split(',')[i].split('//')[0].replace("\n",""))

    # KALMAN
    filename_kalman_init = ""
    filename_kalman_meas_in = ""
    filename_kalman_meas_out = [] # kalman measurement output file: the measruement output here is the same unit as the state. So it's the conversion of the measuremnt input file (filename_kalman_meas_in) to the unit of the state. Ex: if the measurements in filename_kalman_meas_in are r/v in ECEF but the state r/v in ECI then the file filename_kalman_meas_out is the measurements r/v in ECI. 
    if line_start_section_kalman < len(read_file_input) - 1:
        filename_kalman_init = read_file_input[line_start_section_kalman+1].split()[0]
        file_kalman_init = open(filename_kalman_init)
        filename_kalman_meas_in = file_kalman_init.readlines()[0].split()[0]
        for isc in range(nb_sc):
            filename_kalman_meas_out.append( "meas_converted_kalman_" + output_file_name_list[isc])
        file_kalman_init.close()

    # # COVERAGE GROUND STATIONS
    filename_ground_station = ""
    if line_start_section_ground_stations < len(read_file_input) - 1:
        filename_ground_station = read_file_input[line_start_section_ground_stations].replace(" ","").split()[0]

    order_variables = []
    variables = []
    variables.append(date_start)
    order_variables.append("date_start | " + str(len(order_variables)))
    variables.append(date_stop)
    order_variables.append("date_stop | " + str(len(order_variables)))
    variables.append(dt)
    order_variables.append("dt | " + str(len(order_variables)))

    variables.append(dt_output)
    order_variables.append("dt_output | " + str(len(order_variables)))
    variables.append(nb_steps)
    order_variables.append("nb_steps | " + str(len(order_variables)))
    variables.append(nb_sc)
    order_variables.append("nb_sc | " + str(len(order_variables)))
    variables.append(gps_name)
    order_variables.append("gps_name | " + str(len(order_variables)))
    variables.append(output_file_path_list)
    order_variables.append("output_file_path_list | " + str(len(order_variables)))
    variables.append(output_file_name_list)
    order_variables.append("output_file_name_list | " + str(len(order_variables)))
    variables.append(nb_surfaces)
    order_variables.append("nb_surfaces | " + str(len(order_variables)))
    variables.append(nb_storms)
    order_variables.append("nb_storms | " + str(len(order_variables)))
    variables.append(storm_name)
    order_variables.append("storm_name | " + str(len(order_variables)))
    variables.append(nb_ensembles_cd)
    order_variables.append("nb_ensembles_cd | " + str(len(order_variables)))
    variables.append(nb_ensembles_coe)
    order_variables.append("nb_ensembles_coe | " + str(len(order_variables)))
    variables.append(nb_ensembles_attitude)
    order_variables.append("nb_ensembles_attitude | " + str(len(order_variables)))
    variables.append(ensembles_to_output)
    order_variables.append("ensembles_to_output | " + str(len(order_variables)))
    variables.append(filename_kalman_init)
    order_variables.append("filename_kalman_init | " + str(len(order_variables)))
    variables.append(filename_kalman_meas_in)
    order_variables.append("filename_kalman_meas_in | " + str(len(order_variables)))
    variables.append(filename_kalman_meas_out)
    order_variables.append("filename_kalman_meas_out | " + str(len(order_variables)))

    variables.append(filename_ground_station)
    order_variables.append("filename_ground_station | " + str(len(order_variables)))
    variables.append(filename_gps_tle)
    order_variables.append("filename_gps_tle | " + str(len(order_variables)))
    variables.append(nb_ensembles_density)
    order_variables.append("nb_ensembles_density | " + str(len(order_variables)))
    variables.append(collision_filename)
    order_variables.append("collision_filename | " + str(len(order_variables)))
    variables.append(compute_drag)
    order_variables.append("compute_drag | " + str(len(order_variables)))
    variables.append(cdd)
    order_variables.append("cdd | " + str(len(order_variables)))
    variables.append(area)
    order_variables.append("area | " + str(len(order_variables)))
    variables.append(mass)
    order_variables.append("mass | " + str(len(order_variables)))
    
    # variables.append(solar_cell_eff)
    # order_variables.append("solar_cell_eff | " + str(len(order_variables)))
    # variables.append(geometry_filename)
    # order_variables.append("geometry_filename | " + str(len(order_variables)))
    # variables.append(type_orbit_ini)
    # order_variables.append("type_orbit_ini | " + str(len(order_variables)))
    # variables.append(attitude_type)
    # order_variables.append("attitude_type | " + str(len(order_variables))) 
    #print output_file_path_list 
    return variables, order_variables

def find_in_read_input_order_variables(order_variables, var_to_find):
    
    index_of_var_to_find = (int)([s for s in order_variables if var_to_find in s][0].split('|')[1])

    return index_of_var_to_find


# print date_start
# print date_stop
# print dt_output
# print nb_steps
# print nb_sc
# print gps_name
# print output_file_path_list
# print output_file_name_list
# print nb_surfaces
# print nb_storms
# print storm_name
# print nb_ensembles_cd
# print nb_ensembles_coe
# print nb_ensembles_attitude
# print ensembles_to_output
# print filename_ground_station
# print filename_gps_tle
# print nb_ensembles_density
# print collision_filename
# print compute_drag 
  
