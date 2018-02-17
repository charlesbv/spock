# This script creates a main input file for SpOCK
# ASSUMPTIONS:
# - need to be in spokc/srcPython to run this script

import os

def spock_main_input( 
    main_input_filename,
    # for TIME section
    date_start,
    date_end,
    dt,
    # for SPACECRAFT section
    nb_sc,
    gps_tle,
    mass,
    geometry_filename,
    # for ORBIT section
    tle_filename,
    # for FORCES section
    gravity_order,
    forces,
    density_mode,
    # for OUTPUT section
    name_output,
    dt_output,
    # for ATTITUDE section
    attitude,
    # for GROUND_STATIONS section
    filename_ground_stations,
    # for SPICE section
    spice_path,
    # for DENSITY_MOD section
    rho_mod
    ):

    # ################################################################
    # # PARAMETERS TO SET
    # run_dir = "run.cygnss" # no path, just name
    # main_input_filename = "test_main_input.txt"
    # # for TIME section
    # date_start = "2016-08-25T12:00:00" # use this exact same format: YYYY-MM-DDTHH:MM:SS
    # date_end = "2016-08-26T12:00:00" # use this exact same format: YYYY-MM-DDTHH:MM:SS
    # dt = 60
    # # for ORBIT section
    # tle_filename = "ex_tle.txt" # no path
    # # for FORCES section
    # gravity_order = 4
    # forces = "drag"
    # # for OUTPUT section
    # name_output = "out" # no path
    # dt_output = 60

    ################################################################
    # ALGORITHM
    main_input_file  = open( main_input_filename, "w")

    # TIME section
    print >> main_input_file, "#TIME\n" + date_start + "\n" + date_end + "\n" + str(dt) + "\n"

    # SPACECRAFT section
    print >> main_input_file, "#SPACECRAFT\n" + str(nb_sc) + "\n" + gps_tle + "\n" + str(mass) + "\n-1\n" + geometry_filename  + "\n" 

    # ORBIT section
    if type(tle_filename) == str:
        print >> main_input_file, "#ORBIT\ntle\n" + tle_filename + "\n"
    else:
        print >> main_input_file, "#ORBIT\n" + tle_filename[0] 
        if tle_filename[0] == 'collision':
            print >> main_input_file, tle_filename[1]
        else:
            isc = 0
            while isc < nb_sc:
                print >> main_input_file, tle_filename[isc + 1] 
                isc = isc + 1
        print >> main_input_file, ""
    # ATTITUDE section
    print >> main_input_file, "#ATTITUDE\n" + attitude + "\n"

    # FORCES section
    if density_mode == 'dynamic': # by default if dynamic is chosen then F10.7 and Ap are downloaded from swpc -> I then added the option density_mode == swpc but kept this one 'dynamic' because other scripts were using the option dynamic so I didn't want to change all these scripts
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\ndynamic\nswpc\n100\n15\n"
    elif density_mode == 'swpc': # by default if dynamic is chosen then F10.7 and Ap are downloaded from swpc
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\ndynamic\nswpc\n100\n15\n"
    elif density_mode == 'omniweb': # by default if dynamic is chosen then F10.7 and Ap are downloaded from swpc
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\ndynamic\nomniweb\n100\n15\n"
    elif density_mode == 'static':
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\nstatic\n100\n100\n15\n"
    elif ( ( type(density_mode) == list ) & ( density_mode[0] == 'static' ) ):
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\nstatic\n" + density_mode[1] + "\n" + density_mode[2] + "\n" + density_mode[3] + "\n"  # 1 is f107, 2 is f107A, 3 is ap. I kept the option 'static' (right above) because it's from an old version that many scripts were using so I didn't to change all these scripts
    elif ( ( type(density_mode) == list ) & ( density_mode[0] == 'swpc_mod' ) ): 
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\ndynamic\nswpc_mod " + str(density_mode[1]) + "\n"
    else: # files for F10.7 and Ap input by user (first line is F10.7 filename, second line is Ap filename)
        print >> main_input_file, "#FORCES\n" + str(gravity_order) + "\n" + forces + "\ndynamic\n" + density_mode[0] + "\n" + density_mode[1] + "\n"

    
    #GROUND_STATIONS section
    if filename_ground_stations != '0': # if '0' then no section GROUND_STATIONS will be written in the main input file 
        print >> main_input_file, "#GROUND_STATIONS\n" + filename_ground_stations + "\n" 
    # OUTPUT section
    print >> main_input_file, "#OUTPUT\n" + name_output + "\n" + str(dt_output) + "\n"


    #SPICE section
    print >> main_input_file, "#SPICE\n" + spice_path + "\n"

    #DENSITY_MOD section
    print >> main_input_file, "#DENSITY_MOD" 
    # the weird thing I do below is because I didn't want to change all thes ciprt that were calling spock_main_input
    if type(rho_mod) != list:
        print >> main_input_file,  str(rho_mod) + "\n"
    else:
        nb_other_stuff = len(rho_mod)
        for iother in range(nb_other_stuff): # the first one needs to be the rho_mod factor, the rest can be any option section
            print >> main_input_file,  str(rho_mod[iother]) 
        print >> main_input_file,  ""

    main_input_file.close()

    return 0
