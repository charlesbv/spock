interv_dur = 11 # in min, duration of the interval of time to look at. It'sa ctually not exactly in minutes, see comment right below for delta_inter
interv_dur = (int)(interv_dur);
inter_dur_sec = interv_dur * 60 # It' actually not exactly in seconds, see comment right above for delta_inter

delta_inter = 60. # inter_dur_sec#!!!!!! was putting 60 before 04-23-2019 # move forward the interval of time interv_dur delta_inter indices (in theory seconds but sometimes there is missing data in the
#netcdf so you imght jump more than delta_inter seconds) 
delta_inter = (int)(delta_inter)

first_and_second_score = []
first_or_second_score = []
first_score = []
second_score = []
seconds_sampling_start = []
seconds_sampling_stop = []
gap_prn_save = []
gap_first_half_prn_save = []
port_ant_for_prn = []

time_first_score_1st_half_wrong = []
duration_first_score_1st_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_score_1st_half_wrong_list_conc = []
time_first_score_1st_half_wrong_or_port = []
duration_first_score_1st_half_wrong_or_port_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_score_1st_half_wrong_or_port_list_conc = []
time_second_score_1st_half_wrong_or_port = []
duration_second_score_1st_half_wrong_or_port_idate = [] # not really a duration in seconds but a numebr of steps
duration_second_score_1st_half_wrong_or_port_list_conc = []

# -> for each interval iinterv_dur, duration_first_score_1st_half_wrong_idate[iinter_dur] is the number of steps for which the PRN predicted by SpOCK to
# be highest at the first step of iinterv_dur is actually not selected by the onboard algorithm
time_second_score_1st_half_wrong = []
duration_second_score_1st_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
time_first_and_second_score_1st_half_wrong = []
duration_first_and_second_score_1st_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
time_first_or_second_score_1st_half_wrong = []
duration_first_or_second_score_1st_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps

duration_second_score_1st_half_wrong_list_conc = []
duration_first_and_second_score_1st_half_wrong_list_conc = []
duration_first_or_second_score_1st_half_wrong_list_conc = []

time_first_or_second_score_1st_half_wrong_or_port = []
duration_first_or_second_score_1st_half_wrong_or_port_idate = [] # not really a duration in seconds but a numebr of steps
time_first_and_second_score_1st_half_wrong_or_port = []
duration_first_and_second_score_1st_half_wrong_or_port_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_or_second_score_1st_half_wrong_or_port_list_conc = []
duration_first_and_second_score_1st_half_wrong_or_port_list_conc = []

time_first_or_second_score_2nd_half_wrong_or_star = []
duration_first_or_second_score_2nd_half_wrong_or_star_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_or_second_score_2nd_half_wrong_or_star_list_conc = []
time_first_and_second_score_2nd_half_wrong_or_star = []
duration_first_and_second_score_2nd_half_wrong_or_star_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_and_second_score_2nd_half_wrong_or_star_list_conc = []

time_netcdf_two_star_2nd_half = []
duration_netcdf_two_star_2nd_half_idate = [] # not really a duration in seconds but a numebr of steps
duration_netcdf_two_star_2nd_half_list_conc = []
time_netcdf_one_star_2nd_half = []
duration_netcdf_one_star_2nd_half_idate = [] # not really a duration in seconds but a numebr of steps
duration_netcdf_one_star_2nd_half_list_conc = []


time_first_score_2nd_half_wrong = []
duration_first_score_2nd_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_score_2nd_half_wrong_list_conc = []
time_first_score_2nd_half_wrong_or_star = []
duration_first_score_2nd_half_wrong_or_star_idate = [] # not really a duration in seconds but a numebr of steps
duration_first_score_2nd_half_wrong_or_star_list_conc = []
time_second_score_2nd_half_wrong_or_star = []
duration_second_score_2nd_half_wrong_or_star_idate = [] # not really a duration in seconds but a numebr of steps
duration_second_score_2nd_half_wrong_or_star_list_conc = []

# -> for each interval iinterv_dur, duration_first_score_2nd_half_wrong_idate[iinter_dur] is the number of steps for which the PRN predicted by SpOCK to
# be highest at the first step of iinterv_dur is actually not selected by the onboard algorithm
time_second_score_2nd_half_wrong = []
duration_second_score_2nd_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
time_first_and_second_score_2nd_half_wrong = []
duration_first_and_second_score_2nd_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
time_first_or_second_score_2nd_half_wrong = []
duration_first_or_second_score_2nd_half_wrong_idate = [] # not really a duration in seconds but a numebr of steps
duration_second_score_2nd_half_wrong_list_conc = []
duration_first_and_second_score_2nd_half_wrong_list_conc = []
duration_first_or_second_score_2nd_half_wrong_list_conc = []

for idate in range(nb_date):#(1,2):#!!!!!nb_date):
    #idate = 0
    ntot = len(gps_spock_all_date[idate])
    gap_prn_save_idate = []
    first_and_second_score_idate = []
    first_or_second_score_idate = []
    first_score_idate = []
    second_score_idate = []
    seconds_sampling_start_idate = []
    seconds_sampling_stop_idate = []
    port_ant_for_prn_idate = []
    gap_first_half_prn_save_idate = []

    time_first_score_1st_half_wrong_idate = []
    duration_first_score_1st_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_score_1st_half_wrong_or_port_idate = []
    duration_first_score_1st_half_wrong_or_port_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_second_score_1st_half_wrong_or_port_idate = []
    duration_second_score_1st_half_wrong_or_port_list_idate = [] # not really a duration in seconds but a numebr of steps
    
    # -> for each interval iinterv_dur, duration_first_score_1st_half_wrong_list[iinter_dur] is the number of steps for which the PRN predicted by SpOCK to
    # be highest at the first step of iinterv_dur is actually not selected by the onboard algorithm
    time_second_score_1st_half_wrong_idate = []
    duration_second_score_1st_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_and_second_score_1st_half_wrong_idate = []
    duration_first_and_second_score_1st_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_or_second_score_1st_half_wrong_idate = []
    duration_first_or_second_score_1st_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps

    duration_first_or_second_score_1st_half_wrong_or_port_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_or_second_score_1st_half_wrong_or_port_idate = []
    time_first_and_second_score_1st_half_wrong_or_port_idate = []
    duration_first_and_second_score_1st_half_wrong_or_port_list_idate = [] # not really a duration in seconds but a numebr of steps

    duration_first_or_second_score_2nd_half_wrong_or_star_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_or_second_score_2nd_half_wrong_or_star_idate = []
    time_first_and_second_score_2nd_half_wrong_or_star_idate = []
    duration_first_and_second_score_2nd_half_wrong_or_star_list_idate = [] # not really a duration in seconds but a numebr of steps

    duration_netcdf_two_star_2nd_half_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_netcdf_two_star_2nd_half_idate = []
    duration_netcdf_one_star_2nd_half_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_netcdf_one_star_2nd_half_idate = []
    
    time_first_score_2nd_half_wrong_idate = []
    duration_first_score_2nd_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_score_2nd_half_wrong_or_star_idate = []
    duration_first_score_2nd_half_wrong_or_star_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_second_score_2nd_half_wrong_or_star_idate = []
    duration_second_score_2nd_half_wrong_or_star_list_idate = [] # not really a duration in seconds but a numebr of steps
    
    # -> for each interval iinterv_dur, duration_first_score_2nd_half_wrong_list[iinter_dur] is the number of steps for which the PRN predicted by SpOCK to
    # be highest at the first step of iinterv_dur is actually not selected by the onboard algorithm
    time_second_score_2nd_half_wrong_idate = []
    duration_second_score_2nd_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_and_second_score_2nd_half_wrong_idate = []
    duration_first_and_second_score_2nd_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    time_first_or_second_score_2nd_half_wrong_idate = []
    duration_first_or_second_score_2nd_half_wrong_list_idate = [] # not really a duration in seconds but a numebr of steps
    
    for itime in np.arange(0, ntot-inter_dur_sec, delta_inter):
        print itime, ntot,idate
        ######################################
        port_ant_for_prn_itime = []
        port_ant_prn_list = []

        time_first_score_1st_half_wrong_itime = []
        time_first_score_1st_half_wrong_or_port_itime = []
        time_second_score_1st_half_wrong_or_port_itime = []
        time_second_score_1st_half_wrong_itime = []
        time_first_and_second_score_1st_half_wrong_itime = []
        time_first_or_second_score_1st_half_wrong_itime = []
        time_first_score_2nd_half_wrong_itime = []
        time_first_score_2nd_half_wrong_or_star_itime = []
        time_second_score_2nd_half_wrong_or_star_itime = []
        time_second_score_2nd_half_wrong_itime = []
        time_first_and_second_score_2nd_half_wrong_itime = []
        time_first_or_second_score_2nd_half_wrong_itime = []

        time_first_and_second_score_1st_half_wrong_or_port_itime = []
        time_first_or_second_score_1st_half_wrong_or_port_itime = []
        time_first_and_second_score_2nd_half_wrong_or_star_itime = []
        time_first_or_second_score_2nd_half_wrong_or_star_itime = []

        time_netcdf_two_star_2nd_half_itime = []
        time_netcdf_one_star_2nd_half_itime = []
        
        # BLOCK BELOW IF LOOKING AT BINOMIAL SCORE METRIC
        first_score_idate_temp, second_score_idate_temp,  first_ant_idate_temp, second_ant_idate_temp = select_two_prns()
        first_score_idate.append(first_score_idate_temp)
        second_score_idate.append(second_score_idate_temp)
            ########################################
            #######################################
        # end of BLOCK BELOW IF LOOKING AT BINOMIAL SCORE METRIC
        for iin in range(inter_dur_sec):
            which_ant_netcdf_first_score_here = -1; which_ant_netcdf_second_score_here = -1
            first_score_netcdf_port = -1; second_score_netcdf_port = -1
            first_score_netcdf_star = -1; second_score_netcdf_star = -1
            first_score_selected_onboard = -1; second_score_selected_onboard = -1
            if (first_score_idate[-1] in gps_netcdf_all_date[idate][itime+iin]) == True:
                first_score_selected_onboard = 1
                ispec_first_score_netcdf = np.where(gps_netcdf_all_date[idate][itime+iin] == first_score_idate[-1])[0][0]
                which_ant_netcdf_first_score_here = which_ant_netcdf_all_date[idate][itime+iin][ispec_first_score_netcdf]
                if (which_ant_netcdf_first_score_here == 3): # if this PRN was selected onboard and was assigned to the port antenna
                    first_score_netcdf_port = 1
                elif (which_ant_netcdf_first_score_here == 2): # if this PRN was selected onboard and was assigned to the starboard antenna
                    first_score_netcdf_star = 1
            else:
                first_score_selected_onboard = 0
            if (second_score_idate[-1] in gps_netcdf_all_date[idate][itime+iin]) == True:
                second_score_selected_onboard = 1
                ispec_second_score_netcdf = np.where(gps_netcdf_all_date[idate][itime+iin] == second_score_idate[-1])[0][0]
                which_ant_netcdf_second_score_here = which_ant_netcdf_all_date[idate][itime+iin][ispec_second_score_netcdf]
                if (which_ant_netcdf_second_score_here == 3): # if this PRN was selected onboard and was assigned to the port antenna
                    second_score_netcdf_port = 1
                if (which_ant_netcdf_second_score_here == 2): # if this PRN was selected onboard and was assigned to the starboard antenna
                    second_score_netcdf_star = 1
            else:
                second_score_selected_onboard = 0

            if (iin <= (inter_dur_sec / 2)): # first half of the overpass                    
                if (first_score_selected_onboard == 0):
                    time_first_score_1st_half_wrong_itime.append(iin+itime)
                if ((first_score_selected_onboard == 0) | (first_score_netcdf_port == 1)):
                    time_first_score_1st_half_wrong_or_port_itime.append(iin+itime)
                if (second_score_selected_onboard == 0):
                    time_second_score_1st_half_wrong_itime.append(iin+itime)
                if ((second_score_selected_onboard == 0) | (second_score_netcdf_port == 1)):
                    time_second_score_1st_half_wrong_or_port_itime.append(iin+itime)
                if ((first_score_selected_onboard == 0) & (second_score_selected_onboard == 0)):
                    time_first_and_second_score_1st_half_wrong_itime.append(iin+itime)
                if ((first_score_selected_onboard == 0) | (second_score_selected_onboard == 0)):
                    time_first_or_second_score_1st_half_wrong_itime.append(iin+itime)
                if ( ((first_score_selected_onboard == 0) | (first_score_netcdf_port == 1)) & ((second_score_selected_onboard == 0) | (second_score_netcdf_port == 1))): # if both PRNs were not selected onboard or were selected onboard but assigned to the port antenna
                    time_first_and_second_score_1st_half_wrong_or_port_itime.append(iin+itime)
                if ( (first_score_selected_onboard == 0) | (first_score_netcdf_port == 1) | (second_score_selected_onboard == 0) | (second_score_netcdf_port == 1)): # if one of the 2 PRNs was not selected onboard or was selected onboard but assigned to the port antenna
                    time_first_or_second_score_1st_half_wrong_or_port_itime.append(iin+itime)
            elif (iin >= (inter_dur_sec / 2 + 130)): # 130 s after the start of the second half of the overpass
                if (first_score_selected_onboard == 0):
                    time_first_score_2nd_half_wrong_itime.append(iin+itime)
                if ((first_score_selected_onboard == 0) | (first_score_netcdf_star == 1)):
                    time_first_score_2nd_half_wrong_or_star_itime.append(iin+itime)
                if (second_score_selected_onboard == 0):
                    time_second_score_2nd_half_wrong_itime.append(iin+itime)
                if ((second_score_selected_onboard == 0) | (second_score_netcdf_star == 1)):
                    time_second_score_2nd_half_wrong_or_star_itime.append(iin+itime)
                if ((first_score_selected_onboard == 0) & (second_score_selected_onboard == 0)):
                    time_first_and_second_score_2nd_half_wrong_itime.append(iin+itime)
                if ((first_score_selected_onboard == 0) | (second_score_selected_onboard == 0)):
                    time_first_or_second_score_2nd_half_wrong_itime.append(iin+itime)
                if ( ((first_score_selected_onboard == 0) | (first_score_netcdf_star == 1)) & ((second_score_selected_onboard == 0) | (second_score_netcdf_star == 1))): # if both PRNs were not selected onboard or were selected onboard but assigned to the star antenna
                    time_first_and_second_score_2nd_half_wrong_or_star_itime.append(iin+itime)
                if ( (first_score_selected_onboard == 0) | (first_score_netcdf_star == 1) | (second_score_selected_onboard == 0) | (second_score_netcdf_star == 1)): # if one of the 2 PRNs was not selected onboard or was selected onboard but assigned to the star antenna
                    time_first_or_second_score_2nd_half_wrong_or_star_itime.append(iin+itime)
                # figure out the number of PRNs assigned to the port antenna 130 s after the start of the second half of the overpass 
                nb_netcdf_star = 0
                for ispec_netcdf in range(4):
                    if (which_ant_netcdf_all_date[idate][itime+iin][ispec_netcdf] == 3):
                        nb_netcdf_star = nb_netcdf_star + 1
                if nb_netcdf_star == 2:
                    time_netcdf_two_star_2nd_half_itime.append(iin+itime)
                elif nb_netcdf_star == 1:
                    time_netcdf_one_star_2nd_half_itime.append(iin+itime)

        # if len(time_first_or_second_score_2nd_half_wrong_or_star_itime) == 0:
        #     print first_score_idate_temp, second_score_idate_temp,  first_ant_idate_temp, second_ant_idate_temp
        #     ipdb.set_trace()
        time_first_score_1st_half_wrong_idate.append(time_first_score_1st_half_wrong_itime)
        duration_first_score_1st_half_wrong_list_idate.append(np.float(len(time_first_score_1st_half_wrong_itime)))
        time_first_score_1st_half_wrong_or_port_idate.append(time_first_score_1st_half_wrong_or_port_itime)
        duration_first_score_1st_half_wrong_or_port_list_idate.append(np.float(len(time_first_score_1st_half_wrong_or_port_itime)))
        time_second_score_1st_half_wrong_or_port_idate.append(time_second_score_1st_half_wrong_or_port_itime)
        duration_second_score_1st_half_wrong_or_port_list_idate.append(np.float(len(time_second_score_1st_half_wrong_or_port_itime)))

        time_second_score_1st_half_wrong_idate.append(time_second_score_1st_half_wrong_itime)
        duration_second_score_1st_half_wrong_list_idate.append(np.float(len(time_second_score_1st_half_wrong_itime)))
        time_first_and_second_score_1st_half_wrong_idate.append(time_first_and_second_score_1st_half_wrong_itime)
        duration_first_and_second_score_1st_half_wrong_list_idate.append(np.float(len(time_first_and_second_score_1st_half_wrong_itime)))
        time_first_or_second_score_1st_half_wrong_idate.append(time_first_or_second_score_1st_half_wrong_itime)
        duration_first_or_second_score_1st_half_wrong_list_idate.append(np.float(len(time_first_or_second_score_1st_half_wrong_itime)))

        duration_first_score_1st_half_wrong_list_conc.append(np.float(len(time_first_score_1st_half_wrong_itime)))
        duration_first_score_1st_half_wrong_or_port_list_conc.append(np.float(len(time_first_score_1st_half_wrong_or_port_itime)))
        duration_second_score_1st_half_wrong_or_port_list_conc.append(np.float(len(time_second_score_1st_half_wrong_or_port_itime)))
        duration_second_score_1st_half_wrong_list_conc.append(np.float(len(time_second_score_1st_half_wrong_itime)))
        duration_first_and_second_score_1st_half_wrong_list_conc.append(np.float(len(time_first_and_second_score_1st_half_wrong_itime)))
        duration_first_or_second_score_1st_half_wrong_list_conc.append(np.float(len(time_first_or_second_score_1st_half_wrong_itime)))


        time_first_score_2nd_half_wrong_idate.append(time_first_score_2nd_half_wrong_itime)
        duration_first_score_2nd_half_wrong_list_idate.append(np.float(len(time_first_score_2nd_half_wrong_itime)))
        time_first_score_2nd_half_wrong_or_star_idate.append(time_first_score_2nd_half_wrong_or_star_itime)
        duration_first_score_2nd_half_wrong_or_star_list_idate.append(np.float(len(time_first_score_2nd_half_wrong_or_star_itime)))
        time_second_score_2nd_half_wrong_or_star_idate.append(time_second_score_2nd_half_wrong_or_star_itime)
        duration_second_score_2nd_half_wrong_or_star_list_idate.append(np.float(len(time_second_score_2nd_half_wrong_or_star_itime)))
        time_second_score_2nd_half_wrong_idate.append(time_second_score_2nd_half_wrong_itime)
        duration_second_score_2nd_half_wrong_list_idate.append(np.float(len(time_second_score_2nd_half_wrong_itime)))
        time_first_and_second_score_2nd_half_wrong_idate.append(time_first_and_second_score_2nd_half_wrong_itime)
        duration_first_and_second_score_2nd_half_wrong_list_idate.append(np.float(len(time_first_and_second_score_2nd_half_wrong_itime)))
        time_first_or_second_score_2nd_half_wrong_idate.append(time_first_or_second_score_2nd_half_wrong_itime)
        duration_first_or_second_score_2nd_half_wrong_list_idate.append(np.float(len(time_first_or_second_score_2nd_half_wrong_itime)))

        duration_first_score_2nd_half_wrong_list_conc.append(np.float(len(time_first_score_2nd_half_wrong_itime)))
        duration_first_score_2nd_half_wrong_or_star_list_conc.append(np.float(len(time_first_score_2nd_half_wrong_or_star_itime)))
        duration_second_score_2nd_half_wrong_or_star_list_conc.append(np.float(len(time_second_score_2nd_half_wrong_or_star_itime)))
        duration_second_score_2nd_half_wrong_list_conc.append(np.float(len(time_second_score_2nd_half_wrong_itime)))
        duration_first_and_second_score_2nd_half_wrong_list_conc.append(np.float(len(time_first_and_second_score_2nd_half_wrong_itime)))
        duration_first_or_second_score_2nd_half_wrong_list_conc.append(np.float(len(time_first_or_second_score_2nd_half_wrong_itime)))

        duration_first_and_second_score_1st_half_wrong_or_port_list_conc.append(np.float(len(time_first_and_second_score_1st_half_wrong_or_port_itime)))
        duration_first_or_second_score_1st_half_wrong_or_port_list_conc.append(np.float(len(time_first_or_second_score_1st_half_wrong_or_port_itime)))
        time_first_or_second_score_1st_half_wrong_or_port_idate.append(time_first_or_second_score_1st_half_wrong_or_port_itime)
        duration_first_or_second_score_1st_half_wrong_or_port_list_idate.append(np.float(len(time_first_or_second_score_1st_half_wrong_or_port_itime)))
        time_first_and_second_score_1st_half_wrong_or_port_idate.append(time_first_and_second_score_1st_half_wrong_or_port_itime)
        duration_first_and_second_score_1st_half_wrong_or_port_list_idate.append(np.float(len(time_first_and_second_score_1st_half_wrong_or_port_itime)))

        duration_first_and_second_score_2nd_half_wrong_or_star_list_conc.append(np.float(len(time_first_and_second_score_2nd_half_wrong_or_star_itime)))
        duration_first_or_second_score_2nd_half_wrong_or_star_list_conc.append(np.float(len(time_first_or_second_score_2nd_half_wrong_or_star_itime)))
        time_first_or_second_score_2nd_half_wrong_or_star_idate.append(time_first_or_second_score_2nd_half_wrong_or_star_itime)
        duration_first_or_second_score_2nd_half_wrong_or_star_list_idate.append(np.float(len(time_first_or_second_score_2nd_half_wrong_or_star_itime)))
        time_first_and_second_score_2nd_half_wrong_or_star_idate.append(time_first_and_second_score_2nd_half_wrong_or_star_itime)
        duration_first_and_second_score_2nd_half_wrong_or_star_list_idate.append(np.float(len(time_first_and_second_score_2nd_half_wrong_or_star_itime)))

        duration_netcdf_two_star_2nd_half_list_conc.append(np.float(len(time_netcdf_two_star_2nd_half_itime)))
        time_netcdf_two_star_2nd_half_idate.append(time_netcdf_two_star_2nd_half_itime)
        duration_netcdf_two_star_2nd_half_list_idate.append(np.float(len(time_netcdf_two_star_2nd_half_itime)))
        duration_netcdf_one_star_2nd_half_list_conc.append(np.float(len(time_netcdf_one_star_2nd_half_itime)))
        time_netcdf_one_star_2nd_half_idate.append(time_netcdf_one_star_2nd_half_itime)
        duration_netcdf_one_star_2nd_half_list_idate.append(np.float(len(time_netcdf_one_star_2nd_half_itime)))
        
        
    port_ant_for_prn.append(port_ant_for_prn_idate)
    seconds_sampling_start.append(seconds_sampling_start_idate)
    seconds_sampling_stop.append(seconds_sampling_stop_idate)
    first_and_second_score.append( first_and_second_score_idate )
    first_or_second_score.append( first_or_second_score_idate )
    first_score.append( first_score_idate )
    second_score.append( second_score_idate )
    gap_prn_save.append(gap_prn_save_idate)
    
    time_first_score_1st_half_wrong.append( time_first_score_1st_half_wrong_idate )
    duration_first_score_1st_half_wrong_idate.append( np.array(duration_first_score_1st_half_wrong_list_idate) )
    time_first_score_1st_half_wrong_or_port.append( time_first_score_1st_half_wrong_or_port_idate )
    duration_first_score_1st_half_wrong_or_port_idate.append( np.array(duration_first_score_1st_half_wrong_or_port_list_idate) ) 
    time_second_score_1st_half_wrong_or_port.append( time_second_score_1st_half_wrong_or_port_idate )
    duration_second_score_1st_half_wrong_or_port_idate.append( np.array(duration_second_score_1st_half_wrong_or_port_list_idate) ) 
    time_second_score_1st_half_wrong.append( time_second_score_1st_half_wrong_idate )
    duration_second_score_1st_half_wrong_idate.append( np.array(duration_second_score_1st_half_wrong_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_first_and_second_score_1st_half_wrong.append( time_first_and_second_score_1st_half_wrong_idate )
    duration_first_and_second_score_1st_half_wrong_idate.append( np.array(duration_first_and_second_score_1st_half_wrong_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_first_or_second_score_1st_half_wrong.append( time_first_or_second_score_1st_half_wrong_idate )
    duration_first_or_second_score_1st_half_wrong_idate.append( np.array(duration_first_or_second_score_1st_half_wrong_list_idate) ) # not really a duration in seconds but a numebr of steps
    
    time_first_or_second_score_1st_half_wrong_or_port.append( time_first_or_second_score_1st_half_wrong_or_port_idate )
    duration_first_or_second_score_1st_half_wrong_or_port_idate.append( np.array(duration_first_or_second_score_1st_half_wrong_or_port_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_first_and_second_score_1st_half_wrong_or_port.append( time_first_and_second_score_1st_half_wrong_or_port_idate )
    duration_first_and_second_score_1st_half_wrong_or_port_idate.append( np.array(duration_first_and_second_score_1st_half_wrong_or_port_list_idate) ) # not really a duration in seconds but a numebr of steps

    time_first_or_second_score_2nd_half_wrong_or_star.append( time_first_or_second_score_2nd_half_wrong_or_star_idate )
    duration_first_or_second_score_2nd_half_wrong_or_star_idate.append( np.array(duration_first_or_second_score_2nd_half_wrong_or_star_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_first_and_second_score_2nd_half_wrong_or_star.append( time_first_and_second_score_2nd_half_wrong_or_star_idate )
    duration_first_and_second_score_2nd_half_wrong_or_star_idate.append( np.array(duration_first_and_second_score_2nd_half_wrong_or_star_list_idate) ) # not really a duration in seconds but a numebr of steps

    time_netcdf_two_star_2nd_half.append( time_netcdf_two_star_2nd_half_idate )
    duration_netcdf_two_star_2nd_half_idate.append( np.array(duration_netcdf_two_star_2nd_half_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_netcdf_one_star_2nd_half.append( time_netcdf_one_star_2nd_half_idate )
    duration_netcdf_one_star_2nd_half_idate.append( np.array(duration_netcdf_one_star_2nd_half_list_idate) ) # not really a duration in seconds but a numebr of steps
    
    time_first_score_2nd_half_wrong.append( time_first_score_2nd_half_wrong_idate )
    duration_first_score_2nd_half_wrong_idate.append( np.array(duration_first_score_2nd_half_wrong_list_idate) )
    time_first_score_2nd_half_wrong_or_star.append( time_first_score_2nd_half_wrong_or_star_idate )
    duration_first_score_2nd_half_wrong_or_star_idate.append( np.array(duration_first_score_2nd_half_wrong_or_star_list_idate) )
    time_second_score_2nd_half_wrong_or_star.append( time_second_score_2nd_half_wrong_or_star_idate )
    duration_second_score_2nd_half_wrong_or_star_idate.append( np.array(duration_second_score_2nd_half_wrong_or_star_list_idate) ) 
    time_second_score_2nd_half_wrong.append( time_second_score_2nd_half_wrong_idate )
    duration_second_score_2nd_half_wrong_idate.append( np.array(duration_second_score_2nd_half_wrong_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_first_and_second_score_2nd_half_wrong.append( time_first_and_second_score_2nd_half_wrong_idate )
    duration_first_and_second_score_2nd_half_wrong_idate.append( np.array(duration_first_and_second_score_2nd_half_wrong_list_idate) ) # not really a duration in seconds but a numebr of steps
    time_first_or_second_score_2nd_half_wrong.append( time_first_or_second_score_2nd_half_wrong_idate )
    duration_first_or_second_score_2nd_half_wrong_idate.append( np.array(duration_first_or_second_score_2nd_half_wrong_list_idate) ) # not really a duration in seconds but a numebr of steps

duration_first_score_1st_half_wrong_conc = np.array(duration_first_score_1st_half_wrong_list_conc)
duration_first_score_1st_half_wrong_or_port_conc = np.array(duration_first_score_1st_half_wrong_or_port_list_conc)
duration_second_score_1st_half_wrong_or_port_conc = np.array(duration_second_score_1st_half_wrong_or_port_list_conc)
duration_second_score_1st_half_wrong_conc = np.array(duration_second_score_1st_half_wrong_list_conc)
duration_first_and_second_score_1st_half_wrong_conc = np.array(duration_first_and_second_score_1st_half_wrong_list_conc)
duration_first_or_second_score_1st_half_wrong_conc = np.array(duration_first_or_second_score_1st_half_wrong_list_conc)

duration_first_and_second_score_1st_half_wrong_or_port_conc = np.array(duration_first_and_second_score_1st_half_wrong_or_port_list_conc)
duration_first_or_second_score_1st_half_wrong_or_port_conc = np.array(duration_first_or_second_score_1st_half_wrong_or_port_list_conc)

duration_first_and_second_score_2nd_half_wrong_or_star_conc = np.array(duration_first_and_second_score_2nd_half_wrong_or_star_list_conc)
duration_first_or_second_score_2nd_half_wrong_or_star_conc = np.array(duration_first_or_second_score_2nd_half_wrong_or_star_list_conc)

duration_netcdf_two_star_2nd_half_conc = np.array(duration_netcdf_two_star_2nd_half_list_conc)
duration_netcdf_one_star_2nd_half_conc = np.array(duration_netcdf_one_star_2nd_half_list_conc)

duration_first_score_2nd_half_wrong_conc = np.array(duration_first_score_2nd_half_wrong_list_conc)
duration_first_score_2nd_half_wrong_or_star_conc = np.array(duration_first_score_2nd_half_wrong_or_star_list_conc)
duration_second_score_2nd_half_wrong_or_star_conc = np.array(duration_second_score_2nd_half_wrong_or_star_list_conc)
duration_second_score_2nd_half_wrong_conc = np.array(duration_second_score_2nd_half_wrong_list_conc)
duration_first_and_second_score_2nd_half_wrong_conc = np.array(duration_first_and_second_score_2nd_half_wrong_list_conc)
duration_first_or_second_score_2nd_half_wrong_conc = np.array(duration_first_or_second_score_2nd_half_wrong_list_conc)
