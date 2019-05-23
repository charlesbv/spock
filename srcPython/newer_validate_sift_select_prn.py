score_prn = np.zeros([33,33])
score_first_half_prn = np.zeros([33,33])
gap_prn = np.zeros([33, 33]) + 1e6 # set to 1e6 so that prns that are not selected have a very big gap
gap_first_half_prn = np.zeros([33, 33]) + 1e6 # set to 1e6 so that prns that are not selected have a very big gap
single_gap_prn_port_second_half = np.zeros([33]) + 1e6 # record the gap of the PRN taht are port starting at the second half  + 160s
prn_list = [] # record all PRNs during the entire time
prn_port_second_half_list =  [] # record the PRN taht are port starting at the second half  + 160s 
for iin in range(inter_dur_sec):
    for ispec in range(4):
        if ( gps_spock_all_date[idate][itime+iin][ispec] in prn_list ) == False:
            prn_list.append(gps_spock_all_date[idate][itime+iin][ispec])
        if ((iin >= (inter_dur_sec/2 + 130 + 30)) & (which_ant_spock_all_date[idate][itime+iin][ispec] == 3) & (( gps_spock_all_date[idate][itime+iin][ispec] in prn_port_second_half_list ) == False)):
            prn_port_second_half_list.append(gps_spock_all_date[idate][itime+iin][ispec])
prn_list = np.array(prn_list)
nprn = len(prn_list)
prn_list_sort = prn_list[np.argsort(prn_list)]
#  for the prns that are selected, innitialize the gap to 0 
for iii in range(nprn):
    if prn_list_sort[iii] in prn_port_second_half_list: # if this PRN is port starting at the second half  + 160s 
        single_gap_prn_port_second_half[prn_list_sort[iii]] = 0
    for jjj in range(nprn):
        if iii != jjj:
            gap_prn[prn_list_sort[iii], prn_list_sort[jjj]] = 0
            gap_first_half_prn[prn_list_sort[iii], prn_list_sort[jjj]] = 0 


seconds_sampling_start_idate.append(nb_seconds_since_initial_epoch_spock_all_date[idate][itime])
seconds_sampling_stop_idate.append(nb_seconds_since_initial_epoch_spock_all_date[idate][itime+inter_dur_sec-1])
for iin in range(inter_dur_sec): #array([ 7,  8, 11, 16, 18, 27])
    iout = -1
    port_ant_for_prn_iprn = []
    # first, determine if,  starting at the second hald of the overpass + 160 s, the port antenna is assigned to a PRN
    if iin >= (inter_dur_sec/2 + 130 + 30): # starting at the second hald of the overpass + 160 s, look if the port antenna is assigned to a PRN
        for prn_out in prn_list_sort:
            prn_out_is_gap = 0 # used only to figure out the gap for PRN that are port starting at the second hald of the overpass + 160 s  
            if len(np.where(gps_spock_all_date[idate][itime+iin] == prn_out)[0]) > 0: #the prn is selected by SpOCK at this particular time            
                iprn_out = np.where(gps_spock_all_date[idate][itime+iin] == prn_out)[0][0]
                if which_ant_spock_all_date[idate][itime+iin][iprn_out] == 3: # the port antenna is assigned to a PRN
                    port_ant_for_prn_iprn.append([itime, itime+iin, gps_spock_all_date[idate][itime+iin][iprn_out]])
                    if ((gps_spock_all_date[idate][itime+iin][iprn_out] in port_ant_prn_list) == False): # record only once the PRN that's assigned to port
                        port_ant_prn_list.append(gps_spock_all_date[idate][itime+iin][iprn_out]) # actually this is probably the same list as prn_port_second_half_list. oh well
            if (prn_out in prn_port_second_half_list):  # if the PRN is port at some point starting at the second hald of the overpass + 160 s
                if (len(np.where(gps_spock_all_date[idate][itime+iin] == prn_out)[0]) == 0 ): # if there is a gap for this prn
                    prn_out_is_gap = 1
                else:
                    iprn_out = np.where(gps_spock_all_date[idate][itime+iin] == prn_out)[0][0]
                    if ( fom_spock_all_date[idate][itime+iin][iprn_out] == 0): # or if its FOM is 0 (count it as a gap)
                        prn_out_is_gap = 1
                    elif (which_ant_spock_all_date[idate][itime+iin][iprn_out] == 2): # or at this time the antenna is starboard  (this prn is port at some point starting at the second hald of the overpass + 160 s but not neceessarily at all times starting at the second hald of the overpass + 160 s)
                        prn_out_is_gap = 1
                if prn_out_is_gap == 1:
                    single_gap_prn_port_second_half[prn_out] = single_gap_prn_port_second_half[prn_out] + 1

                        
    # now compute the score and gap for every combination
    for prn_out in prn_list_sort[:-1]: # no need to look at the last element since all combinations ahve already been considered #array([ 7,  8, 11, 16, 18, 27])
        prn_out_is_gap = 0
        #ipdb.set_trace()
        iout = iout + 1
        ant_out = -1
        if len(np.where(gps_spock_all_date[idate][itime+iin] == prn_out)[0]) > 0: #the prn is selected by SpOCK at this particular time            
            iprn_out = np.where(gps_spock_all_date[idate][itime+iin] == prn_out)[0][0]            
            gain_out = fom_spock_all_date[idate][itime+iin][iprn_out]
            ant_out = which_ant_spock_all_date[idate][itime+iin][iprn_out]
            if gain_out == 0:# if the gain is 0, count it as a gap (since SpOCK is very different from onboard for gains of 0) 
                prn_out_is_gap = 1
            elif ((iin <= inter_dur_sec/2) & (ant_out == 3)): # if the antenna is port during the first half of the overpass, count it as a gap since the beacon will only transmit to the starboard antenna
                prn_out_is_gap = 1
                
                # TOTO # if one of the two PRNs is assigned to port during the first hallf of the overpass then count it as a gap
        else:
            gain_out = 0 # !!!!!! used ot be -1 to penalize non selected prn
            prn_out_is_gap = 1

        for prn_in in prn_list_sort[iout+1:]:
            prn_in_is_gap = 0
            ant_in = -1
            if len(np.where(gps_spock_all_date[idate][itime+iin] == prn_in)[0]) > 0: #the prn is selected by SpOCK at this particular time
                iprn_in = np.where(gps_spock_all_date[idate][itime+iin] == prn_in)[0][0]
                gain_in = fom_spock_all_date[idate][itime+iin][iprn_in]
                ant_in = which_ant_spock_all_date[idate][itime+iin][iprn_in]
                max_gain_out_in = np.max([gain_out, gain_in])
                if gain_in == 0: # if the gain is 0, count it as a gap (since SpOCK is very different from onboard for gains of 0)
                    prn_in_is_gap = 1
                elif ((iin <= inter_dur_sec/2) & (ant_in == 3)): # if the antenna is port during the first half of the overpass, count it as a gap since the beacon will only transmit to the starboard antenna
                    prn_in_is_gap = 1
            else:
                gain_in = 0 # !!!!!! used ot be -1 to penalize non selected prn
                max_gain_out_in = np.max([gain_out, gain_in])
                prn_in_is_gap = 1
            score_prn[prn_out, prn_in] = score_prn[prn_out, prn_in] + max_gain_out_in
            if iin <= inter_dur_sec / 2:
                score_first_half_prn[prn_out, prn_in] = score_first_half_prn[prn_out, prn_in] + max_gain_out_in
            if ((prn_out_is_gap == 1) | (prn_in_is_gap == 1)):
                gap_prn[prn_out, prn_in] = gap_prn[prn_out, prn_in] + 1
                if (iin <= inter_dur_sec / 2):
                    gap_first_half_prn[prn_out, prn_in] = gap_first_half_prn[prn_out, prn_in] + 1
            if ((prn_out_is_gap == 1) & (prn_in_is_gap == 1)):
                gap_prn[prn_out, prn_in] = gap_prn[prn_out, prn_in] + 10000 # we won't to exclude the possiblity of choosing this combination since both prn have the gap at the same time
                if iin <= inter_dur_sec	/ 2:
                    gap_first_half_prn[prn_out, prn_in] = gap_first_half_prn[prn_out, prn_in] + 10000 # we won't to exclude the possiblity of choosing this combination since both prn have the gap at the same time
    if len(port_ant_prn_list) > 0:
        port_ant_for_prn_itime.append(port_ant_for_prn_iprn)

port_ant_for_prn_idate.append(port_ant_for_prn_itime)
ncomb = len(np.where(gap_first_half_prn != 1e6)[0]) # should be equal to (nprn*(nprn-1))/2# total number of combinaiton. /2 because combinaiton [X,Y] is the same as [Y,X].
gap_first_half_prn_symm = gap_first_half_prn + np.transpose(gap_first_half_prn)
score_first_half_prn_symm = score_first_half_prn + np.transpose(score_first_half_prn)
if len(port_ant_prn_list) != 0: # This means that one of the PRN was assigneed to the port antenna starting at the second hald of the overpass + 160 s
    if len(port_ant_prn_list) == 1: # exactly onr PRN was assigned to port
        prn_that_is_port = port_ant_prn_list[0]
        second_score_temp = prn_that_is_port
        second_ant_temp = 3
        # determine which other PRN monimized the gap with this PRN during the firsrt half of the overpass
        ## there cuold be more than one combination that minmizes the gap, in which case choose the combination that has the highest score_prn
        list_gap = np.argsort(gap_first_half_prn_symm[:, prn_that_is_port])
        nmin = len(list_gap)
        imin = 0
        found_imin = 0
        min_gap = gap_first_half_prn_symm[list_gap[imin], prn_that_is_port]
        list_min_gap = []
        score_prn_list_min_gap = []
        list_min_gap.append(list_gap[imin])
        score_prn_list_min_gap.append(score_first_half_prn_symm[list_gap[imin], prn_that_is_port])
        while ((imin < nmin) & (found_imin == 0)):
            imin = imin + 1
            if gap_first_half_prn_symm[list_gap[imin], prn_that_is_port] == min_gap: # more than one combination minimze the gap
                list_min_gap.append(list_gap[imin])
                score_prn_list_min_gap.append(score_first_half_prn_symm[list_gap[imin], prn_that_is_port])
            else:
                found_imin = 1
        score_prn_list_min_gap = np.array(score_prn_list_min_gap)
        if len(list_min_gap) > 1: # there are more than one combination that minmizes the gap, in which case choose the combination that has the highest score_first_half_prn
            iprn_min_gap_with_this_port_prn = np.where(score_prn_list_min_gap == np.max(score_prn_list_min_gap))[0][0]
            prn_min_gap_with_this_port_prn = list_min_gap[iprn_min_gap_with_this_port_prn]
        else:
            prn_min_gap_with_this_port_prn = list_min_gap[0]
        first_score_temp = prn_min_gap_with_this_port_prn
        first_ant_temp = 2
    elif len(port_ant_prn_list) >= 2: # two or more PRNs were assigned to port. Select the PRN port that has the min gap start at the scond half + 160s. For the second PRN, select so that the combinatino minimizes the gap with this port prn
        prn_port_min_gap_second_half = np.where(single_gap_prn_port_second_half == np.min(single_gap_prn_port_second_half))[0][0] # Select the PRN port that has the min gap start at the scond half + 160s
        second_score_temp = prn_port_min_gap_second_half
        second_ant_temp = 3
        list_gap = np.argsort(gap_first_half_prn_symm[:, prn_port_min_gap_second_half])
        nmin = len(list_gap)
        imin = 0
        found_imin = 0
        min_gap = gap_first_half_prn_symm[list_gap[imin], prn_port_min_gap_second_half]
        list_min_gap = []
        score_prn_list_min_gap = []
        list_min_gap.append(list_gap[imin])
        score_prn_list_min_gap.append(score_first_half_prn_symm[list_gap[imin], prn_port_min_gap_second_half])
        while ((imin < nmin) & (found_imin == 0)):
            imin = imin + 1
            if gap_first_half_prn_symm[list_gap[imin], prn_port_min_gap_second_half] == min_gap: # more than one combination minimze the gap
                list_min_gap.append(list_gap[imin])
                score_prn_list_min_gap.append(score_first_half_prn_symm[list_gap[imin], prn_port_min_gap_second_half])
            else:
                found_imin = 1
        score_prn_list_min_gap = np.array(score_prn_list_min_gap)
        if len(list_min_gap) > 1: # there are more than one combination that minmizes the gap, in which case choose the combination that has the highest score_first_half_prn
            iprn_min_gap_with_this_port_prn = np.where(score_prn_list_min_gap == np.max(score_prn_list_min_gap))[0][0]
            prn_min_gap_with_this_port_prn = list_min_gap[iprn_min_gap_with_this_port_prn]
        else:
            prn_min_gap_with_this_port_prn = list_min_gap[0]
        first_score_temp = prn_min_gap_with_this_port_prn
        first_ant_temp = 2

else: # this means that the port antenna was never assigned to any PRN starting at the second hald of the overpass + 160 s
    score_index_sort_temp = np.dstack(np.unravel_index(np.argsort(score_first_half_prn.ravel()), score_first_half_prn.shape))
    score_index_sort = score_index_sort_temp[0, :, :] # sorted array of combinations that give the ghihest score (ascending order)
    icomb = -1
    #found_comb_without_gap = 0
    gap_prn_here = np.zeros([ncomb])
    while icomb >= -ncomb: # go through score_index_sort from the combination that gives the higeshest score (score_index_sort[-1,:]) to the combination that gives the lowest score (score_index_sort[-ncomb, :])
        comb_now = score_index_sort[icomb, :]
        prn_out_here = comb_now[0]
        prn_in_here = comb_now[1]
        gap_prn_here[icomb] = gap_first_half_prn_symm[prn_out_here, prn_in_here]
        icomb = icomb-1
    #if found_comb_without_gap == 0: # if none of the combinations had no gap (ie if all combinations had at least one second gap) then take the combination with the smallest amount of gap
    list_gap = np.where(gap_prn_here == np.min(gap_prn_here))[0] # score_index_sort[-ncomb + list_gap, :] is the list of combinations that minimize the gap
    min_gap = np.min(gap_prn_here)
    nmin = len(list_gap)
    imin = 0
    found_imin = 0
    list_min_gap = []
    score_prn_list_min_gap = []
    combi_now = score_index_sort[-ncomb + list_gap[imin]]
    list_min_gap.append(combi_now)
    score_prn_list_min_gap.append(score_first_half_prn_symm[combi_now[0], combi_now[1]])
    imin = imin + 1
    while ((imin < nmin) & (found_imin == 0)): # more than one combination minimze the gap
        combi_now =	score_index_sort[-ncomb + list_gap[imin]]
        if gap_first_half_prn_symm[combi_now[0], combi_now[1]] == min_gap: # more than one combination minimze the gap
            list_min_gap.append(combi_now)
            score_prn_list_min_gap.append(score_first_half_prn_symm[combi_now[0], combi_now[1]])
        else:
            found_imin = 1
        imin = imin + 1
    score_prn_list_min_gap = np.array(score_prn_list_min_gap)
    iprn_min_gap_with_this_port_prn = np.where(score_prn_list_min_gap == np.max(score_prn_list_min_gap))[0][0]
    optim_comb = list_min_gap[iprn_min_gap_with_this_port_prn]

    first_score_temp = optim_comb[0]
    second_score_temp = optim_comb[1]
    first_ant_temp = 2
    second_ant_temp = 2
