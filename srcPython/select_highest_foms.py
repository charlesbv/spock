# Based on the body elevation and azimuths and the antenna FOM maps, determine the FOM for each SP
# and select the SPs that have the 4 highest FOMs (with a few other small tricks that the onboard algo does)

import numpy as np
import ipdb
def select_highest_foms(elev, azim, elev_not_int, azim_not_int, elev_gps_from_cyg_body, prn, fom_antenna, prn_previous_step, which_ant_previous_step, reproduce_onboard_bug_antenna_selection):
    # reproduce_onboard_bug_antenna_selection: if set to 1 then the PRN and antenna selection algorithm reproduces the onboard bug which doesn't properly switch from fore and aft antenna when the FM is yawed by +-90 deg
    # prn_previous_step: only needed if reproduce_onboard_bug_antenna_selection == 1. list of PRNs selected at the previous step
    # which_ant_previous_step: only needed if reproduce_onboard_bug_antenna_selection == 1. list of antennas assigned to each PRN at the previous step
    nprn_max = 4
    nprn = len(elev)
    antenna_gain_map_azimuth_step = 15
    antenna_gain_map_elevation_step = 5
    antenna_gain_map_width = (360 / antenna_gain_map_azimuth_step )
    refl_elevation_cut_off = 5.0
    fom = np.zeros(nprn); which_ant = np.zeros(nprn);
    # For each PRN Determine the FOM, (highest FOM between the two antennas
    for iprn in range(nprn):
        azimIndex = (azim[iprn] + 180)/antenna_gain_map_azimuth_step
        if (azimIndex == antenna_gain_map_width): # in theory, azimIndex can only vary from 0 to antenna_gain_map_width - 1. But if azim_not_int is between 179.00001 and 180 then azim will be 180 so azimIndex will be antenna_gain_map_width, whicile it should still be antenna_gain_map_width -1, ie the last element of the elev row.
            azimIndex = antenna_gain_map_width - 1
        elevIndex = elev[iprn] / antenna_gain_map_elevation_step
        fom_star = fom_antenna[0, elevIndex, azimIndex ]
        fom_port = fom_antenna[1, elevIndex,  azimIndex ]
        if ((elev_not_int[iprn] > (90-28)) & (azim_not_int[iprn] > 0)): # any elevation greater than 28 deg and on the starboard side cannot be seen by the port antenna so the fom port is asigned a gain of -100 (that's what the onboard algoirthm does).
            fom_port = -100
        if ((elev_not_int[iprn] > (90-28)) & (azim_not_int[iprn] < 0)): # any elevation greater than 28 deg and on the port side cannot be seen by the starboard antenna so the fom star is asigned a gain of -100 (that's what the onboard algoirthm does)
            fom_star = -100
        # take the highest fom among the two antennas
        max_fom = fom_star
        which_ant_here = 2 # starboard antenna
        # if prn[iprn] == 20:
        #     ipdb.set_trace()
        if fom_port > fom_star: # port is higher FOM
            max_fom = fom_port
            which_ant_here = 3 # port antenna
        if (elev_gps_from_cyg_body[iprn] < refl_elevation_cut_off):
            max_fom = -100
            
        fom[iprn] = max_fom
        which_ant[iprn] = which_ant_here
    index_for_sort = np.argsort(-fom)
    index_highest_fom = index_for_sort[:nprn_max]
    prn_selected = np.array(prn)[index_highest_fom]
    fom_selected = fom[index_highest_fom]
    which_ant_selected = which_ant[index_highest_fom]
    elev_selected = np.array(elev)[index_highest_fom]
    azim_selected = np.array(azim)[index_highest_fom]
    elev_not_int_selected = np.array(elev_not_int)[index_highest_fom]
    azim_not_int_selected = np.array(azim_not_int)[index_highest_fom]
    for iselected in range(nprn_max):
        if azim_selected[iselected] < 0:
            azim_selected[iselected] = 360 + azim_selected[iselected]
        if azim_not_int_selected[iselected] < 0:
            azim_not_int_selected[iselected] = 360 + azim_not_int_selected[iselected]

            
        if (reproduce_onboard_bug_antenna_selection == 1):
            if (prn_selected[iselected] in prn_previous_step): # this means that the GPS was selected in the previous step
                iprn_previous = np.where(prn_previous_step == prn_selected[iselected])[0][0]
                which_ant_previous = which_ant_previous_step[iprn_previous]
                if ((which_ant_previous == 2) & (which_ant_selected[iselected] == 3)): # if the GPS was selected in the previous step and that for the current step it is assigned to the port antenna while it was assigned the starboard antenna in the previous step, then change the assigned antenna to the starboard antenna (ie, same as previous step. This is to reproduce the bug in the CYGNSS onboard algorithm
                    which_ant_selected[iselected] = 2


            
    return prn_selected, fom_selected, which_ant_selected, elev_selected, azim_selected, elev_not_int_selected, azim_not_int_selected
