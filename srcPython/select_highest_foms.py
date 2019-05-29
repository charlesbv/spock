# Based on the body elevation and azimuths and the antenna FOM maps, determine the FOM for each SP
# and select the SPs that have the 4 highest FOMs (with a few other small tricks that the onboard algo does)

import numpy as np
import ipdb
def select_highest_foms(elev, azim, elev_not_int, azim_not_int, elev_gps_from_cyg_body, prn, fom_antenna):
    nprn = len(elev)
    antenna_gain_map_azimuth_step = 15
    antenna_gain_map_elevation_step = 5
    antenna_gain_map_width = (360 / antenna_gain_map_azimuth_step )
    refl_elevation_cut_off = 5.0
    fom = np.zeros(nprn)
    # For each PRN Determine the FOM, (highest FOM between the two antennas
    for iprn in range(nprn):
        azimIndex = (azim[iprn] + 180)/antenna_gain_map_azimuth_step
        if (azimIndex == antenna_gain_map_width): # in theory, azimIndex can only vary from 0 to antenna_gain_map_width - 1. But if azim_not_int is between 179.00001 and 180 then azim will be 180 so azimIndex will be antenna_gain_map_width, whicile it should still be antenna_gain_map_width -1, ie the last element of the elev row.
            azimIndex = antenna_gain_map_width - 1
        elevIndex = elev[iprn] / antenna_gain_map_elevation_step
        fom_star = fom_antenna[0, elevIndex, azimIndex ]
        fom_port = fom_antenna[1, elevIndex,  azimIndex ]
        if ((elev_not_int > (90-28)) & (azim_not_int > 0)): # any elevation greater than 28 deg and on the starboard side cannot be seen by the port antenna so the fom port is asigned a gain of -100 (that's what the onboard algoirthm does).
            fom_port = -100
        if ((elev_not_int > (90-28)) & (azim_not_int < 0)): # any elevation greater than 28 deg and on the port side cannot be seen by the starboard antenna so the fom star is asigned a gain of -100 (that's what the onboard algoirthm does)
            fom_star = -100
        # take the highest fom among the two antennas
        max_fom = fom_star
        if fom_port > fom_star:
            max_fom = fom_port
        if (elev_gps_from_cyg_body < refl_elevation_cut_off):
            max_fom = -100
        fom[iprn] = max_fom
    
        
    return 0
