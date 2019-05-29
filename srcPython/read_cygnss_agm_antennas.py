# THis script reqads the two antenna FOM coarsemaps (agm fiels)
import numpy as np
def read_cygnss_agm_antennas():
    filename_gain_list = ['ant_1_starboard_ddmi_v1.agm', 'ant_1_port_ddmi_v1.agm']
    nant = len(filename_gain_list)
    for iant in range(nant):
        nheader = 0
        filename_gain = '/Users/cbv/cspice/data/' + filename_gain_list[iant]
        file_gain = open(filename_gain, "r")
        read_file_gain = file_gain.readlines()
        numEl = len(read_file_gain) - nheader
        numAz = len(read_file_gain[0+nheader].split(','))
        # Elevation
        el_start_deg_file = 0. 
        el_inc_deg = 5. # FileFormatNotes.txt: step is 0.1 deg
        el_stop_deg_file = el_start_deg_file + el_inc_deg * (numEl-1)
        # in FileFormatNotes.txt, elev goes in decreasing order
        el_deg_file = np.linspace(el_start_deg_file, el_stop_deg_file, numEl)
        # Azimuth
        az_start_deg = -180. 
        az_inc_deg =  15. # FileFormatNotes.txt: step is 0.1 deg
        az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
        az_deg_file = np.linspace(az_start_deg, az_stop_deg, numAz)
        # Gain
        if iant == 0: # numEl, numAz same for both ants
            gain_file = np.zeros([nant, numEl, numAz])

        for itheta in range(numEl):
            for iphi in range(numAz):
                gain_file[iant,itheta, iphi] = np.float( read_file_gain[itheta+nheader].split(',')[iphi] )
    return gain_file
