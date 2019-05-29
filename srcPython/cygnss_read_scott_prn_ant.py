# This script reads the file CYGNSS_yaw_test_DMR_satselection.dat prepared by Scott Gleason. This file incldues the PRNs and the assigned antennas on Sep 26 2018 for FM02

import numpy as np
filename_scott = '/Users/cbv/cygnss/scott/CYGNSS_yaw_test_DMR_satselection.dat'
file_scott = open(filename_scott)
read_file_scott = file_scott.readlines()
nscott = len(read_file_scott)
nb_seconds_since_start_scott = np.zeros([nscott])
gps_scott = np.zeros([nscott, 4])
which_ant_scott = np.zeros([nscott, 4])
for i in range(nscott):
    nb_seconds_since_start_scott[i] = np.float(read_file_scott[i].split()[0])
    for ispec in range(4):
        gps_scott[i, ispec] = np.float(read_file_scott[i].split()[3+ispec*2])
        which_ant_scott[i, ispec] = np.float(read_file_scott[i].split()[4+ispec*2])
