# This script reads in the port and starboard antenna patterns
# and merge them into one single antenna pattern.
# INPUTS:
# - res_map: resolution of the map (set to 'coarse' or 'fine')
# - starboard_antenna_filename: name of the starboard antenna pattern file
# - port_antenna_filename: name of the port antenna pattern file
# OUTPUT:
# - the merged antenna gain pattern is created as follow: for each
# elev/azim bin, the gain is the max gain between the gain of the starboard
# antenna and the gain of the port antenna for this bin
# ASSUMPTION:
# - the starboard and port antenna files have to have the same format as the
# on-board files (if res_map is 'coarse' (.agm files)) or as the files
# created by Scott Gleason (and subequently by Darren McKague) (if res_map
# is 'fine' (.txt files)).
# For a complete description of the format of these files, read the header of 
# the script cygnss_antenna_pattern.py


# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
res_map = 'coarse' # coarse or fine

# # Coarse
starboard_antenna_filename ='/Users/cbv/cspice/data/ant_1_starboard_ddmi_v1.agm'
port_antenna_filename = '/Users/cbv/cspice/data/ant_1_port_ddmi_v1.agm'

# Fine
# starboard_antenna_filename ='/Users/cbv/work/spockOut/beacon/ant_data/ant_1_starboard_v6.txt'
# port_antenna_filename = '/Users/cbv/work/spockOut/beacon/ant_data/ant_1_port_v6.txt'

# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import numpy as np
import ipdb
if ((res_map != 'coarse') & (res_map != 'fine')):
    sys.exit("!***!\n res_map must be set to 'coarse' or 'fine'\
            \n!***!")
                    
print 'Map resolution: ' + res_map

# Read starboard antenna gain
starboard_antenna_file = open(starboard_antenna_filename)
read_starboard = starboard_antenna_file.readlines()
if res_map == 'coarse':
    numEl = len(read_starboard) # should be 18
    numAz = len(read_starboard[0].split(',')) # should be 24
    el_start_deg = 0.
    el_inc_deg = 5.
    el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
    el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)
    az_start_deg = -180.
    az_inc_deg = 15.
    az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
    az_deg = np.linspace(az_start_deg, az_stop_deg, numAz)
    starboard_gain = np.zeros([numEl, numAz])
    for itheta in range(numEl):
        for iphi in range(numAz):
            starboard_gain[itheta, iphi] = np.float(
                read_starboard[itheta].split(',')[iphi] )
else:
    nheader = 2
    numAz = len(read_starboard)-nheader # should be 3601
    numEl = len(read_starboard[nheader].split(',')) # should be 901
    el_start_deg = 90.
    el_inc_deg = -0.1
    el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
    el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)
    az_start_deg = 0.
    az_inc_deg = 0.1
    az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
    az_deg = np.linspace(az_start_deg, az_stop_deg, numAz)
    starboard_gain = np.zeros([numEl, numAz])
    # in the fine formar, rows are azimuths, columns are elevations
    for iphi in range(numAz):
        print iphi, numAz-1
        for itheta in range(numEl):
            starboard_gain[itheta, iphi] = np.float( read_starboard[iphi+
            nheader].split(',')[itheta] )

starboard_antenna_file.close()


# Read port antenna gain
port_antenna_file = open(port_antenna_filename)
read_port = port_antenna_file.readlines()
if res_map == 'coarse':
    numEl = len(read_port) # should be 18
    numAz = len(read_port[0].split(',')) # should be 24
    el_start_deg = 0.
    el_inc_deg = 5.
    el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
    el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)
    az_start_deg = -180.
    az_inc_deg = 15.
    az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
    az_deg = np.linspace(az_start_deg, az_stop_deg, numAz)
    port_gain = np.zeros([numEl, numAz])
    for itheta in range(numEl):
        for iphi in range(numAz):
            port_gain[itheta, iphi] = np.float(
                read_port[itheta].split(',')[iphi] )
else:
    nheader = 2
    numAz = len(read_port)-nheader # should be 3601
    numEl = len(read_port[nheader].split(',')) # should be 901
    el_start_deg = 90.
    el_inc_deg = -0.1
    el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
    el_deg = np.linspace(el_start_deg, el_stop_deg, numEl)
    az_start_deg = 0.
    az_inc_deg = 0.1
    az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
    az_deg = np.linspace(az_start_deg, az_stop_deg, numAz)
    port_gain = np.zeros([numEl, numAz])
    # in the fine formar, rows are azimuths, columns are elevations
    for iphi in range(numAz):
        print iphi, numAz-1
        for itheta in range(numEl):
            port_gain[itheta, iphi] = np.float( read_port[iphi+
            nheader].split(',')[itheta] )

port_antenna_file.close()



# Calculate the merged gain
merged_gain =  np.zeros([numEl, numAz])

starboard_antenna_filename_no_path = starboard_antenna_filename.split('/')[-1]
path_to_starboard_file = starboard_antenna_filename.replace(
    starboard_antenna_filename_no_path, '')
port_antenna_filename_no_path = port_antenna_filename.split('/')[-1]
if res_map == 'coarse':
    merged_antenna_filename = path_to_starboard_file + 'merged_' + \
        starboard_antenna_filename_no_path.replace('.agm', '_with_') + \
        port_antenna_filename_no_path
else:
    merged_antenna_filename = path_to_starboard_file + 'merged_' + \
        starboard_antenna_filename_no_path.replace('.txt', '_with_') + \
        port_antenna_filename_no_path

merged_antenna_file = open(merged_antenna_filename, "w")
if res_map == 'coarse': # rows are elevations, columns are azimuths
    for itheta in range(numEl):
        row = ''
        for iphi in range(numAz):
            if starboard_gain[itheta, iphi] >= port_gain[itheta, iphi]:
                merged_gain[itheta, iphi] = starboard_gain[itheta, iphi]
            else:
                merged_gain[itheta, iphi] =	port_gain[itheta, iphi]
            if iphi != numAz-1:
                row = row + str((int)(merged_gain[itheta, iphi])) + ','
            else:
                row = row + str((int)(merged_gain[itheta, iphi]))
        print >> merged_antenna_file, row
        
else: # rows are azimuths, columns are elevations
    print >> merged_antenna_file,'99,0,0' # first line is the rotation angle
    # of the antenna frame with respect to the body frame but is actually
    # ignore (according to FileFormatNotes.txt) so set to 99,0,0
    print >> merged_antenna_file, '0,3601,0.1,0,901,0.1' # this is always the
    # first line of
    # the fine files (although the start el/az and numEl/Az are set earlier in
    # the script, they are actually ignored (actually I don't know why but the
    # start value of elevation is 90 but set to 0 at this line (in the files
    # written by Scott and Darren)
    for iphi in range(numAz):
        print iphi, numAz-1
        row = ''
        for itheta in range(numEl):
            if starboard_gain[itheta, iphi] >= port_gain[itheta, iphi]:
                merged_gain[itheta, iphi] = starboard_gain[itheta, iphi]
            else:
                merged_gain[itheta, iphi] =	port_gain[itheta, iphi]
            if itheta != numEl-1:
                row = row + str(merged_gain[itheta, iphi]) + ','
            else:
                row = row + str(merged_gain[itheta, iphi])
        print >> merged_antenna_file, row
        
merged_antenna_file.close()
