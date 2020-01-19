# this script converts a spp3 file from https://swarm-diss.eo.esa.int/#swarm%2FLevel1b%2FEntire_mission_data%2FGPSxNAV into ECI coordinates in a format similar to the output format of cygnss_read_netcdf_to_eci_observation.py
# attitude information is availble at: https://swarm-diss.eo.esa.int/#swarm%2FLevel1b%2FEntire_mission_data%2FSTRxATT. The script that reads the attitude files and convert them into quaternion files for SpOCK in on Big in ~/cde/swarm/swarm_attitude_to_quaternion_spock.py

import numpy as np
from datetime import datetime, timedelta
import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from ecef2eci import *
import os

dir_sp3 = './20170901_20170910'

filename_list = []
for file in os.listdir(dir_sp3):
    if file.endswith(".sp3"):
        filename_list.append(os.path.join(dir_sp3, file))

nfile = len(filename_list)
first_date = []
first_date_raw = []
for ifile in range(nfile):
    date_temp = filename_list[ifile].split('/')[-1].split('_1B_')[1].split('_')[0]
    first_date_raw.append(date_temp)
    date_temp_date = datetime.strptime(date_temp, "%Y%m%dT%H%M%S")
    first_date.append(date_temp_date)
first_date_arr = np.array(first_date)
date_sort_index =  np.argsort(first_date_arr)

filename_eci  = filename_list[ifile].split('/')[-1].split('_1B_')[0] + '_1B_' + first_date_raw[date_sort_index[0]] + '_' + first_date_raw[date_sort_index[-1]] + '_eci.txt'
file_eci = open(filename_eci, 'w')
print >> file_eci,'#Date position(km/s) velocity(km/s)'
print >> file_eci,'#START'

for ifile in date_sort_index:
    filename_sp3 = filename_list[ifile]
    print filename_sp3
    file_sp3 = open(filename_sp3)
    read_file_sp3 = file_sp3.readlines()
    #skip header
    iheader = 0
    while read_file_sp3[iheader][:5] != '*  20':
        iheader = iheader + 1
    nsp3 = len(read_file_sp3) - iheader - 1
    ntime = nsp3 / 3

    for itime in range(0, ntime):
        # figure out date
        date_sp3_raw = read_file_sp3[itime*3 + iheader].split('*')[1].replace('\n','').replace('\r','')
        yy = (date_sp3_raw.split()[0])
        mm = (date_sp3_raw.split()[1])
        dd = (date_sp3_raw.split()[2])
        hr = (date_sp3_raw.split()[3])
        minute = (date_sp3_raw.split()[4])
        sec = (date_sp3_raw.split()[5].split('.')[0])
        mil = (date_sp3_raw.split()[5].split('.')[1])
        bla = 0
        if sec == '60':
            bla = 1
            minute = str((int)(minute) + 1)
            sec = '0'
            #print date_sp3_raw
            #raise Exception
        date_sp3_str_temp = yy + '-' + mm + '-' + dd + 'T' + hr + ':' + minute + ':' + sec

        date_sp3_temp = datetime.strptime(date_sp3_str_temp, "%Y-%m-%dT%H:%M:%S")
        if mil[:5] == '99999':
            date_sp3 = date_sp3_temp + timedelta(seconds = 1)
        elif mil[:5] == '00000':
            date_sp3 = date_sp3_temp
        else:
            print read_file_sp3[itime*3 + iheader]
            raise Exception
        if bla == 1:
            print date_sp3_raw, date_sp3, 'XXXXXXXXX'
        # else:
        #     print date_sp3_raw, date_sp3

        date_sp3_str = datetime.strftime(date_sp3, "%Y-%m-%dT%H:%M:%S")
        # figure out ecef position
        r_sp3_ecef = np.array([np.float(i) for i in read_file_sp3[itime*3 + 1 + iheader].split()[1:4]]) # km

        # figure out ecef velocity
        v_sp3_ecef = np.array([np.float(i) for i in read_file_sp3[itime*3 + 2 + iheader].split()[1:4]]) # decimeter/s (!!!)
        v_sp3_ecef = v_sp3_ecef / 10000. # km/s
        # convert ECEF to ECI coordinates
        if itime == 0:
            load_spice = 1
        else:
            load_spice = 0
        r_sp3_eci, v_sp3_eci = ecef2eci(r_sp3_ecef, v_sp3_ecef, date_sp3_str, load_spice)
        print >> file_eci, date_sp3_str, r_sp3_eci[0], r_sp3_eci[1], r_sp3_eci[2], v_sp3_eci[0], v_sp3_eci[1], v_sp3_eci[2]
        #file_eci.close() # !!!!!! remove
        #raise Exception
file_eci.close() 
