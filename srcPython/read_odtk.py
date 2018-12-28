# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
# This script reads the output create by ODTK and sent by Kyle Nave's team. Inspired from find_rho.py

import numpy as np
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
from datetime import datetime, timedelta

def read_odtk(filename_odtk, date_start_str, date_stop_str):
# IN
# filename_odtk = '/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/eclipse/cyg_data/F7_20170824_145117_STKdefpred_v001.e'
# date_start_str = "2017-08-20T12:00:00" # date to start saving r_odtk and v_odtk 
# date_stop_str = "2017-08-20T18:00:00" # date to stop saving r_odtk and v_odtk 


    # ALGORITHM
    while True:
        try:
            date_start = datetime.strptime(date_start_str, "%Y-%m-%dT%H:%M:%S")
            break
        except ValueError:
            try: 
                date_start = datetime.strptime(date_start_str, "%Y-%m-%dT%H:%M:%S.%f")
                break
            except ValueError:
                print '***!The start date has to follow the format "%Y-%m-%dT%H:%M:%S" or "%Y-%m-%dT%H:%M:%S.%f". The program will stop.!***'
                raise Exception
    while True:
        try:
            date_stop = datetime.strptime(date_stop_str, "%Y-%m-%dT%H:%M:%S")
            break
        except ValueError:
            try: 
                date_stop = datetime.strptime(date_stop_str, "%Y-%m-%dT%H:%M:%S.%f")
                break
            except ValueError:
                print '***!The stop date has to follow the format "%Y-%m-%dT%H:%M:%S" or "%Y-%m-%dT%H:%M:%S.%f". The program will stop.!***'
                raise Exception

    file = open(filename_odtk)
    read_file = file.readlines()
    iline = 0
    found_epoch = 0
    while (found_epoch == 0):
        if ('YYDDD' in read_file[iline]):
            epochyyddd = np.float(read_file[iline].split(":")[1])
            found_epoch = 1
        iline = iline + 1

    found_hdr = 0
    while (found_hdr == 0):
        if (len(read_file[iline].split()) > 0):
            if read_file[iline].split()[0] == "EphemerisTimePosVel":
                found_hdr = 1
        iline = iline + 1
    nb_header = iline + 1

    n = len(read_file) - nb_header
    date_odtk = []
    date_odtk_raw = []
    r_odtk = []
    v_odtk = []

    epochyyddd_no_decimal = (int)(epochyyddd)
    epochyyddd_decimal_only = epochyyddd - epochyyddd_no_decimal
    epochyyddd_date = datetime.strptime( str(epochyyddd_no_decimal), "%y%j" ) + timedelta( hours = epochyyddd_decimal_only*24 )

    iline = 0
    while (iline < n):
        if len(date_odtk) > 0:
            if date_odtk[-1] >= date_stop:
                break
        if len(read_file[iline+nb_header].split()) == 0:
            break
        date_odtk_temp = np.float( read_file[iline+nb_header].split()[0] )
        if epochyyddd_date + timedelta( seconds = date_odtk_temp ) >= date_start:
            date_odtk_raw.append(read_file[iline+nb_header].split()[0] )
            date_odtk.append( epochyyddd_date + timedelta( seconds = date_odtk_temp ) )
            r_odtk.append( [np.float(read_file[iline+nb_header].split()[1]), np.float(read_file[iline+nb_header].split()[2]), np.float(read_file[iline+nb_header].split()[3])] )
            v_odtk.append( [np.float(read_file[iline+nb_header].split()[4]), np.float(read_file[iline+nb_header].split()[5]), np.float(read_file[iline+nb_header].split()[6])] )
        iline = iline + 1

    r_odtk = np.array(r_odtk)/1000.
    v_odtk = np.array(v_odtk)/1000.
    date_odtk = np.array(date_odtk)
    date_odtk_raw = np.array(date_odtk_raw)
    
    return r_odtk, v_odtk, date_odtk, date_odtk_raw # date_odtk_raw is the date as it appears in the odtk file (in seconds sicne epoch)
