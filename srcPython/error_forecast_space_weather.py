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



import os
from os import listdir
import matplotlib.colors as colors
import matplotlib.cm as cmx
from os.path import isfile, join
from datetime import datetime, timedelta
import numpy as np
from matplotlib  import pyplot as plt
import sys
import ipdb
# format = 'swpc_f107' # swpc_f107, swpc_ap
def read_obs(format, filename, date_start_str, date_stop_str):
    print 'Format of the observation file: ' + format
    print 'Filename: ' + filename
    print 'Start date: ' + date_start_str
    print 'End date: ' + date_stop_str
    date_start = datetime.strptime(date_start_str, "%Y-%m-%d")
    date_stop = datetime.strptime(date_stop_str, "%Y-%m-%d")
    file_obs = open(filename)
    read_file_obs = file_obs.readlines()
    n_hdr_obs = 0
    while (read_file_obs[n_hdr_obs][0] != date_start_str[0]):
        n_hdr_obs = n_hdr_obs + 1
    nobs = (int)((date_stop - date_start).total_seconds() / 3600. / 24) + 1 # should be equal to n_pred
    obs = np.zeros([nobs])
    ## Skip observations prior to the start date of the historical predictions
    iobs = 0
    date_obs_temp = read_file_obs[n_hdr_obs + iobs].split()[0] + read_file_obs[n_hdr_obs + iobs].split()[1] + read_file_obs[n_hdr_obs + iobs].split()[2]
    date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
    while ( date_obs != date_start ):
        iobs = iobs + 1
        try:
            date_obs_temp = read_file_obs[n_hdr_obs + iobs].split()[0] + read_file_obs[n_hdr_obs + iobs].split()[1] + read_file_obs[n_hdr_obs + iobs].split()[2]
        except IndexError:
            ipdb.set_trace()
        date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
    obs_start = iobs
    date_obs = [] # !!!!!! if everything goes well then: date_obs = date_pred. Otherwise, there is a problem somewhere so debug!
    skip_lines = 0
    if (format == 'swpc_f107'):
        col_start = 12; col_stop = 15
    elif (format == 'swpc_ap'):
        col_start = 59; col_stop = 62
    for iobs in range( obs_start, obs_start + nobs):
        if ((format == 'swpc_f107') | (format == 'swpc_ap')):
            if ( read_file_obs[skip_lines + n_hdr_obs + iobs].split()[0] != ":Product:" ):
                date_obs_temp = read_file_obs[skip_lines + n_hdr_obs + iobs].split()[0] + read_file_obs[skip_lines + n_hdr_obs + iobs].split()[1] + read_file_obs[skip_lines + n_hdr_obs + iobs].split()[2]
                date_obs.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) )
                obs[iobs-obs_start] = np.float( read_file_obs[skip_lines + n_hdr_obs + iobs][col_start:col_stop] )
            else:
                skip_lines = skip_lines + n_hdr_obs
                date_obs_temp = read_file_obs[skip_lines + n_hdr_obs + iobs].split()[0] + read_file_obs[skip_lines + n_hdr_obs + iobs].split()[1] + read_file_obs[skip_lines + n_hdr_obs + iobs].split()[2]
                date_obs.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
                obs[iobs-obs_start] = np.float( read_file_obs[skip_lines + n_hdr_obs + iobs].split()[col_start:col_stop] )

    return date_obs, obs

# format: swpc_3d
def read_pred(format, filename):
    print 'Format of the prediction file: ' + format
    print 'Filename: ' + filename
    bug = 0
    if format == 'swpc_3d':
        file_pred = open(filename)
        read_file_pred = file_pred.readlines()
        nday = 3
        yy_pred = read_file_pred[read_file_pred.index([s for s in read_file_pred if ":issued:" in s.lower()][0])].split()[1]
        # F10.7
        pred_f107 = np.zeros([nday])
        try:
            iline = read_file_pred.index([s for s in read_file_pred if "iv.  penticton" in s.lower()][0]) + 2 # f10.7 prediction line
            pred_f107[0] = read_file_pred[iline].split()[4].split('/')[0]
            pred_f107[1] = read_file_pred[iline].split()[4].split('/')[1]
            pred_f107[2] = read_file_pred[iline].split()[4].split('/')[2]
            date_pred_temp_start = read_file_pred[iline].split()[2].split('-')[0] + read_file_pred[iline].split()[1]
            if date_pred_temp_start == 'Jan01':
                yy_pred = str((int)(yy_pred) + 1)
            date_pred_start = datetime.strptime(yy_pred + date_pred_temp_start, "%Y%b%d")
        except IndexError:
            pred_f107 = np.zeros([nday]) - 999
            bug = 1
        # Ap
        pred_ap = np.zeros([nday])
        try:
            iline = read_file_pred.index([s for s in read_file_pred if "v.  geomagnetic" in s.lower()][0]) + 3 # ap prediction line
            try:
                pred_ap[0] = read_file_pred[iline].split()[5].split('-')[0].split('/')[1]
                pred_ap[1] = read_file_pred[iline].split()[5].split('-')[1].split('/')[1]
                pred_ap[2] = read_file_pred[iline].split()[5].split('-')[2].split('/')[1]
            except ValueError:
                pred_ap[:] = np.zeros([nday])-999
        except IndexError:
            pred_ap = np.zeros([nday]) - 999
    if format == 'swpc_27':
        nday = 27
        # F10.7
        pred_f107 = np.zeros([nday])
        file_pred = open(filename, 'rb')
        fileReader = PyPDF2.PdfFileReader(file_pred)

    date_pred = []
    if bug == 0:
        for iday in range(nday):
            date_pred.append(date_pred_start  + timedelta(days = iday))


    return date_pred, pred_f107, pred_ap
nday_pred = 3
year_start = 1997
year_stop = 2018
nyear = year_stop + 1 - year_start
leap_year = [1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020]
error_f107 = np.zeros([nyear, nday_pred])
error_ap = np.zeros([nyear, nday_pred])
iyear = -1
error_f107_all_day = []
error_ap_all_day = []
date_bug = []
date_ok = []
year_arr = np.zeros([nyear])
for year in range(year_start, year_stop + 1):
    iyear = iyear + 1
    # Observations
    format = 'swpc_f107' # swpc_f107, swpc_ap
    filename = '/Users/cbv/work/spaceWeather/swpcObs/' + str(year) + '_DSD.txt'
    date_start_str = str(year) + '-01-01'
    date_stop_str = str(year) + '-12-31'
    date_obs_f107, f107_obs = read_obs(format, filename, date_start_str, date_stop_str)
 
    format = 'swpc_ap' # swpc_f107, swpc_ap
    filename = '/Users/cbv/work/spaceWeather/swpcObs/' + str(year) + '_DGD.txt'
    date_obs_ap, ap_obs = read_obs(format, filename, date_start_str, date_stop_str)

    if year in leap_year:
        nday = 366
    else:
        nday = 365
    error_f107_temp = np.zeros([nday-nday_pred, nday_pred])
    error_ap_temp = np.zeros([nday-nday_pred, nday_pred])
    for iday in range(nday-nday_pred):
        print date_obs_f107[iday]
        date_now = str(datetime.strptime(str(year) + '-01-01', '%Y-%m-%d') + timedelta(days = iday))[:10].replace('-', '')
        # # Predictions
        format = 'swpc_3d' # swpc_3d
        filename = '/Users/cbv/work/spaceWeather/swpcPred/3day/' + str(year) + '_RSGA/' + date_now + 'RSGA.txt'
        try:
            date_pred, f107_pred, ap_pred = read_pred(format, filename)
        except IOError:
            date_bug.append(date_now)
            continue
        for ihor in range(nday_pred):
            error_f107_temp[iday, ihor] = f107_pred[ihor] - f107_obs[iday+1+ihor]
            error_ap_temp[iday, ihor] = ap_pred[ihor] - ap_obs[iday+1+ihor]
        date_ok.append(date_now)
    error_f107_all_day.append(error_f107_temp)
    error_ap_all_day.append(error_ap_temp)
    error_f107[iyear, :] = np.mean(error_f107_temp, axis = 0)
    error_ap[iyear, :] = np.mean(error_ap_temp, axis = 0)
    year_arr[iyear] = year


# Figures
