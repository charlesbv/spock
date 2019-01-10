# This script compares the historical predictions of F10.7 and Ap to observations. It derives the NRMS difference between both sets. This NRMS difference gives the uncertainty of the predictions of F10.7 and Ap.
# The directory that has the historical predicions of F10.7 and Ap is: /home/cbv/solar_wind_prediction/f107_ap/swpc/predictions/historical
# The directory that has the observations of F10.7 for 2016 (until Sep 21) is: /home/cbv/solar_wind_prediction/f107_ap/swpc/observations/f107
# The directory that has the observations of Ap for 2016 (until Sep 21) is: /home/cbv/solar_wind_prediction/f107_ap/swpc/observations/ap
# This script writes the results (sigma F10.7 and sigma Ap VS TIME) in a file called sigma_f107_ap.txt. This file is put in "../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_run + ".txt" and in "run_dir/input/density/density_NRLMSIS00e/sigma_f107_ap.txt" where run_dir is the first argument of the command line when running this script, and date_run is the date when running this script to create the files (format 'MMDDYY') (not that this date is not in "run_dir/input/density/density_NRLMSIS00e/sigma_f107_ap.txt". This is because the propagator look for a file called "sigma_f107_ap.txt" in the directory run_dir/input/density/density_NRLMSIS00e/ (file that does not include the date for the propagator))
# Note: since the historical forecast data that is available only starts on Sep 1st 2016, we actually don't need the observations of Ap and F10.7 for dates prior to Sep 1st 2016.

# ASSUMPTIONS: 
# - the predictions are 45-day
# - the date of the most recent prediction file is older than or equal to the date as the most recent observation (last line of ap_2016Q3.txt and ap_2016Q3.txt)

import os
from os import listdir
from os.path import isfile, join
from datetime import datetime, timedelta
import numpy as np
from matplotlib  import pyplot as plt
import sys

plt.ion()

#run_dir = sys.argv[1]

# Historical predictions
nb_days_pred = 45 # !!! if you want to change this number, you probably need to change the statements that read the predictions files
path_hist_pred = '/home/cbv/solar_wind_prediction/f107_ap/swpc/predictions/historical/'
list_filename_pred_temp = [path_hist_pred + f for f in listdir(path_hist_pred) if isfile(join(path_hist_pred, f))]
n_pred = len(list_filename_pred_temp)
date_pred_not_sort = []
for ipred in range(n_pred):
    date_pred_temp = "16" + list_filename_pred_temp[ipred].split('/')[-1].split('.')[0]
    date_pred_not_sort.append( datetime.strptime(date_pred_temp, "%y%m%d") ) 

date_pred = sorted(date_pred_not_sort)
date_start_pred = date_pred[0]
list_filename_pred = np.array(list_filename_pred_temp)
list_filename_pred = list_filename_pred[sorted(range(len(date_pred_not_sort)),key=lambda x:date_pred_not_sort[x])]

f107_pred = np.zeros([n_pred, nb_days_pred])
ap_pred = np.zeros([n_pred, nb_days_pred])
for ipred in range(n_pred):
    filename_pred = list_filename_pred[ipred]
    file_pred = open(filename_pred)
    read_file_pred = file_pred.readlines()
    n_header_pred = 0
    while (read_file_pred[n_header_pred].split(' ')[0].replace("\n", "") != '45-DAY'):
        n_header_pred = n_header_pred + 1
    n_header_pred = n_header_pred + 1
    ncol = 5
    nline = 9
    # Ap predictions
    for iline in range(nline):# 9 lines of predictions per file
        for icol in range(ncol):# 5 predictions per line
            ap_pred[ipred, iline*ncol + icol] = np.float( read_file_pred[n_header_pred + iline].split()[2*icol+1] )
    # F10.7 predictions
    for iline in range(nline):# 9 lines of predictions per file
        for icol in range(ncol):# 5 predictions per line
            f107_pred[ipred, iline*ncol + icol] = np.float( read_file_pred[n_header_pred + nline + 1 + iline].split()[2*icol+1] ) # nline to skip the Ap block, +1 to skip the header of f10.7

# F10.7 observations
date_today = datetime.today()
filename_f107_obs1 = "/home/cbv/solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q3.txt"
filename_f107_obs2 = "/home/cbv/solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q4.txt"
filename_f107_obs = "/home/cbv/solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q3_concatenated_with_2016Q4.txt"
os.system("rm -f " + filename_f107_obs)
os.system("cat " + filename_f107_obs1 + " " + filename_f107_obs2 + " >> " + filename_f107_obs)
file_f107_obs = open(filename_f107_obs)
read_file_f107_obs = file_f107_obs.readlines()
n_hdr_f107_obs = 0
while (read_file_f107_obs[n_hdr_f107_obs].split(' ')[0][0:4] != '#---'):
    n_hdr_f107_obs = n_hdr_f107_obs + 1
n_hdr_f107_obs = n_hdr_f107_obs + 1
n_f107_obs = len(read_file_f107_obs) - n_hdr_f107_obs
f107_obs = np.zeros([n_pred])
## Skip observations prior to the start date of the historical predictions
iobs = 0
date_obs_temp = read_file_f107_obs[n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[2]
date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
while ( date_obs != date_start_pred ):
    iobs = iobs + 1
    date_obs_temp = read_file_f107_obs[n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[2]
    date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
obs_start = iobs
date_obs_f107 = [] # !!!!!! if everything goes well then: date_obs_f107 = date_pred. Otherwise, there is a problem somewhere so debug!
skip_lines = 0
for iobs in range( obs_start, obs_start + n_pred):
    if ( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] != ":Product:" ):
        date_obs_temp = read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[2]
        date_obs_f107.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        f107_obs[iobs-obs_start] = np.float( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[3] )
    else:
        skip_lines = n_hdr_f107_obs
        date_obs_temp = read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[2]
        date_obs_f107.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        f107_obs[iobs-obs_start] = np.float( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[3] )
        

# # Ap observations
filename_ap_obs1 = "/home/cbv/solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q3.txt"
filename_ap_obs2 = "/home/cbv/solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q4.txt"
filename_ap_obs = "/home/cbv/solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q3_concatenated_with_2016Q4.txt"
os.system("rm -f " + filename_ap_obs)
os.system("cat " + filename_ap_obs1 + " " + filename_ap_obs2 + " >> " + filename_ap_obs)
file_ap_obs = open(filename_ap_obs)
read_file_ap_obs = file_ap_obs.readlines()
n_hdr_ap_obs = 0
while (read_file_ap_obs[n_hdr_ap_obs].split('   ')[0][0:7] != '#  Date'):
    n_hdr_ap_obs = n_hdr_ap_obs + 1
n_hdr_ap_obs = n_hdr_ap_obs + 1
n_ap_obs = len(read_file_ap_obs) - n_hdr_ap_obs
ap_obs = np.zeros([n_pred])
## Skip observations prior to the start date of the historical predictions
iobs = 0
date_obs_temp = read_file_ap_obs[n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[2]
date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
while ( date_obs != date_start_pred ):
    iobs = iobs + 1
    date_obs_temp = read_file_ap_obs[n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[2]
    date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
obs_start = iobs
date_obs_ap = [] # !!!!!! if everything goes well then: date_obs_ap = date_pred. Otherwise, there is a problem somewhere so debug!
skip_lines = 0
for iobs in range( obs_start, obs_start + n_pred):
    if ( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] != ":Product:" ):
        date_obs_temp = read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[2]
        date_obs_ap.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        ap_obs[iobs-obs_start] = np.float( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[21] )
    else:
        skip_lines = n_hdr_ap_obs
        date_obs_temp = read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[2]
        date_obs_ap.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        ap_obs[iobs-obs_start] = np.float( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[21] )

# RMS difference between observations and predictions
## F10.7
new_nb_days_pred = np.min([nb_days_pred, n_pred])-1
rms_error_f107 = np.zeros( new_nb_days_pred )
for iforecast in range(1, new_nb_days_pred + 1):
    difference_pred_obs = np.zeros([n_pred - iforecast ])
    for ipred in range(n_pred-iforecast):
        difference_pred_obs[ipred] = f107_pred[ipred, iforecast-1] - f107_obs[ipred + iforecast]
    rms_error_f107[iforecast-1] = np.sqrt(  np.mean( ( difference_pred_obs )**2 ) ) 

## Ap
rms_error_ap = np.zeros( new_nb_days_pred )
for iforecast in range(1, new_nb_days_pred + 1):
    difference_pred_obs = np.zeros([n_pred - iforecast ])
    for ipred in range(n_pred-iforecast):
        difference_pred_obs[ipred] = ap_pred[ipred, iforecast-1] - ap_obs[ipred + iforecast]
    rms_error_ap[iforecast-1] = np.sqrt(  np.mean( ( difference_pred_obs )**2 ) ) 

# # Write uncertainties in a file
date_today_str = datetime.strftime( datetime.today(), "%m%d%y")
file_out = open("../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_today_str + ".txt", "w+")
print >> file_out, "#1-sigma uncertainty on F10.7 and Ap at a function of the forecast time\n#Time(day) sigma_F10.7 sigma_Ap"
# # !!!!!!!!!!! TEMPORARY ERASE
# rms_error_ap = np.arange(1,nb_days_pred+2)*1.5
# rms_error_f107 = np.arange(1,nb_days_pred+2)*4
# # !!!!!!!!!!! END OF TEMPORARY ERASE

for ipred in range(1,new_nb_days_pred+1):
    print >> file_out, ipred, rms_error_f107[ipred-1], rms_error_ap[ipred-1]

file_out.close()

# Copy this file to the density_NRLMSIS00e directory of the run_dir chosen in the command line
raise Exception
os.system("cp ../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_today_str + ".txt ../" + run_dir + "/input/density/density_NRLMSIS00e/sigma_f107_ap.txt" )
