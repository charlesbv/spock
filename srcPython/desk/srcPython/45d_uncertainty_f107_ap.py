# copy of  new_uncertainty_f107_ap_bis.py on Aug 1. Changed to llook at 45 day predictions of 107 ap, not 3d. This was done fo rthe paper CA3 on CYGNSS collision acvoidance
# COmments below are still from new_uncertainty_f107_ap_bis.py so probably not all consistent with this new sscript 45d_uncertainty_f107_ap.py

# This script compares the historical predictions of F10.7 and Ap to observations. It derives the NRMS difference between both sets. This NRMS difference gives the uncertainty of the predictions of F10.7 and Ap. It is based on a comparison between observation and historical 3d predictions by SWPC since Jan 1, 2016.
# The directory that has the historical predicions of F10.7 and Ap is: /home/cbv/solar_wind_prediction/f107_ap/swpc/predictions/3d/historical/
# The directory that has the observations of F10.7 for 2016 is: /home/cbv/solar_wind_prediction/f107_ap/swpc/observations/f107
# The directory that has the observations of Ap for 2016 is: /home/cbv/solar_wind_prediction/f107_ap/swpc/observations/ap
# This script writes the results (sigma F10.7 and sigma Ap VS TIME) in a file called sigma_f107_ap.txt. This file is put in "../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_run + ".txt" and in "run_dir/input/density/density_NRLMSIS00e/sigma_f107_ap.txt" where run_dir is the first argument of the command line when running this script, and date_run is the date when running this script to create the files (format 'MMDDYY') (not that this date is not in "run_dir/input/density/density_NRLMSIS00e/sigma_f107_ap.txt". This is because the propagator look for a file called "sigma_f107_ap.txt" in the directory run_dir/input/density/density_NRLMSIS00e/ (file that does not include the date for the propagator))


import ipdb
import os
from os import listdir
import matplotlib.colors as colors
import matplotlib.cm as cmx
from os.path import isfile, join
from datetime import datetime, timedelta
import numpy as np
from matplotlib  import pyplot as plt
import sys



plt.ion()

#run_dir = sys.argv[1]

# Historical predictions
nb_days_pred = 45
path_hist_pred = 'f107_ap_data/solar_wind_prediction/predictions/historical/' # with last slash

date_start = datetime.strptime("20170120", "%Y%m%d")
date_stop = datetime.strptime("20180731", "%Y%m%d")#datetime.today()
n_pred = (date_stop - date_start).days

f107_pred_temp = []
ap_pred_temp = []
nb_na = 0
date_pred_all = []
date_pred_temp = []
for ipred in range(n_pred):
    filename_pred = str(date_start + timedelta(days = ipred))[:10].replace("-","") + ".txt"
    file_pred = open(path_hist_pred + filename_pred)
    read_file_pred = file_pred.readlines()
    f107_pred_temp_sub = []
    ap_pred_temp_sub = []
    if len(read_file_pred) !=0 :
        nline = 9 # 45 day prediciton, 5 days per line
        # skip header to ap data
        nheader = 0
        while read_file_pred[nheader].split()[0] != '45-DAY':
            nheader = nheader + 1
        nheader = nheader + 1
        date_pred_sub = []
        if (datetime.strptime(read_file_pred[nheader].split()[0], "%d%b%y") in date_pred_temp) == False: # sometimes two files in the data base are the same
            for iline in range(nline):
                for icol in range(5):
                    date_pred_sub.append( read_file_pred[iline+nheader].split()[2*icol] )
                    ap_pred_temp_sub.append( np.float( read_file_pred[iline+nheader].split()[2*icol + 1]) )

            # skip header to f107 data
            nheader = 0
            found_header_twice = 0
            while (read_file_pred[nheader].replace("\r", "").replace("\n","") != '45-DAY F10.7 CM FLUX FORECAST'):
                nheader = nheader + 1
            nheader = nheader + 1
            for iline in range(nline):
                for icol in range(5):
                    f107_pred_temp_sub.append( np.float( read_file_pred[iline+nheader].split()[2*icol + 1]) )


            date_pred_all.append(date_pred_sub)
            date_pred_temp.append(datetime.strptime(date_pred_sub[0], "%d%b%y"))
            f107_pred_temp.append(f107_pred_temp_sub)
            ap_pred_temp.append(ap_pred_temp_sub)
date_start_pred = date_start
# F10.7 observations
date_today = datetime.today()
filename_f107_obs = "f107_ap_data/solar_wind_prediction/observations/f107/2017_concatenated_2018Q1_concatenated_2018Q2_concatenated_2018Q3_DSD.txt"
# os.system("rm -f " + filename_f107_obs)
# os.system("cat " + filename_f107_obs1 + " " + filename_f107_obs2 + " " + filename_f107_obs3 + " " + filename_f107_obs4 + " >> " + filename_f107_obs)
file_f107_obs = open(filename_f107_obs)
read_file_f107_obs = file_f107_obs.readlines()
n_hdr_f107_obs = 0
while (read_file_f107_obs[n_hdr_f107_obs].split(' ')[0][0:4] != '#---'):
    n_hdr_f107_obs = n_hdr_f107_obs + 1
n_hdr_f107_obs = n_hdr_f107_obs + 1
n_f107_obs = len(read_file_f107_obs) - n_hdr_f107_obs
f107_obs_temp = np.zeros([n_pred])
## Skip observations prior to the start date of the historical predictions
iobs = 0
date_obs_temp = read_file_f107_obs[n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[2]
date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
while ( date_obs != date_start_pred ):
    iobs = iobs + 1
    date_obs_temp = read_file_f107_obs[n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[n_hdr_f107_obs + iobs].split()[2]
    date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
obs_start = iobs
date_obs_f107_temp = [] # !!!!!! if everything goes well then: date_obs_f107_temp = date_pred. Otherwise, there is a problem somewhere so debug!
skip_lines = 0
for iobs in range( obs_start, obs_start + n_pred):
    if ( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] != ":Product:" ):
        date_obs_temp = read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[2]
        date_obs_f107_temp.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        f107_obs_temp[iobs-obs_start] = np.float( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[3] )
    else:
        skip_lines = skip_lines + n_hdr_f107_obs
        date_obs_temp = read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[2]
        date_obs_f107_temp.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        f107_obs_temp[iobs-obs_start] = np.float( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[3] )
        
# # Ap observations
filename_ap_obs = "f107_ap_data/solar_wind_prediction/observations/ap/2017_concatenated_2018Q1_concatenated_2018Q2_concatenated_2018Q3_DGD.txt"
# os.system("rm -f " + filename_ap_obs)
# os.system("cat " + filename_ap_obs1 + " " + filename_ap_obs2 + " " + filename_ap_obs3 + " " + filename_ap_obs4 + " >> " + filename_ap_obs)
file_ap_obs = open(filename_ap_obs)
read_file_ap_obs = file_ap_obs.readlines()
n_hdr_ap_obs = 0
while (read_file_ap_obs[n_hdr_ap_obs].split('   ')[0][0:7] != '#  Date'):
    n_hdr_ap_obs = n_hdr_ap_obs + 1
n_hdr_ap_obs = n_hdr_ap_obs + 1
n_ap_obs = len(read_file_ap_obs) - n_hdr_ap_obs
ap_obs_temp = np.zeros([n_pred])
## Skip observations prior to the start date of the historical predictions
iobs = 0
date_obs_temp = read_file_ap_obs[n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[2]
date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
while ( date_obs != date_start_pred ):
    iobs = iobs + 1
    date_obs_temp = read_file_ap_obs[n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[n_hdr_ap_obs + iobs].split()[2]
    date_obs = datetime.strptime( date_obs_temp, "%Y%m%d" )
obs_start = iobs
date_obs_ap_temp = [] # !!!!!! if everything goes well then: date_obs_ap_temp = date_pred. Otherwise, there is a problem somewhere so debug!
skip_lines = 0
for iobs in range( obs_start, obs_start + n_pred):
    if ( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] != ":Product:" ):
        date_obs_temp = read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[2]
        date_obs_ap_temp.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        ap_obs_temp[iobs-obs_start] = np.float( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split('   ')[3].split()[0] )
    else:
        skip_lines = skip_lines + n_hdr_ap_obs
        date_obs_temp = read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[2]
        date_obs_ap_temp.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        ap_obs_temp[iobs-obs_start] = np.float( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split('   ')[3].split()[0] )


# make sure that date_pred is the same as date_obs_ap and date_obs_f107_temp for all pred.
if date_obs_ap_temp != date_obs_f107_temp:
    print "***! Weird, f107 and ap observations files should have the same dates. The program will stop. !***"; raise Exception
## first make sure all elements of ap and f107 obs are in pred
n_pred = len(date_pred_temp)
ipred_count = -1
date_f107_ap_removed = []
ap_obs = np.zeros([n_pred])
date_obs_ap = []
f107_obs = np.zeros([n_pred])
date_obs_f107 = []

for ipred in range(len(date_obs_ap_temp)):
    if (date_obs_ap_temp[ipred] in date_pred_temp) == True:
        ipred_count = ipred_count + 1
        ap_obs[ipred_count] = ap_obs_temp[ipred]
        date_obs_ap.append(date_obs_ap_temp[ipred])
        f107_obs[ipred_count] = f107_obs_temp[ipred]
        date_obs_f107.append(date_obs_f107_temp[ipred])
    else:
        date_f107_ap_removed.append([date_obs_ap_temp[ipred],f107_obs_temp[ipred],  ap_obs_temp[ipred]])
        

## then make sure all pred are in obs
date_pred = []
date_pred_removed = []
f107_pred = []
ap_pred = []
for ipred in range(len(date_pred_temp)):
    if ( date_pred_temp[ipred] in date_obs_f107 ) == True:
        date_pred.append(date_pred_temp[ipred])
        f107_pred.append(f107_pred_temp[ipred])
        ap_pred.append(ap_pred_temp[ipred])
    else:
        date_pred_removed.append( [date_pred_temp[ipred], f107_pred_temp[ipred], ap_pred_temp[ipred]])
    
f107_pred_arr = np.array(f107_pred)
ap_pred_arr = np.array(ap_pred)
n_pred = len(date_pred)
# RMS difference between observations and predictions
nbins = 10
## F10.7
new_nb_days_pred = nb_days_pred
rms_error_f107 = np.zeros( new_nb_days_pred )
f107_bin_pred_array = np.zeros([new_nb_days_pred, nbins])
percentage_in_f107_bin_pred_array = np.zeros([new_nb_days_pred, nbins])
quartile_array = np.arange(10, 91, 10)
nquartile = len(quartile_array)
f107_quartile = np.zeros([new_nb_days_pred, nquartile])
difference_pred_obs_f107_array = np.zeros([new_nb_days_pred, n_pred])
print 'F107'
for iforecast in range(1, new_nb_days_pred + 1):
    difference_pred_obs_f107 = np.zeros([n_pred - iforecast ])
    for ipred in range(n_pred-iforecast):
        difference_pred_obs_f107[ipred]= f107_pred_arr[ipred, iforecast-1] - f107_obs[ipred + iforecast-1]
    if iforecast < 10:
        filename_out = 'f107_ap_forecast_error_day_' + str(iforecast) + '.txt'
        file_out = open(filename_out, "w")
        print >> file_out, "bin\nf107_dist\nap_dist"
        fig,ax = plt.subplots()
        fig.suptitle('f107 - ' + str(iforecast))
        hist_f107 = ax.hist(difference_pred_obs_f107, bins = 7, range = (-15,15))[0] / len(difference_pred_obs_f107)
        print >> file_out,  ' '.join(map(str, ax.hist(difference_pred_obs_f107, bins = 7, range = (-15,15))[1]))
        print >> file_out, ' '.join(map(str,hist_f107))
#         print iforecast
#         print hist_f107
#         print ''


    difference_pred_obs_ap = np.zeros([n_pred - iforecast ])
    for ipred in range(n_pred-iforecast):
        difference_pred_obs_ap[ipred]= ap_pred_arr[ipred, iforecast-1] - ap_obs[ipred + iforecast-1]
    if iforecast < 10:
        fig,ax = plt.subplots()
        fig.suptitle('ap - ' + str(iforecast))
        hist_ap = ax.hist(difference_pred_obs_ap, bins = 7, range = (-15,15))[0] / len(difference_pred_obs_ap)
        print >> file_out, ' '.join(map(str,hist_ap))
#         print iforecast
#         print hist_ap
#         print ''

    file_out.close()
