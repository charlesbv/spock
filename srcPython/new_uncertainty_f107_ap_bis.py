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

# This script compares the historical predictions of F10.7 and Ap to observations. It derives the NRMS difference between both sets. This NRMS difference gives the uncertainty of the predictions of F10.7 and Ap. It is based on a comparison between observation and historical 3d predictions by SWPC since Jan 1, 2016.
# The directory that has the historical predicions of F10.7 and Ap is: /home/cbv/solar_wind_prediction/f107_ap/swpc/predictions/3d/historical/
# The directory that has the observations of F10.7 for 2016 is: /home/cbv/solar_wind_prediction/f107_ap/swpc/observations/f107
# The directory that has the observations of Ap for 2016 is: /home/cbv/solar_wind_prediction/f107_ap/swpc/observations/ap
# This script writes the results (sigma F10.7 and sigma Ap VS TIME) in a file called sigma_f107_ap.txt. This file is put in "../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_run + ".txt" and in "run_dir/input/density/density_NRLMSIS00e/sigma_f107_ap.txt" where run_dir is the first argument of the command line when running this script, and date_run is the date when running this script to create the files (format 'MMDDYY') (not that this date is not in "run_dir/input/density/density_NRLMSIS00e/sigma_f107_ap.txt". This is because the propagator look for a file called "sigma_f107_ap.txt" in the directory run_dir/input/density/density_NRLMSIS00e/ (file that does not include the date for the propagator))



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
nb_days_pred = 3 
path_hist_pred = './solar_wind_prediction/f107_ap/swpc/predictions/3d/historical/'

date_start = datetime.strptime("20160101", "%Y%m%d")
date_stop = datetime.strptime("20161120", "%Y%m%d")#datetime.today()
n_pred = (date_stop - date_start).days

f107_pred = np.zeros([n_pred, nb_days_pred])
ap_pred = np.zeros([n_pred, nb_days_pred])
nb_na = 0
for ipred in range(n_pred):
    filename_pred = str(date_start + timedelta(days = ipred))[:10].replace("-","") + ".txt"
    file_pred = open("./solar_wind_prediction/f107_ap/swpc/predictions/3d/historical/" + filename_pred)
    read_file_pred = file_pred.readlines()
    nline = len(read_file_pred)
    for iline in range(nline):
        if read_file_pred[iline].split(' ')[0] == 'IV.':
            f107_pred[ipred, 0] = np.float(read_file_pred[iline+2].split('/')[0].split(' ')[-1])
            f107_pred[ipred, 1] = np.float(read_file_pred[iline+2].split('/')[1])
            f107_pred[ipred, 2] = np.float(read_file_pred[iline+2].split('/')[2].replace("\n",""))

        if read_file_pred[iline].split(' ')[0] == 'V.':   
            if read_file_pred[iline+3].split('/')[2].split('-')[0].replace(" ", "") == 'NA': # if no prediction for ap then take pred of afr
                ap_pred[ipred, 0] = np.float( read_file_pred[iline+3].split('/')[1].split(' ')[-1] )
                ap_pred[ipred, 1] = np.float( read_file_pred[iline+3].split('/')[2].split('-')[1] )
                ap_pred[ipred, 2] = np.float( read_file_pred[iline+3].split('/')[3].split('-')[1] )
                nb_na = nb_na + 1
            else:                
                ap_pred[ipred, 0] = np.float( read_file_pred[iline+3].split('/')[2].split('-')[0] )
                ap_pred[ipred, 1] = np.float( read_file_pred[iline+3].split('/')[3].split('-')[0] )
                ap_pred[ipred, 2] = np.float( read_file_pred[iline+3].split('/')[4].replace("\n","") )
            
date_start_pred = date_start
# F10.7 observations
date_today = datetime.today()
filename_f107_obs1 = "./solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q1.txt"
filename_f107_obs2 = "./solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q2.txt"
filename_f107_obs3 = "./solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q3.txt"
filename_f107_obs4 = "./solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q4.txt"
filename_f107_obs = "./solar_wind_prediction/f107_ap/swpc/observations/f107/f107_2016Q1_concatenated_with_2016Q2_concatenated_with_2016Q3_concatenated_with_2016Q4.txt"
# os.system("rm -f " + filename_f107_obs)
# os.system("cat " + filename_f107_obs1 + " " + filename_f107_obs2 + " " + filename_f107_obs3 + " " + filename_f107_obs4 + " >> " + filename_f107_obs)
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
        skip_lines = skip_lines + n_hdr_f107_obs
        date_obs_temp = read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[0] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[1] + read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[2]
        date_obs_f107.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        f107_obs[iobs-obs_start] = np.float( read_file_f107_obs[skip_lines + n_hdr_f107_obs + iobs].split()[3] )
        

# # Ap observations
filename_ap_obs1 = "./solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q1.txt"
filename_ap_obs2 = "./solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q2.txt"
filename_ap_obs3 = "./solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q3.txt"
filename_ap_obs4 = "./solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q4.txt"
filename_ap_obs = "./solar_wind_prediction/f107_ap/swpc/observations/ap/ap_2016Q1_concatenated_with_2016Q2_concatenated_with_2016Q3_concatenated_with_2016Q4.txt"
# os.system("rm -f " + filename_ap_obs)
# os.system("cat " + filename_ap_obs1 + " " + filename_ap_obs2 + " " + filename_ap_obs3 + " " + filename_ap_obs4 + " >> " + filename_ap_obs)
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
        ap_obs[iobs-obs_start] = np.float( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split('   ')[3].split()[0] )
    else:
        skip_lines = skip_lines + n_hdr_ap_obs
        date_obs_temp = read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[0] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[1] + read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split()[2]
        date_obs_ap.append( datetime.strptime( date_obs_temp, "%Y%m%d" ) ) 
        ap_obs[iobs-obs_start] = np.float( read_file_ap_obs[skip_lines + n_hdr_ap_obs + iobs].split('   ')[3].split()[0] )


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
for iforecast in range(1, new_nb_days_pred + 1):
    difference_pred_obs_f107 = np.zeros([n_pred - iforecast ])
    for ipred in range(n_pred-iforecast):
        difference_pred_obs_f107[ipred] = f107_pred[ipred, iforecast-1] - f107_obs[ipred + iforecast]
    rms_error_f107[iforecast-1] = np.sqrt(  np.mean( ( difference_pred_obs_f107 )**2 ) ) 
    for ibin in range(nbins):
        f107_bin_pred_array[iforecast-1, ibin] = ( np.histogram(difference_pred_obs_f107, bins = nbins, density = 1)[1][ibin+1] + np.histogram(difference_pred_obs_f107, bins = nbins, density = 1)[1][ibin]  ) / 2.
    if iforecast == 2:
        fig,ax = plt.subplots()
        hist_f107 = ax.hist(difference_pred_obs_f107, bins = 7, range = (-15,15))[0] / len(difference_pred_obs_f107)
        print ''
        print hist_f107[1]
        print hist_f107[0]/len(difference_pred_obs_f107)*100
    #percentage_in_f107_bin_pred_array[iforecast-1, :] = np.histogram(difference_pred_obs_f107, bins = nbins, density = 1)[0] # if density is True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1. 
    #Note on preivous line:  the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; it is not a probability mass function.
    for iquartile in range(nquartile):
        f107_quartile[iforecast-1, iquartile] = np.percentile(difference_pred_obs_f107, quartile_array[iquartile])
    difference_pred_obs_f107_array[iforecast-1, 0:n_pred-iforecast] = difference_pred_obs_f107[0:n_pred-iforecast]
#     if iforecast == 1:
#         hist_f107_ave =  hist_f107[0]
#     else:
#         hist_f107_ave = hist_f107_ave + hist_f107[0]
    #print np.percentile(difference_pred_obs_f107,75) - np.percentile(difference_pred_obs_f107,25)

#hist_f107_ave = hist_f107_ave / new_nb_days_pred
#print hist_f107_ave    
# GENERATE DISCTINCT COLORS
# fig,ax = plt.subplots()
# for ipred in range(n_pred-new_nb_days_pred):
#     ax.plot(difference_pred_obs_f107_array[:, ipred])

NCURVES = nquartile
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# for iquartile in range(nquartile):
#     colorVal = scalarMap.to_rgba(values[iquartile])
#     ax.plot(f107_quartile[:, iquartile], color = colorVal)
# raise Exception
## Ap
rms_error_ap = np.zeros( new_nb_days_pred )
ap_bin_pred_array = np.zeros([new_nb_days_pred, nbins])
percentage_in_ap_bin_pred_array = np.zeros([new_nb_days_pred, nbins])
nquartile = len(quartile_array)
ap_quartile = np.zeros([new_nb_days_pred, nquartile])


for iforecast in range(1, new_nb_days_pred + 1):
    difference_pred_obs_ap = np.zeros([n_pred - iforecast ])
    for ipred in range(n_pred-iforecast):
        difference_pred_obs_ap[ipred] = ap_pred[ipred, iforecast-1] - ap_obs[ipred + iforecast]
    rms_error_ap[iforecast-1] = np.sqrt(  np.mean( ( difference_pred_obs_ap )**2 ) ) 
    for ibin in range(nbins):
        ap_bin_pred_array[iforecast-1, ibin] = ( np.histogram(difference_pred_obs_ap, bins = nbins, density = 1)[1][ibin+1] + np.histogram(difference_pred_obs_ap, bins = nbins, density = 1)[1][ibin]  ) / 2.
    if iforecast == 2:
        fig,ax = plt.subplots()
        hist_ap = ax.hist(difference_pred_obs_ap, bins = 7, range = (-15,15))[0] / len(difference_pred_obs_ap)
        print hist_ap[1]
        print hist_ap[0] / len(difference_pred_obs_ap) * 100

    percentage_in_ap_bin_pred_array[iforecast-1, :] = np.histogram(difference_pred_obs_ap, bins = nbins, density = 1)[0] # if density is True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1. 
    # Note on previous line:  the sum of the histogram values will not be equal to 1 unless bins of unity width are chosen; it is not a probability mass function.
    for iquartile in range(nquartile):
        ap_quartile[iforecast-1, iquartile] = np.percentile(difference_pred_obs_ap, quartile_array[iquartile])
    #print np.percentile(difference_pred_obs_ap,75) - np.percentile(difference_pred_obs_ap,25)
# for iquartile in range(nquartile):
#     colorVal = scalarMap.to_rgba(values[iquartile])
#     ax.plot(ap_quartile[:, iquartile], color = colorVal)


# # Write uncertainties in a file
raise Exception
date_today_str = datetime.strftime( datetime.today(), "%m%d%y")
file_out = open("../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_today_str + ".txt", "w+")
print >> file_out, "#1-sigma uncertainty on F10.7 and Ap at a function of the forecast time\n#Time(day) sigma_F10.7 sigma_Ap"
# # !!!!!!!!!!! TEMPORARY ERASE
# rms_error_ap = np.arange(1,nb_days_pred+2)*1.5
# rms_error_f107 = np.arange(1,nb_days_pred+2)*4
# # !!!!!!!!!!! END OF TEMPORARY ERASE
new_nb_days_pred = 45
for ipred in range(1,new_nb_days_pred+1):
    if ipred <= 3:
        print >> file_out, ipred, (int)(round(rms_error_f107[ipred-1])), (int)(round(rms_error_ap[ipred-1]))
    else:
        print >> file_out, ipred, 1, 1

file_out.close()
raise Exception
# Copy this file to the density_NRLMSIS00e directory of the run_dir chosen in the command line
raise Exception
os.system("cp ../Examples_more/Cygnss/density_NRLMSIS00e/sigma_f107_ap_" + date_today_str + ".txt ../" + run_dir + "/input/density/density_NRLMSIS00e/sigma_f107_ap.txt" )
