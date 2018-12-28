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
from matplotlib import pyplot as plt
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from datetime import datetime, timedelta
import numpy as np
import os

plt.ion()

day_stop_analysis = "20160815"
day_stop_analysis = datetime.strptime(day_stop_analysis, "%Y%m%d")

# Predictions by SWPC
day_start_analysis = "20160707"
day_start_analysis = datetime.strptime(day_start_analysis, "%Y%m%d")
nb_days_analysis = 60
day_now_filename = day_start_analysis
ap_pred_swpc = [] 
date_pred_swpc = []
for iday in range(nb_days_analysis):
    if iday != 15:
        if ( day_now_filename != day_stop_analysis ):
            day_now_filename_str =  datetime.strftime(day_now_filename, "%m%d")
            filename_wget = "ftp://ftp.swpc.noaa.gov/pub/forecasts/geomag_forecast/" + day_now_filename_str + "geomag_forecast.txt"
        #    os.system("wget " + filename_wget + " -O kp_ap_files/temp_ap_pred/" + day_now_filename_str + ".txt")
            file_pred_in = open("kp_ap_files/temp_ap_pred/" + day_now_filename_str + ".txt")
            read_file_pred_in = file_pred_in.readlines()
            line_pred = 7
            date_pred_swpc.append( day_now_filename )
            ap_pred_swpc_sub = []
            ap_pred_swpc_sub.append(np.float(read_file_pred_in[line_pred-1].split()[4])) # first index is estimated
            ap_pred_swpc_sub.append(np.float(read_file_pred_in[line_pred].split()[5].split('-')[0])) # three last index are predicted
            ap_pred_swpc_sub.append(np.float(read_file_pred_in[line_pred].split()[5].split('-')[1]))
            ap_pred_swpc_sub.append(np.float(read_file_pred_in[line_pred].split()[5].split('-')[2]))
            ap_pred_swpc.append(ap_pred_swpc_sub)
        else:
            break
    day_now_filename = day_now_filename + timedelta(days = 1)

ap_pred_swpc = np.array(ap_pred_swpc)
n_swpc = len(date_pred_swpc)

# Read OMNIWEB
# Omniweb
filename_omni = 'kp_ap_files/all_info_2_months.txt'
file_omni = open(filename_omni, "r")
read_file_omni = file_omni.readlines()
n_header_omni = 0
n_elts_omni = len(read_file_omni) - n_header_omni
kp_omni = []; ap_omni = []; 
bx = []; by_gse = []; bz_gse = []; 
by_gsm = []; bz_gsm = []; 
t = []; n = []; v = []; 
f107 = [];
date_omni = []
for i in range(n_elts_omni):
    date_omni_temp = read_file_omni[i + n_header_omni].split()[0]+ '-' + read_file_omni[i + n_header_omni].split()[1]+ '-' + read_file_omni[i + n_header_omni].split()[2]
    if datetime.strptime(date_omni_temp, "%Y-%j-%H")  != day_stop_analysis:
        if datetime.strptime(date_omni_temp, "%Y-%j-%H")  != datetime.strptime("20160722", "%Y%m%d"):
            date_omni.append( datetime.strptime(date_omni_temp, "%Y-%j-%H") )
            by_gse.append( np.float(read_file_omni[i + n_header_omni].split()[5]) )
            bz_gse.append( np.float(read_file_omni[i + n_header_omni].split()[6]) )
            by_gsm.append( np.float(read_file_omni[i + n_header_omni].split()[7]) )
            bz_gsm.append( np.float(read_file_omni[i + n_header_omni].split()[8]) )
            t.append( np.float(read_file_omni[i + n_header_omni].split()[9]) )
            n.append( np.float(read_file_omni[i + n_header_omni].split()[10]) )
            v.append( np.float(read_file_omni[i + n_header_omni].split()[11]) ) 
            kp_omni.append( np.float(read_file_omni[i + n_header_omni].split()[12]) / 10. ) # /10 because omniweb gives kp * 10
            ap_omni.append( np.float(read_file_omni[i + n_header_omni].split()[13]) )
            f107.append( np.float(read_file_omni[i + n_header_omni].split()[14]) )
            if ((bz_gsm[-1] == 0) & (by_gsm[-1] == 0)):
                print 'GSM', date_omni
                bz_gsm[-1] = 0.1 # to avoid 0/0, we set bz to this so that theta = 0
            if ((bz_gse[-1] == 0) & (by_gse[-1] == 0)):
                print 'GSE', date_omni
                bz_gse[-1] = 0.1 # to avoid 0/0, we set bz to this so that theta = 0
    else:
        break

n_omni = len(date_omni)

kp_omni = np.array(kp_omni); ap_omni = np.array(ap_omni); 
bx = np.array(bx); by_gse = np.array(by_gse); bz_gse = np.array(bz_gse); 
by_gsm = np.array(by_gsm); bz_gsm = np.array(bz_gsm); 
t = np.array(t); n = np.array(n); v = np.array(v); 
f107 = np.array(f107);


# 3 hourly average of these quantities


# Newell 2008  GSM
bt_gsm = np.sqrt(by_gsm**2 + bz_gsm**2)
theta_gsm = np.arctan(by_gsm / bz_gsm)
f1_gsm = v**(4./3) * bt_gsm**(2./3) * ((np.sin(theta_gsm/2.))**8)**(1./3)
f2_gsm = np.sqrt(n)*v**2
kp_newell_2008_gsm = 0.05 + 2.244 * 10**(-4) * f1_gsm + 2.844*10**(-6) * f2_gsm

## Conversion to Ap
ap_newell_2008_gsm = np.zeros([n_omni])
ap_array_convert = np.array([0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,80,94,111,132,154,179,207,236,300,400])
kp_array_convert = np.arange(0,9.33,1/3.)
for iday in range(n_omni):
    j = 0
    if kp_newell_2008_gsm[iday] < kp_array_convert[j]:
        ap_newell_2008_gsm[iday] = ap_array_convert[j]
    while (j + 1 < len(kp_array_convert)):
        if ( ( kp_newell_2008_gsm[iday] < kp_array_convert[j+1] ) & ( kp_newell_2008_gsm[iday] >= kp_array_convert[j] ) ):
            ap_newell_2008_gsm[iday] = ap_array_convert[j]
        j = j + 1
    if  kp_newell_2008_gsm[iday] > kp_array_convert[j]:
        ap_newell_2008_gsm[iday] = ap_array_convert[j]


# Newell 2008  GSE
bt_gse = np.sqrt(by_gse**2 + bz_gse**2)
theta_gse = np.arctan(by_gse / bz_gse)
f1_gse = v**(4./3) * bt_gse**(2./3) * ((np.sin(theta_gse/2.))**8)**(1./3)
f2_gse = np.sqrt(n)*v**2
kp_newell_2008_gse = 0.05 + 2.244 * 10**(-4) * f1_gse + 2.844*10**(-6) * f2_gse

## Conversion to Ap
ap_newell_2008_gse = np.zeros([n_omni])
ap_array_convert = np.array([0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,80,94,111,132,154,179,207,236,300,400])
kp_array_convert = np.arange(0,9.33,1/3.)
for iday in range(n_omni):
    j = 0
    if kp_newell_2008_gse[iday] < kp_array_convert[j]:
        ap_newell_2008_gse[iday] = ap_array_convert[j]
    while (j + 1 < len(kp_array_convert)):
        if ( ( kp_newell_2008_gse[iday] < kp_array_convert[j+1] ) & ( kp_newell_2008_gse[iday] >= kp_array_convert[j] ) ):
            ap_newell_2008_gse[iday] = ap_array_convert[j]
        j = j + 1
    if  kp_newell_2008_gse[iday] > kp_array_convert[j]:
        ap_newell_2008_gse[iday] = ap_array_convert[j]


# # Newell lin reg
time_lin_reg = 2
nb_days_pred = 4
rms_error_newell_lin_reg = []
for inow in range(time_lin_reg, n_omni-4):
    kp_newell_lin_reg = np.zeros([nb_days_pred])
    f1_for_lin_reg = f1_gse[inow - time_lin_reg: inow+1]
    f2_for_lin_reg = f2_gse[inow - time_lin_reg: inow+1]
    clf = linear_model.LinearRegression()
    clf.fit(np.array([f1_for_lin_reg,f2_for_lin_reg]).transpose(), kp_omni[inow - time_lin_reg: inow+1])
    a_f1 = clf.coef_[0]
    b_f2 = clf.coef_[1]
    c_f1f2 = clf.intercept_
    for i in range(nb_days_pred):
        kp_newell_lin_reg[i] = a_f1 * f1_gse[inow + i] + b_f2 * f2_gse[inow + i] + c_f1f2
        if kp_newell_lin_reg[i] >= 9:
            kp_newell_lin_reg[i] = 9

    rms_error_newell_lin_reg.append( np.sqrt(mean_squared_error(ap_omni[inow:inow+nb_days_pred], kp_newell_lin_reg))  )
    raise Exception

# RMS error. ASSUMPTION: date_omni and date_pred_swpc are the same arrays
## SWPC predictions with OMNI
rms_error_spwc = []
for ipred in range(n_swpc-4):
    rms_error_spwc.append ( np.sqrt(mean_squared_error(ap_omni[ipred+1:ipred+4], ap_pred_swpc[ipred, 1:])) )
        
## newell 2008 gsm with OMNI
rms_error_newell_gsm = []
for ipred in range(n_swpc-4):
    rms_error_newell_gsm.append( np.sqrt(mean_squared_error(ap_omni[ipred:ipred+4], ap_newell_2008_gsm[ipred:ipred+4])) )

## newell 2008 gse with OMNI
rms_error_newell_gse = []
for ipred in range(n_swpc-4):
    rms_error_newell_gse.append(  np.sqrt(mean_squared_error(ap_omni[ipred:ipred+4], ap_newell_2008_gse[ipred:ipred+4])) )

print np.mean(rms_error_spwc), np.mean(rms_error_newell_gsm), np.mean(rms_error_newell_gse), np.mean(rms_error_newell_lin_reg)

fig, ax = plt.subplots()
ap_swpc_one_day_ahead = np.zeros([n_swpc-1])
for i in range(1,n_swpc):
    ap_swpc_one_day_ahead[i-1] = ap_pred_swpc[i-1,1]
ax.plot(range(0, n_omni),ap_omni, 'k')
ax.plot(range(1, n_omni), ap_swpc_one_day_ahead, 'b')
ax.plot(ap_newell_2008_gse, 'r')

raise Exception
fig, ax = plt.subplots()
ax.plot(rms_error_spwc, 'k')
ax.plot(rms_error_newell_gsm, 'b')
ax.plot(rms_error_newell_gse, 'r')
ax.plot(rms_error_newell_lin_reg, 'g')
plt.show();plt.show();
