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
import matplotlib.gridspec as gridspec
import scipy
from sklearn.metrics import mean_squared_error
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from sklearn import linear_model


plt.ion()

# Omniweb
filename_omni = 'kp_ap_files/kp_ap_omni.txt'
file_omni = open(filename_omni, "r")
read_file_omni = file_omni.readlines()
n_header_omni = 0
while (read_file_omni[n_header_omni].split(' ')[0].replace("\n", "") != 'YEAR' ):
    n_header_omni = n_header_omni + 1
n_header_omni = n_header_omni + 1
n_elts_omni = len(read_file_omni) - n_header_omni
kp_omni = np.zeros([n_elts_omni]); ap_omni = np.zeros([n_elts_omni]); 
date_omni = []
for i in range(n_elts_omni):
    date_omni_temp = read_file_omni[i + n_header_omni].split()[0]+ '-' + read_file_omni[i + n_header_omni].split()[1]+ '-' + read_file_omni[i + n_header_omni].split()[2]
    date_omni.append( datetime.strptime(date_omni_temp, "%Y-%j-%H") )
    kp_omni[i] = np.float(read_file_omni[i + n_header_omni].split()[3]) / 10. # /10 because omniweb gives kp * 10
    ap_omni[i] = np.float(read_file_omni[i + n_header_omni].split()[3])



# Read GITM file to get IMF, v, and n
filename_gitm = '/raid3/Gitm/Runs/20150317/Imf_Files/IMF_full_180_14_input.dat'
file_gitm = open(filename_gitm, "r")
read_file_gitm= file_gitm.readlines()
n_header = 0
while (read_file_gitm[n_header].split(' ')[0].replace("\n", "") != '#START'):
    n_header= n_header + 1
n_header= n_header + 1
n_elts = len(read_file_gitm) - n_header
bx = np.zeros([n_elts]); by = np.zeros([n_elts]); bz = np.zeros([n_elts])
vx = np.zeros([n_elts]); vy = np.zeros([n_elts]); vz = np.zeros([n_elts])
n = np.zeros([n_elts])
t = np.zeros([n_elts])
date_gitm = []
for i in range(n_elts):
    date_gitm_temp = read_file_gitm[i+n_header].split()[0] + "-" + read_file_gitm[i+n_header].split()[1] + "-" + read_file_gitm[i+n_header].split()[2] + "T" + read_file_gitm[i+n_header].split()[3] + ":" + read_file_gitm[i+n_header].split()[4] 
    date_gitm.append( datetime.strptime(date_gitm_temp, "%Y-%m-%dT%H:%M") )
    bx[i] = np.float(read_file_gitm[i+n_header].split()[7])
    by[i] = np.float(read_file_gitm[i+n_header].split()[8])
    bz[i] = np.float(read_file_gitm[i+n_header].split()[9])
    vx[i] = np.float(read_file_gitm[i+n_header].split()[10])
    vy[i] = np.float(read_file_gitm[i+n_header].split()[11])
    vz[i] = np.float(read_file_gitm[i+n_header].split()[12])
    n[i] = np.float(read_file_gitm[i+n_header].split()[13])
    t[i] = np.float(read_file_gitm[i+n_header].split()[14])


# Define start and end time: start time is the newer time and end time is the older time between GITM and OMNI data
delta_time_start = ( date_gitm[0] - date_omni[0] ).days*24*3600 + ( date_gitm[0] - date_omni[0] ).seconds
if ( delta_time_start > 0 ):
    newer_time_start = date_gitm[0]
    older_time_start = date_omni[0]
else:
    newer_time_start = date_omni[0]
    older_time_start = date_gitm[0]

delta_time_stop = ( date_gitm[-1] - date_omni[-1] ).days*24*3600 + ( date_gitm[-1] - date_omni[-1] ).seconds
if ( delta_time_stop > 0 ): 
    older_time_stop = date_omni[-1]
else:
    older_time_stop = date_gitm[-1]

istart_omni = 0
while ( date_omni[istart_omni] < newer_time_start ):
    istart_omni = istart_omni + 1
istop_omni = 0
while ( date_omni[-1 - istop_omni] > older_time_stop ):
    istop_omni = istop_omni + 1
#print date_omni[istart_omni], date_omni[-1 - istop_omni] 

istart_gitm = 0
while ( date_gitm[istart_gitm] < newer_time_start ):
    istart_gitm = istart_gitm + 1
istop_gitm = 0
while ( date_gitm[-1 - istop_gitm] > older_time_stop ):
    istop_gitm = istop_gitm + 1
#print date_gitm[istart_gitm], date_gitm[-1 - istop_gitm] 

date_gitm_ok = date_gitm[istart_gitm:n_elts-istop_gitm]
bx_gitm_ok = bx[istart_gitm:n_elts-istop_gitm]
by_gitm_ok = by[istart_gitm:n_elts-istop_gitm]
bz_gitm_ok = bz[istart_gitm:n_elts-istop_gitm]
vx_gitm_ok = vx[istart_gitm:n_elts-istop_gitm]
vy_gitm_ok = vy[istart_gitm:n_elts-istop_gitm]
vz_gitm_ok = vz[istart_gitm:n_elts-istop_gitm]
n_gitm_ok = n[istart_gitm:n_elts-istop_gitm]
t_gitm_ok = t[istart_gitm:n_elts-istop_gitm]
n_elts_gitm_ok = len(n_gitm_ok)


# Hourly average imf, velocity, and n. The hourly average is from HH:01 to HH+1:00 included
i = 0
bx_gitm_ok_hourly_average = []
by_gitm_ok_hourly_average = []
bz_gitm_ok_hourly_average = []
vy_gitm_ok_hourly_average = []
vz_gitm_ok_hourly_average = []
n_gitm_ok_hourly_average = [] 
t_gitm_ok_hourly_average = [] 
vx_gitm_ok_hourly_average = []
date_gitm_ok_hourly_average = [] # date of the start of the bin in which the hourly average is derived
while ( i < n_elts_gitm_ok ):
    bx_gitm_ok_hourly_average_sum = bx_gitm_ok[i]
    by_gitm_ok_hourly_average_sum = by_gitm_ok[i]
    bz_gitm_ok_hourly_average_sum = bz_gitm_ok[i]
    vy_gitm_ok_hourly_average_sum = vy_gitm_ok[i]
    vz_gitm_ok_hourly_average_sum = vz_gitm_ok[i]
    n_gitm_ok_hourly_average_sum =  n_gitm_ok[i]
    t_gitm_ok_hourly_average_sum =  t_gitm_ok[i]
    vx_gitm_ok_hourly_average_sum =  vx_gitm_ok[i]
    n_sum = 1
    if i < n_elts_gitm_ok -1:
        date_gitm_ok_hourly_average.append(date_gitm_ok[i])
#        print date_gitm_ok[i]
        i = i + 1
        while ( (date_gitm_ok[i].minute !=0) | (date_gitm_ok[i].second !=0) ):
            bx_gitm_ok_hourly_average_sum = bx_gitm_ok_hourly_average_sum + bx_gitm_ok[i]
            by_gitm_ok_hourly_average_sum = by_gitm_ok_hourly_average_sum + by_gitm_ok[i]
            bz_gitm_ok_hourly_average_sum = bz_gitm_ok_hourly_average_sum + bz_gitm_ok[i]
            vy_gitm_ok_hourly_average_sum = vy_gitm_ok_hourly_average_sum + vy_gitm_ok[i]
            vz_gitm_ok_hourly_average_sum = vz_gitm_ok_hourly_average_sum + vz_gitm_ok[i]
            n_gitm_ok_hourly_average_sum =  n_gitm_ok_hourly_average_sum + n_gitm_ok[i]
            t_gitm_ok_hourly_average_sum = t_gitm_ok_hourly_average_sum +  t_gitm_ok[i]
            vx_gitm_ok_hourly_average_sum = vx_gitm_ok_hourly_average_sum + vx_gitm_ok[i]
            i = i + 1
            n_sum = n_sum + 1
            if i == n_elts_gitm_ok:
                break
        # print date_gitm_ok[i-1]
        # print n_sum
        # print ""
        bx_gitm_ok_hourly_average.append( bx_gitm_ok_hourly_average_sum / n_sum )
        by_gitm_ok_hourly_average.append( by_gitm_ok_hourly_average_sum / n_sum )
        bz_gitm_ok_hourly_average.append( bz_gitm_ok_hourly_average_sum / n_sum )
        vy_gitm_ok_hourly_average.append( vy_gitm_ok_hourly_average_sum / n_sum )
        vz_gitm_ok_hourly_average.append( vz_gitm_ok_hourly_average_sum / n_sum )
        n_gitm_ok_hourly_average.append( n_gitm_ok_hourly_average_sum / n_sum ) 
        t_gitm_ok_hourly_average.append( t_gitm_ok_hourly_average_sum / n_sum ) 
        vx_gitm_ok_hourly_average.append( vx_gitm_ok_hourly_average_sum / n_sum )

    else:
        break



# Do another round of time start and end: we don't want to take values from gitm if their corresponding time did not start at HH:00
if ( (date_gitm_ok_hourly_average[0].minute != 0) | (date_gitm_ok_hourly_average[0].second != 0) ):
    date_gitm_ok_hourly_average.remove(date_gitm_ok_hourly_average[0])
    bx_gitm_ok_hourly_average.remove(bx_gitm_ok_hourly_average[0])
    by_gitm_ok_hourly_average.remove(by_gitm_ok_hourly_average[0])
    bz_gitm_ok_hourly_average.remove(bz_gitm_ok_hourly_average[0])
    vy_gitm_ok_hourly_average.remove(vy_gitm_ok_hourly_average[0])
    vz_gitm_ok_hourly_average.remove(vz_gitm_ok_hourly_average[0])
    n_gitm_ok_hourly_average.remove(n_gitm_ok_hourly_average[0]) 
    t_gitm_ok_hourly_average.remove(t_gitm_ok_hourly_average[0]) 
    vx_gitm_ok_hourly_average.remove(vx_gitm_ok_hourly_average[0])
if ( (date_gitm_ok_hourly_average[-1].minute != 0) | (date_gitm_ok_hourly_average[-1].second != 0) ):
    date_gitm_ok_hourly_average.remove(date_gitm_ok_hourly_average[-1])
    bx_gitm_ok_hourly_average.remove(bx_gitm_ok_hourly_average[-1])
    by_gitm_ok_hourly_average.remove(by_gitm_ok_hourly_average[-1])
    bz_gitm_ok_hourly_average.remove(bz_gitm_ok_hourly_average[-1])
    vy_gitm_ok_hourly_average.remove(vy_gitm_ok_hourly_average[-1])
    vz_gitm_ok_hourly_average.remove(vz_gitm_ok_hourly_average[-1])
    n_gitm_ok_hourly_average.remove(n_gitm_ok_hourly_average[-1]) 
    t_gitm_ok_hourly_average.remove(t_gitm_ok_hourly_average[-1]) 
    vx_gitm_ok_hourly_average.remove(vx_gitm_ok_hourly_average[-1])



## Define start and end time: start time is the newer time and end time is the older time between GITM and OMNI data
delta_time_start = ( date_gitm_ok_hourly_average[0] - date_omni[0] ).days*24*3600 + ( date_gitm_ok_hourly_average[0] - date_omni[0] ).seconds
if ( delta_time_start > 0 ):
    newer_time_start = date_gitm_ok_hourly_average[0]
    older_time_start = date_omni[0]
else:
    newer_time_start = date_omni[0]
    older_time_start = date_gitm_ok_hourly_average[0]

delta_time_stop = ( date_gitm_ok_hourly_average[-1] - date_omni[-1] ).days*24*3600 + ( date_gitm_ok_hourly_average[-1] - date_omni[-1] ).seconds
if ( delta_time_stop > 0 ): 
    older_time_stop = date_omni[-1]
else:
    older_time_stop = date_gitm_ok_hourly_average[-1]

istart_omni = 0
while ( date_omni[istart_omni] < newer_time_start ):
    istart_omni = istart_omni + 1
istop_omni = 0
while ( date_omni[-1 - istop_omni] > older_time_stop ):
    istop_omni = istop_omni + 1
#print date_omni[istart_omni], date_omni[-1 - istop_omni] 

istart_gitm = 0
while ( date_gitm_ok[istart_gitm] < newer_time_start ):
    istart_gitm = istart_gitm + 1
istop_gitm = 0
while ( date_gitm_ok[-1 - istop_gitm] > older_time_stop ):
    istop_gitm = istop_gitm + 1

date_omni_ok = date_omni[istart_omni:n_elts_omni - istop_omni]
kp_omni_ok = kp_omni[istart_omni:n_elts_omni - istop_omni]
ap_omni_ok = kp_omni[istart_omni:n_elts_omni - istop_omni]
n_elts_omni_ok = len(kp_omni_ok)

date_gitm_ok = date_gitm[istart_gitm:n_elts-istop_gitm]
bx_gitm_ok = bx[istart_gitm:n_elts-istop_gitm]
by_gitm_ok = by[istart_gitm:n_elts-istop_gitm]
bz_gitm_ok = bz[istart_gitm:n_elts-istop_gitm]
vx_gitm_ok = vx[istart_gitm:n_elts-istop_gitm]
vy_gitm_ok = vy[istart_gitm:n_elts-istop_gitm]
vz_gitm_ok = vz[istart_gitm:n_elts-istop_gitm]
n_gitm_ok = n[istart_gitm:n_elts-istop_gitm]
t_gitm_ok = t[istart_gitm:n_elts-istop_gitm]
n_elts_gitm_ok = len(n_gitm_ok)

v_gitm_ok = np.sqrt(vx_gitm_ok**2 + vy_gitm_ok**2 + vz_gitm_ok**2)
bt_gitm_ok = np.sqrt(by_gitm_ok**2 + bz_gitm_ok**2 )

## Newell 2008 using not hourly average parameters
bt_gitm_ok = np.sqrt(by_gitm_ok**2 + bz_gitm_ok**2)
theta_gitm_ok = np.arctan(by_gitm_ok / bz_gitm_ok)
f1 = v_gitm_ok**(4./3) * bt_gitm_ok**(2./3) * ((np.sin(theta_gitm_ok/2.))**8)**(1./3)
f2 = np.sqrt(n_gitm_ok)*v_gitm_ok**2
kp_newell_ok = 0.05 + 2.244 * 10**(-4) * f1 + 2.844*10**(-6) * f2




bx_gitm_ok_hourly_average = np.array(bx_gitm_ok_hourly_average)
by_gitm_ok_hourly_average = np.array(by_gitm_ok_hourly_average)
bz_gitm_ok_hourly_average = np.array(bz_gitm_ok_hourly_average)
vy_gitm_ok_hourly_average = np.array(vy_gitm_ok_hourly_average)
vz_gitm_ok_hourly_average = np.array(vz_gitm_ok_hourly_average)
n_gitm_ok_hourly_average = np.array(n_gitm_ok_hourly_average) 
t_gitm_ok_hourly_average = np.array(t_gitm_ok_hourly_average) 
vx_gitm_ok_hourly_average = np.array(vx_gitm_ok_hourly_average)


n_elts_gitm_ok_hourly_average = len(date_gitm_ok_hourly_average)
v_gitm_ok_hourly_average = np.sqrt(vx_gitm_ok_hourly_average**2 + vy_gitm_ok_hourly_average**2 + vz_gitm_ok_hourly_average**2 )
b_gitm_ok_hourly_average = np.sqrt(bx_gitm_ok_hourly_average**2 + by_gitm_ok_hourly_average**2 + bz_gitm_ok_hourly_average**2 )

# Newell 2008 using hourly_average parameters
bt_gitm_ok_hourly_average = np.sqrt(by_gitm_ok_hourly_average**2 + bz_gitm_ok_hourly_average**2)
theta_gitm_ok_hourly_average = np.arctan(by_gitm_ok_hourly_average / bz_gitm_ok_hourly_average)
f1 = v_gitm_ok_hourly_average**(4./3) * bt_gitm_ok_hourly_average**(2./3) * ((np.sin(theta_gitm_ok_hourly_average/2.))**8)**(1./3)
f2 = np.sqrt(n_gitm_ok_hourly_average)*v_gitm_ok_hourly_average**2
kp_newell_ok_hourly_average = 0.05 + 2.244 * 10**(-4) * f1 + 2.844*10**(-6) * f2



# Linear regression to determine the optimum coefficients a, b, and c for f1 and f2 so that: kp = a*f1 + b*f2 + c. In Newell 2008: a = 2.244 * 10**(-4), b = 2.844*10**(-6), and c = 0.05
####### IMPORTANT: date_omni_ok and date_gitm_ok_hourly_average need to be equal. Indeed, f1 and f2 have to be computed at the same time (date_gitm_ok_hourly_average) as kp is computed (date_omni_ok)
now = 75
time_lin_reg = 58
f1_for_lin_reg = f1[now - time_lin_reg: now+1]
f2_for_lin_reg = f2[now - time_lin_reg: now+1]
clf = linear_model.LinearRegression()
clf.fit(np.array([f1_for_lin_reg,f2_for_lin_reg]).transpose(), kp_omni_ok[now - time_lin_reg: now+1])
a_f1 = clf.coef_[0]
b_f2 = clf.coef_[1]
c_f1f2 = clf.intercept_
kp_newell_lin_reg = np.zeros([n_elts_gitm_ok_hourly_average])
for i in range(n_elts_gitm_ok_hourly_average):
    kp_newell_lin_reg[i] = a_f1 * f1[i] + b_f2 * f2[i] + c_f1f2
    if kp_newell_lin_reg[i] >= 9:
        kp_newell_lin_reg[i] = 9
clf2 = linear_model.LinearRegression()
clf2.fit(np.array([f1,f2]).transpose(), kp_omni_ok)
kp_newell_lin_reg2 = clf2.predict(np.array([f1,f2]).transpose())

rms_newell_lin_reg_with_omni = np.sqrt(mean_squared_error(kp_omni_ok, kp_newell_lin_reg))
rms_newell_2008_with_omni = np.sqrt(mean_squared_error(kp_omni_ok, kp_newell_ok_hourly_average))

correlation_newell_lin_reg_with_omni = scipy.stats.pearsonr(kp_omni_ok, kp_newell_lin_reg)[0]
correlation_newell_2008_with_omni = scipy.stats.pearsonr(kp_omni_ok, kp_newell_ok_hourly_average)[0]

# print rms_newell_lin_reg_with_omni, correlation_newell_lin_reg_with_omni
# print rms_newell_2008_with_omni, correlation_newell_2008_with_omni

# Plots

width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax = fig.add_subplot(gs[0, 0])
ax.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure

# the three lines below work because date_omni_ok[0] = date_gitm_ok[0] = date_gitm_ok_hourly_average[0] (which is always true)
omni_seconds_since_start = np.array( [ (t-date_omni_ok[0]).total_seconds() for t in date_omni_ok ])
gitm_seconds_since_start = np.array( [ (t-date_gitm_ok[0]).total_seconds() for t in date_gitm_ok ])
gitm_seconds_since_start_hourly_average = np.array( [ (t-date_gitm_ok_hourly_average[0]).total_seconds() for t in date_gitm_ok_hourly_average ])
ax.plot(omni_seconds_since_start/3600., kp_omni_ok, linewidth = 2, color = 'r')
#ax.plot(gitm_seconds_since_start, kp_newell_ok, linewidth = 2, color = 'b')
ax.plot(gitm_seconds_since_start_hourly_average/3600., kp_newell_ok_hourly_average, linewidth = 2, color = 'g', label = '2008 - ' + '{0:.2f}'.format(rms_newell_2008_with_omni) + ' - ' + '{0:.3f}'.format(correlation_newell_2008_with_omni))
ax.plot(gitm_seconds_since_start_hourly_average/3600., kp_newell_lin_reg2, linewidth = 2, color = 'b',linestyle = '', marker = '.')#, label = 'Lin reg - ' + '{0:.2f}'.format(rms_newell_lin_reg_with_omni) + ' - ' + '{0:.3f}'.format(correlation_newell_lin_reg_with_omni))
ax.plot(gitm_seconds_since_start_hourly_average/3600., kp_newell_lin_reg, linewidth = 2, color = 'b', label = 'Lin reg 24h - ' + '{0:.2f}'.format(rms_newell_lin_reg_with_omni) + ' - ' + '{0:.3f}'.format(correlation_newell_lin_reg_with_omni))
ax.plot([now,now],[0,ax.get_ylim()[1]], linestyle = 'dotted', linewidth = 2, color = 'k')
ax.legend(loc = 2, fontsize = fontsize_plot)
ax.set_xlabel('Time (hours)', fontsize = fontsize_plot, weight = 'bold')
ax.set_ylabel('Kp', fontsize = fontsize_plot, weight = 'bold')
ax.text(now, ax.get_ylim()[0] + (ax.get_ylim()[1]-ax.get_ylim()[0])/40., 'Now', fontsize = fontsize_plot)
ax.axvspan( now - 6, now, color = 'b', alpha = 0.2)
ax.text(now-3, ax.get_ylim()[1] - (ax.get_ylim()[1]-ax.get_ylim()[0])/6., 'Lin\nreg',color = 'b', horizontalalignment = 'center', fontsize = fontsize_plot)
plt.show();plt.show();

# fig, ax = plt.subplots()
# dref = date_gitm_ok_hourly_average[0]
# ax.plot([ (t-dref).total_seconds() for t in date_gitm_ok ], v_gitm_ok, color = 'b')
# ax.plot([ (t-dref).total_seconds() for t in date_gitm_ok_hourly_average ], v_gitm_ok_hourly_average, color = 'r', linewidth = 4)


# fig, ax = plt.subplots()
# dref = date_gitm_ok_hourly_average[0]
# ax.plot([ (t-dref).total_seconds() for t in date_gitm_ok ], bt_gitm_ok, color = 'b')
# ax.plot([ (t-dref).total_seconds() for t in date_gitm_ok_hourly_average ], bt_gitm_ok_hourly_average, color = 'r', linewidth = 4)








# # Borovosky 2013b equation 79
# ## read F10.7 from omniweb
# filename_omni = 'kp_ap_files/f107_omni.txt'
# file_omni = open(filename_omni, "r")
# read_file_omni = file_omni.readlines()
# n_header_omni = 0
# while (read_file_omni[n_header_omni].split(' ')[0].replace("\n", "") != 'YEAR' ):
#     n_header_omni = n_header_omni + 1
# n_header_omni = n_header_omni + 1
# n_elts_omni = len(read_file_omni) - n_header_omni
# f107_omni = np.zeros([n_elts_omni]); ap_omni = np.zeros([n_elts_omni]); 
# for i in range(n_elts_omni):
#     f107_omni[i] = np.float(read_file_omni[i + n_header_omni].split()[3])
# file_omni.close()

# ## Calculate kp. Note: typical values are brom borovosky 2008
# mu_0 = 4*np.pi * 10**(-7)
# mp = 1.6726 * 10**(-27) # proton masss in kg
# rho = n*100*100*100*mp # total mass density of the charged plasma particles. 100*100*100 because n is expressed in cm^-3 and we want rho in kg * m^-3. Typical value of rho*v**2: 1.33 to 5.34 nPa
# va = b*10**(-9) / np.sqrt(mu_0 * rho) # 10**9 because b is given in nT. Typical value: 50 km/s
# ma = v*10**3 / va # Typical value: ma = 500 / 50 = 10
# f1_bor = -np.sqrt(ma/3.18)
# sigma_p = 0.77 * np.sqrt(np.mean(f107_omni)) # !!!!! should be np.sqrt(f107_omni) but since f107_omni does not have the same dimensions as va and that I don't want to do a linear interpolation, we take the mean. Typical value: sigma_p = 5 mho
# q = va*10**(-3) * sigma_p / 796. # va must be in km/s here
# abohm = 1.72 * 10**(-4) # table 12
# c = ( 2.44*10**(-4) + ( 1 + 1.38*np.log(ma) )**(-6) )**(-1./6)
# beta_s = (ma/6.)**1.92
# gamma = 5./3 # !!!!!!!! no idea if this is the correct value for the solar wind
# kb = 1.38064852 * 10**(-23)
# cs = np.sqrt( gamma * kb * 2*t / mp ) # !!!!!!!  It assumees t_e = t_i AND that the temperature in the file is t_e. Typical value: cs= 20 km/s
# mms  = v*10**3 / np.sqrt( va**2 + cs**2 ) # Typical value: mms = 7.5 to 15
# w = beta_s * c**(-1./2) * np.sqrt( 1+0.5*mms**(-2) ) * (1+beta_s)**(-1./2)
# b_bor =  np.sqrt((n)) * (v)**(5./2) * c**(-1./2) * w**(1./2)
# g_app = (np.sin(theta/2.))**2 * (n)**(1./3) * (v)**(5./3) * ma**(-0.7084) * (np.exp(f1_bor))/(1+q)
# kp_borovsky =  g_app + abohm * b_bor

# kp_borovsky_eq_66a =  np.zeros([n_elts])
# for i in range(n_elts):
#     kp_borovsky_eq_66a[i] =  14.6 * (g_app[i] + b_bor[i])**2.01 / ( 16.6 + 0.877*(g_app[i]+b_bor[i])**2.01 )

# #kp_borovsky14_eq_8 = 0.4 * np.sqrt(mu_0 * mp) * (np.sin(theta/2.))**2 * c**(-1./2) * np.sqrt(n) * v**2 * (1+beta_s)**(-3./4)

