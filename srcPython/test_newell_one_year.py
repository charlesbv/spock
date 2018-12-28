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
filename_omni = 'kp_ap_files/all_info_one_year.txt'
file_omni = open(filename_omni, "r")
read_file_omni = file_omni.readlines()
n_header_omni = 0
n_elts_omni = len(read_file_omni) - n_header_omni
kp_omni = np.zeros([n_elts_omni]); ap_omni = np.zeros([n_elts_omni]); 
bx = np.zeros([n_elts_omni]); by = np.zeros([n_elts_omni]); bz = np.zeros([n_elts_omni]); 
t = np.zeros([n_elts_omni]); n = np.zeros([n_elts_omni]); v = np.zeros([n_elts_omni]); 
f107 = np.zeros([n_elts_omni]);
date_omni = []
for i in range(n_elts_omni):
    date_omni_temp = read_file_omni[i + n_header_omni].split()[0]+ '-' + read_file_omni[i + n_header_omni].split()[1]+ '-' + read_file_omni[i + n_header_omni].split()[2]
    date_omni.append( datetime.strptime(date_omni_temp, "%Y-%j-%H") )
    bx[i] = np.float(read_file_omni[i + n_header_omni].split()[3])
    by[i] = np.float(read_file_omni[i + n_header_omni].split()[4])
    bz[i] = np.float(read_file_omni[i + n_header_omni].split()[5])
    t[i] = np.float(read_file_omni[i + n_header_omni].split()[6])
    n[i] = np.float(read_file_omni[i + n_header_omni].split()[7])
    v[i] = np.float(read_file_omni[i + n_header_omni].split()[8])
    kp_omni[i] = np.float(read_file_omni[i + n_header_omni].split()[9]) / 10. # /10 because omniweb gives kp * 10
    ap_omni[i] = np.float(read_file_omni[i + n_header_omni].split()[10])
    f107[i] = np.float(read_file_omni[i + n_header_omni].split()[11])
    if ((bz[i] == 0) & (by[i] == 0)):
        bz[i] = 0.1 # to avoid 0/0, we set bz to this so that theta = 0
        

# 3 hourly average of these quantities





## Newell 2008 
bt = np.sqrt(by**2 + bz**2)
theta = np.arctan(by / bz)
f1 = v**(4./3) * bt**(2./3) * ((np.sin(theta/2.))**8)**(1./3)
f2 = np.sqrt(n)*v**2
kp_newell_2008 = 0.05 + 2.244 * 10**(-4) * f1 + 2.844*10**(-6) * f2


# Linear regression to determine the optimum coefficients a, b, and c for f1 and f2 so that: kp = a*f1 + b*f2 + c. In Newell 2008: a = 2.244 * 10**(-4), b = 2.844*10**(-6), and c = 0.05
time_lin_reg = 24
nb_days_pred = 3
nb_days_pred = 3 * 24 # in hours
step_analysis = 24
istep = -1
array_analysis = np.arange(time_lin_reg, n_elts_omni-nb_days_pred-time_lin_reg, step_analysis)
rms_newell_lin_reg_with_omni = np.zeros([len(array_analysis)])
rms_newell_2008_with_omni = np.zeros([len(array_analysis)])
correlation_newell_lin_reg_with_omni = np.zeros([len(array_analysis)])
correlation_newell_2008_with_omni = np.zeros([len(array_analysis)])
for inow in array_analysis:
    istep = istep + 1
    kp_newell_lin_reg = np.zeros([nb_days_pred+time_lin_reg])
    f1_for_lin_reg = f1[inow - time_lin_reg: inow+1]
    f2_for_lin_reg = f2[inow - time_lin_reg: inow+1]
    clf = linear_model.LinearRegression()
    clf.fit(np.array([f1_for_lin_reg,f2_for_lin_reg]).transpose(), kp_omni[inow - time_lin_reg: inow+1])
    a_f1 = clf.coef_[0]
    b_f2 = clf.coef_[1]
    c_f1f2 = clf.intercept_
    for i in range(nb_days_pred+time_lin_reg):
        kp_newell_lin_reg[i] = a_f1 * f1[inow + i - time_lin_reg] + b_f2 * f2[inow + i - time_lin_reg] + c_f1f2
        if kp_newell_lin_reg[i] >= 9:
            kp_newell_lin_reg[i] = 9

    rms_newell_lin_reg_with_omni[istep] = np.sqrt(mean_squared_error(kp_omni[inow: inow + nb_days_pred], kp_newell_lin_reg[time_lin_reg:time_lin_reg+nb_days_pred]))
    rms_newell_2008_with_omni[istep] = np.sqrt(mean_squared_error(kp_omni[inow: inow + nb_days_pred], kp_newell_2008[inow: inow + nb_days_pred]))

    # correlation_newell_lin_reg_with_omni[istep] = scipy.stats.pearsonr(kp_omni[inow: inow + nb_days_pred], kp_newell_lin_reg)[0]
    # correlation_newell_2008_with_omni[istep] = scipy.stats.pearsonr(kp_omni[inow: inow + nb_days_pred], kp_newell_2008[inow: inow + nb_days_pred])[0]


fig, ax = plt.subplots()
ax.plot(rms_newell_2008_with_omni, color = 'r')
ax.plot(rms_newell_lin_reg_with_omni, color = 'b')
raise Exception

fig, ax = plt.subplots()
ax.plot(kp_omni[inow-time_lin_reg: inow + nb_days_pred], color = 'k')
ax.plot(kp_newell_2008[inow-time_lin_reg: inow + nb_days_pred], color = 'r')
ax.plot(kp_newell_lin_reg, color = 'b')

