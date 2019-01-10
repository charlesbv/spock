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

# this scirpt crates the f107 and ap prediction error files that are usd in SpOCK to calculate the deviations of the prediction of f107 and Ap from the swpc predictions: swpc_mod deviation_f107_X_ap_Y where X and Y are the cvalues of the centers of the bin of the distribution of the difference between the historical prediciton by SWPC and the observations (see script new_uncertainty_f107_ap_bis.py to get these distirbutinos). The bins for F107 and Ap are both: [-15., -10.71428571,  -6.42857143,  -2.14285714,   2.14285714,   6.42857143,  10.71428571,  15.        ] so the centers of the bins, which are what we put in the deviation_f107_X_ap_Y files are: [-12.85714285,  -8.57142857,  -4.28571429,   0., 4.28571429,   8.57142857,  12.85714285]. So here we create 7*7 files, each file being a combinations of two values tgaken from these distrbutions (one for f107, one for ap). These values are constant over the 3 days of predictions. The files are put in the subfolder f107_ap


import numpy as np


f107_bin_range = np.array([-15.        , -10.71428571,  -6.42857143,  -2.14285714,         2.14285714,   6.42857143,  10.71428571,  15.        ])
ap_bin_range = f107_bin_range

f107_bin_center = (f107_bin_range[:-1] + np.roll(f107_bin_range, -1)[:-1])/2
ap_bin_center = (ap_bin_range[:-1] + np.roll(ap_bin_range, -1)[:-1])/2

nb_f107 = len(f107_bin_center)
nb_ap = len(ap_bin_center)
nb_comb = nb_f107 * nb_ap

for if107 in range(nb_f107):
    for iap in range(nb_ap):
        filename = 'f107_ap/deviation_f107_' + format(f107_bin_center[if107], ".0f") + '_ap_' + format(ap_bin_center[iap], ".0f") + '.txt'
        file = open(filename, "w")
        print >> file, '#1-sigma uncertainty on F10.7 and Ap at a function of the forecast time\n#Time(day) sigma_F10.7 sigma_Ap'
        for iday in range(1,46):
            print >> file, str(iday) + ' ' + format(f107_bin_center[if107], ".0f") + ' ' +  format(ap_bin_center[iap], ".0f")
        file.close()

