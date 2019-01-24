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

# This script plots the antenna patten. It reads in the antenna gain filesx
# ASSUMPTIONS:
# - see # PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
# - the files are ASCII files or binary. rows are theta, columns are phi
# - all antenna files mus thave the same format (over theta and phi)

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from struct import *



# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

filename_gain_list = ['/Users/cbv/cspice/data/ant_1_starboard_ddmi_v1_test.bin']
#filename_gain_list = ['/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/beacon/bruce/tds-bop V1.2.3/tds_antennaMap1_coarse.bin']


# [
#                       '/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/beacon/bruce/tds-bop V1.2.3/tds_antennaMap1_fine.bin',
#                       '/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/beacon/bruce/tds-bop V1.2.3/tds_antennaMap1_coarse.bin']
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

nb_file = len(filename_gain_list)
gain = []
phi_arr = []
theta_arr = []
for ifile in range(nb_file):
    #ifile = 0
    filename_gain = filename_gain_list[ifile]
    extension = filename_gain.split('.')[-1]
    if extension == 'bin':
        file_gain = open(filename_gain, "rb")
        x1 = unpack('i', file_gain.read(4))[0] # don t know what this is ....
        x2 = unpack('i', file_gain.read(4))[0] # don t know what this is ....
        numAz = unpack('i', file_gain.read(4))[0]
        numEl = unpack('i', file_gain.read(4))[0]
        az_start_deg = unpack('d', file_gain.read(8))[0]
        el_start_deg = unpack('d', file_gain.read(8))[0]
        az_inc_deg = unpack('d', file_gain.read(8))[0]
        el_inc_deg = unpack('d', file_gain.read(8))[0]
        #!!!! for some reason, a negative sign is needed for this file
        if 'tds_antennaMap1_coarse.bin' in filename_gain: 
            el_inc_deg = -el_inc_deg
        az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
        az_deg = np.arange(az_start_deg, az_stop_deg+az_inc_deg, az_inc_deg)
        el_stop_deg = el_start_deg + el_inc_deg * (numEl-1)
        el_deg = np.arange(el_start_deg, el_stop_deg+el_inc_deg, el_inc_deg)
        phi_arr_file = az_deg
        theta_arr_file = el_deg
        Z = np.zeros([numEl,numAz])
        for iaz in range(numAz):
            for iel in range(numEl):
                Z[iel, iaz] = unpack('d', file_gain.read(8))[0]
        file_gain.close
        gain_file = Z
        phi_arr.append(phi_arr_file)
        theta_arr.append(theta_arr_file)

    else: # if antenna gain file is not binary
        file_gain = open(filename_gain, "r")
        read_file_gain = file_gain.readlines()
        if ifile == 0:
            nb_theta = len(read_file_gain)
            nb_phi = len(read_file_gain[0].split(','))
            theta_max = 90. 
            dtheta = (int)( theta_max/nb_theta )
            theta_arr = np.arange(0, theta_max, dtheta)
            phi_max = 180. 
            dphi = (int)( phi_max * 2/nb_phi ) 
            phi_arr = np.arange(-180, phi_max, dphi)
            gain_file = np.zeros([nb_theta, nb_phi])
        for itheta in range(nb_theta):
            for iphi in range(nb_phi):
                gain_file[itheta, iphi] = np.float( read_file_gain[itheta].split(',')[iphi] )
    gain.append(gain_file)
max_gain = 15#np.max(gain)
# BEFORE JAN 22 2019
# for ifile in range(nb_file):
#     #ifile = 0
#     filename_gain = filename_gain_list[ifile]
#     file_gain = open(filename_gain, "r")
#     read_file_gain = file_gain.readlines()
#     if ifile == 0:
#         nb_theta = len(read_file_gain)
#         nb_phi = len(read_file_gain[0].split(','))
#         theta_max = 90. 
#         dtheta = (int)( theta_max/nb_theta )
#         theta_arr = np.arange(0, theta_max, dtheta)
#         phi_max = 180. 
#         dphi = (int)( phi_max * 2/nb_phi ) 
#         phi_arr = np.arange(-180, phi_max, dphi)
#         gain = np.zeros([nb_file, nb_theta, nb_phi])
#     for itheta in range(nb_theta):
#         for iphi in range(nb_phi):
#             gain[ifile, itheta, iphi] = np.float( read_file_gain[itheta].split(',')[iphi] )
# max_gain = 15#np.max(gain)
# end of  BEFORE JAN 22 2019
#raise Exception
# # Out[31]: (31.268536, 42.866177, [44, 0])
# theta = 31.268536
# phi = 42.866177
# itheta = np.where((theta_arr > theta))[0][0] - 1
# iphi = np.where((phi_arr > phi))[0][0] - 1

# PLOT
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
color_arr = ['b', 'r','cornflowerblue','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']


# Contour of delta Pc for uncertainty VS median f107 vs different altitudes
fig_title = ''
y_label = 'Theta '  + u'(\N{DEGREE SIGN})'
x_label = 'Phi '  + u'(\N{DEGREE SIGN})'


for ifile in range(nb_file):
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

    ax = fig.add_subplot(gs[0, 0])

    ax.set_title(filename_gain_list[ifile].split('/')[-1], weight = 'bold', fontsize  = fontsize_plot)
    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold


    origin = 'lower'

    x = phi_arr[ifile]
    y = theta_arr[ifile]
    X, Y = np.meshgrid(x, y)
    if el_start_deg < el_stop_deg: # file from low to high elevation
        # top of graph is first row of Z,
        # which should be highest elevation (since y axis is
        # low (bottom) to high (top) elevation) but since
        # el_start_deg < el_stop_deg then first row of gain[ifile]
        # is low elevation so invert it so that first row of Z is
        # highest elevation
        Z = np.zeros([numEl, numAz])
        for iel in range(numEl):
            Z[iel, :] =  gain[ifile][numEl-1-iel, :]
    else: # first row of gain[ifile] is highest elevation
        Z = gain[ifile]
        


    nr, nc = Z.shape

    Z = np.ma.array(Z)

    # We are using automatic selection of contour levels;
    # this is usually not such a good idea, because they don't
    # occur on nice boundaries, but we do it here for purposes
    # of illustration.
    levels = np.arange(0, max_gain + 1, 1.)
    # CS1 = ax.contourf(X, Y, Z, levels,
    #                   #[-1, -0.1, 0, 0.1],
    #                   #alpha=0.5,
    #                   cmap = plt.cm.get_cmap("jet"), # rainbow
    #                   origin=origin,
    #                  extend='both')

    # CS2 = ax.contour(X, Y, Z, levels,
    #                   colors=('k',),
    #                   linewidths=(1,),
    #                   origin=origin)

    # CS1.cmap.set_under('white')
    # CS1.cmap.set_over('white')
    CS1 = ax.imshow(Z, extent=[np.min(x), np.max(x), np.min(y), np.max(y)], cmap = 'jet', aspect = 'auto',
                    vmin = 0, vmax = max_gain, origin='upper')
    cbar = plt.colorbar(CS1, ax = ax)
    cbar.ax.set_ylabel('RCG', fontsize = fontsize_plot, weight = 'bold')
    filename_gain = filename_gain_list[ifile]
    extension = filename_gain.split('.')[-1]
    fig_save_name = filename_gain.replace(extension, 'pdf')
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



theta = 5.000000
phi = 240.000000
where_theta = np.where(theta_arr == theta)[0]
where_phi = np.where(phi_arr == phi)[0]

