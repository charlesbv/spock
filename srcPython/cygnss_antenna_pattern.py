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

# This script plots the antenna patten. It reads in the antenna gain files

# INPUTS:
# - filename_gain_list: name (including path) of the antenna gain pattern file
# - res_map: resolution of the map (set to 'coarse' or 'fine')
# ASSUMPTIONS:
# - the files are ASCII files (e.g., .txt, .agm) or binary
# - if the file is binary, it must be the same format as for tds-bop.exe:
#   - format of header:
#     - 2 first elements (1-2) are uint32 but we don't know what they are
#     - elements 3-4 are uint32: numAz and numEl
#     - elements 5-6 are double: az_start_deg and el_start_deg
#     - elements 7-8 are double: az_inc_deg and el_inc_deg
#   - the subsequent elemnts (9 to end) are the gain values
#   - if the resolution res_map is set to "coarse" then, for an unknow
#   reason, we need to force a minus sign on the el_inc_deg
# - if the file is not binary:
#   - if the resolution res_map is set to "fine" then
#   it must be in the same format as the files created by Scott Gleason (and
#   subequently by Darren McKague). FileFormatNotes.txt describes the format in
#   details:
#     - elev 90 to 0, step 0.1 deg (901 values) (decreasing order)
#     - azim 0 to 360, step 0.1 deg (3601 values) (ascending order, 0 to 180 is
#       starboard, 180 to 360 is port)
#     - first 2 lines of files are header, then line 3 to end are gain values
#   - if res_map is set to "coarse", then it must be in the same format as the
#    on-board algorithm maps:
#     - elev 0 to 85 deg, step 5 deg (18 values) (ascending order)
#     - azim -180 to 165, step 15 deg (24 values) (ascending order, -180 to 0 is
#       port, 0 to 180 is starboard)
#     - no header line, directly start with the gain values
#   - be careful:
#     - in the coarse case, the first azim elements are for the port side (-180
#       to 0), while in the fine case, the first azim elements are for the
#       starboard side (0 to 180)
#     - in the coarse case, rows are elevations and columns are azimuths, while
#       in the fine case, rows are azimuths and columns are elevations
#     - the fine maps have one extra value in elev and azim





# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
res_map = 'fine' # coarse or fine. To set only if the file is ont binary

filename_gain_list = ['/Users/cbv/work/spockOut/beacon/ant_data/ant_1_starboard_v6.txt',
'/Users/cbv/work/spockOut/beacon/ant_data/ant_1_port_v6.txt']


# Coarse .agm
# ['/Users/cbv/cspice/data/ant_1_port_ddmi_v1.agm',
# '/Users/cbv/cspice/data/ant_1_starboard_ddmi_v1.agm'] 

# Fine .txt
# ['/Users/cbv/work/spockOut/beacon/ant_data/ant_1_starboard_v6.txt',
# '/Users/cbv/work/spockOut/beacon/ant_data/ant_1_port_v6.txt']

# Coarse .bin
# ['/Users/cbv/cspice/data/ant_1_port_ddmi_v1_test.bin']   
# ['/Users/cbv/cspice/data/merged_ant_1_starboard_ddmi_v1_with_ant_1_port_ddmi_v1_test.bin']
#['/Users/cbv/cspice/data/ant_1_port_ddmi_v1_test.bin']

# Coarse .bin TDS
# ['/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/beacon/bruce/tds-bop V1.2.3/tds_antennaMap1_coarse.bin']

# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from struct import *
import sys
import ipdb


if ((res_map != 'coarse') & (res_map != 'fine')):
    sys.exit("!***!\n res_map must be set to 'coarse' or 'fine'\
            \n!***!")
                    
print 'Map resolution: ' + res_map

nb_file = len(filename_gain_list)
gain = []
az_deg = []
el_deg = []
el_start_deg = []
el_stop_deg = []
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
        el_start_deg_file = unpack('d', file_gain.read(8))[0]
        az_inc_deg = unpack('d', file_gain.read(8))[0]
        el_inc_deg = unpack('d', file_gain.read(8))[0]
        #!!!! for some reason, a negative sign is needed for this file
        if res_map == 'coarse':
            el_inc_deg = -el_inc_deg
        az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
        az_deg_file = np.linspace(az_start_deg, az_stop_deg, numAz)
        el_stop_deg_file = el_start_deg_file + el_inc_deg * (numEl-1)
        el_deg_file = np.linspace(el_start_deg_file, el_stop_deg_file, numEl)
        Z = np.zeros([numEl,numAz])
        for iaz in range(numAz):
            for iel in range(numEl):
                Z[iel, iaz] = unpack('d', file_gain.read(8))[0]
        file_gain.close
        gain_file = Z

    else: # if antenna gain file is not binary
        if res_map == 'fine':
            nheader = 2
            file_gain = open(filename_gain, "r")
            read_file_gain = file_gain.readlines()
            numAz = len(read_file_gain) - nheader
            numEl = len(read_file_gain[0+nheader].split(','))
            # Elevation
            el_start_deg_file = 90. 
            el_inc_deg = -0.1 # FileFormatNotes.txt: step is 0.1 deg
            el_stop_deg_file = el_start_deg_file + el_inc_deg * (numEl-1)
            # in FileFormatNotes.txt, elev goes in decreasing order
            el_deg_file = np.linspace(el_start_deg_file, el_stop_deg_file, numEl)
            # Azimuth
            az_start_deg = 0. 
            az_inc_deg =  0.1 # FileFormatNotes.txt: step is 0.1 deg
            az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
            az_deg_file = np.linspace(az_start_deg, az_stop_deg, numAz)
            # Gain
            gain_file = np.zeros([numEl, numAz])
            for iphi in range(numAz):
                print iphi, numAz-1
                for itheta in range(numEl):
                    gain_file[itheta, iphi] = np.float( read_file_gain[iphi+
                    nheader].split(',')[itheta] )
            
        elif res_map == 'coarse':
            nheader = 0
            file_gain = open(filename_gain, "r")
            read_file_gain = file_gain.readlines()
            numEl = len(read_file_gain) - nheader
            numAz = len(read_file_gain[0+nheader].split(','))
            # Elevation
            el_start_deg_file = 0. 
            el_inc_deg = 5. # FileFormatNotes.txt: step is 0.1 deg
            el_stop_deg_file = el_start_deg_file + el_inc_deg * (numEl-1)
            # in FileFormatNotes.txt, elev goes in decreasing order
            el_deg_file = np.linspace(el_start_deg_file, el_stop_deg_file, numEl)
            # Azimuth
            az_start_deg = -180. 
            az_inc_deg =  15. # FileFormatNotes.txt: step is 0.1 deg
            az_stop_deg = az_start_deg + az_inc_deg * (numAz-1)
            az_deg_file = np.linspace(az_start_deg, az_stop_deg, numAz)
            # Gain
            gain_file = np.zeros([numEl, numAz])

            for itheta in range(numEl):
                print itheta, numEl-1
                for iphi in range(numAz):
                    gain_file[itheta, iphi] = np.float( read_file_gain[itheta+
                    nheader].split(',')[iphi] )

    gain.append(gain_file)
    az_deg.append(az_deg_file)
    el_deg.append(el_deg_file)
    el_start_deg.append(el_start_deg_file)
    el_stop_deg.append(el_stop_deg_file)

    
max_gain = 15#np.max(gain)

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

    x = az_deg[ifile]
    y = el_deg[ifile]
    X, Y = np.meshgrid(x, y)
    extension = filename_gain.split('.')[-1]
    #    if extension == 'bin':
    if el_start_deg[ifile] < el_stop_deg[ifile]: # file from low to high elevation
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
    # else:
    #     if res_map  == 'coarse': 
    #         numEl = gain[ifile].shape[0]
    #         numAz = gain[ifile].shape[1]
    #         Z = np.zeros([numEl, numAz])
    #         for iel in range(numEl):
    #             Z[iel, :] =  gain[ifile][numEl-1-iel, :]

    #     else:
    #         Z = gain[ifile]
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
    fig_save_name = filename_gain.replace(extension, 'pdf')
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

