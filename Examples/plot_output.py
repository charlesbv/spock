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

import sys
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from read_input_file import *
from read_output_file import *
from datetime import datetime
import numpy as np

# User inputs the filename
input_filename = sys.argv[1]

# Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3
color = ['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
x_label = 'Real time'

# Read the outputs
var_in, var_in_order = read_input_file(input_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 

dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];

var_to_read = ['longitude','latitude','altitude']
nvar = len(var_to_read)
for ivar in range(nvar):
    var = var_to_read[ivar]
    for isc in range(nb_sc):
        if ivar == 0:
            var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
            if isc == 0: # same date for all sc
                date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
                nb_seconds_since_start = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
                nb_steps_new = len(date) # in case the sc reentered the atmosphere before the end of the run
                date_ref = datetime.strptime(date[0], "%Y/%m/%d %H:%M:%S.%f")
                nb_seconds_in_simu = ( nb_steps_new - 1 ) * dt
                altitude = np.zeros([nb_sc, nb_steps_new]) # all output files of one simulation have the same number of steps
                latitude = np.zeros([nb_sc, nb_steps_new]) 
                longitude = np.zeros([nb_sc, nb_steps_new]) 
                x_axis = np.arange(0, nb_seconds_in_simu + 1, dt) 
            altitude[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'altitude')]
            latitude[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'latitude')]
            longitude[isc, :nb_steps_new] = var_out[find_in_read_input_order_variables(var_out_order, 'longitude')]

        # Plot 
        if ( isc == 0 ):
            if var == 'altitude':
                fig_title = 'Altitude as a function of time'
                y_label = 'Altitude (km)'
                fig_save_name = 'altitude.pdf'
            elif var == 'latitude':
                fig_title = 'Latitude as a function of time'
                y_label = 'Latitude (deg)'
                fig_save_name = 'latitude.pdf'
            elif var == 'longitude':
                fig_title = 'Longitude as a function of time'
                y_label = 'Longitude (deg)'
                fig_save_name = 'longitude.pdf'

            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

            fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
            plt.rc('font', weight='normal') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
            ax = fig.add_subplot(gs[0, 0])

            ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
            ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

            [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
            ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
            plt.rc('font', weight='normal') ## make the labels of the ticks in bold

        if var == 'altitude':
            y_axis = altitude[isc,:nb_steps_new]
        elif var == 'latitude':
            y_axis = latitude[isc,:nb_steps_new]
        elif var == 'longitude':
            y_axis = longitude[isc,:nb_steps_new]


        ax.plot(x_axis, y_axis, linewidth = 2, color = color[isc], label = 'SC #' + str(isc + 1))

        if isc == nb_sc - 1:
            legend = ax.legend(loc='upper left', bbox_to_anchor=(0, 1), numpoints = 1, fontsize = fontsize_plot)
            legend.get_title().set_fontsize(str(fontsize_plot))
            # x axis label is in real time
            nb_ticks_xlabel = 8
            dt_xlabel =  nb_seconds_in_simu / nb_ticks_xlabel # dt for ticks on x axis (in seconds)
            xticks = np.arange(0, nb_seconds_in_simu+1, dt_xlabel)
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x-xticks[0]) for x in xticks]
            for i in range(len(xticks)):
                if dt_xlabel >= 3*24*3600:
                    date_list_str.append( str(date_list[i])[5:10] )
                else:
                    date_list_str.append( str(date_list[i])[5:10] + "\n" + str(date_list[i])[11:16] )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)
            ax.margins(0,0); ax.set_xlim([min(xticks), max(xticks)])
            fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


