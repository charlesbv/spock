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

from convert_tle_date_to_date import *
from norad_id_to_cygnss_name import *
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
import pickle
import sys
import fileinput
import time
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
from get_name_mission import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta

######### PARAMETERS TO SET #########
tle_filename = "./cygnss/tle_12-18_to_01-30.txt"

## Parameters for the figure
root_save_fig_name = './cygnss/' # folder where figures are saved
save_results = 0 # set to 1 to save the results
show_plots = 1 # set to 1 to show_plot

height_fig = 9.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
nb_ticks_xlabel = 6.

######### ALGORITHM #########
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius        = 6378.137; # mean equatorial radius (km)

if show_plots == 1:
    plt.ion()


tle_file = open(tle_filename, "r")
read_tle_file = tle_file.readlines()

ecc = []
sma = []
dsmadday = []
alt_perigee = []
alt_apogee = []
date = []
sc_name = []

n = len(read_tle_file) # number of lines in the TLE file
iline = 0
while iline < n: # go over each TLE in the file
    ecc_per_sc = []
    sma_per_sc = []
    dsmadday_per_sc = [] 
    alt_apogee_per_sc = []
    alt_perigee_per_sc = []
    date_per_sc = []
    sc_name.append( norad_id_to_cygnss_name( read_tle_file[iline+1].split()[1] ) )
    sc_name_temp = sc_name[-1]
    while ( ( sc_name_temp == sc_name[-1] ) & ( iline < n ) ): # the file includes blocks of TLE for each sc. Here go over a block for a given sc
        ecc_per_sc.append( np.float( '0.' +  read_tle_file[iline+1].split()[4] ) )
        mean_motion = np.float( read_tle_file[iline+1].split()[7] ) 
        period = 24 * 3600. / mean_motion # mean_motion is in rev/day
        sma_per_sc.append( ( ( period / ( 2*np.pi ) )**2 * earth_mu ) ** (1./3) )
        alt_perigee_per_sc.append( sma_per_sc[-1] * ( 1 - ecc_per_sc[-1] ) - earth_radius )
        alt_apogee_per_sc.append( sma_per_sc[-1] * ( 1 + ecc_per_sc[-1] ) - earth_radius )
        date_per_sc.append( read_tle_file[iline].split()[3] )
        if len(sma_per_sc) > 1:
            current_date = convert_tle_date_to_date( date_per_sc[-1] )
            previous_date = convert_tle_date_to_date( date_per_sc[-2] )
            if current_date != previous_date:
                delta_seconds = ( current_date - previous_date ).total_seconds()
                dsmadday_per_sc.append( ( sma_per_sc[-1] - sma_per_sc[-2] ) / delta_seconds * 3600 * 24 / sma_per_sc[0])
            else: # ometimes 2 TLEs have the same epoch. in that case the slope should be the same as the previous slope
                dsmadday_per_sc.append( dsmadday_per_sc[-1] )
        if len(sma_per_sc) > 2:
            dsmadday_per_sc[-1] = ( dsmadday_per_sc[-1] + dsmadday_per_sc[-2] ) / 2
                
        iline = iline + 2
        if iline < n:
            sc_name_temp = norad_id_to_cygnss_name( read_tle_file[iline+1].split()[1] )
    ecc.append( ecc_per_sc )
    sma.append( sma_per_sc )
    dsmadday.append( dsmadday_per_sc )
    alt_apogee.append( alt_apogee_per_sc )
    alt_perigee.append( alt_perigee_per_sc )
    date.append( date_per_sc )
#raise Exception
nb_sc = len(ecc)

# For plots, generate disctinct colors
NCURVES = nb_sc
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)


date_ref = '2016-12-18T00:00'
date_ref = datetime.strptime( date_ref, "%Y-%m-%dT%H:%M" )

save_last_tle_date = np.zeros([nb_sc])

# Eccentricity plot
if "eccentricity" in sys.argv:
    for isc in range(nb_sc):
        nb_tle_for_sc = len(ecc[isc])
        nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
        for itle in range(nb_tle_for_sc):
            nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
        save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
        if isc == 0:
            fig_title = 'Eccentricity from TLEs for each CYGNSS VS time'
            y_label = 'Eccentricity'
            x_label = 'Real time'

        x_axis = nb_sec_since_date_ref 
        y_axis = ecc[isc]

        if isc == 0:
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle(fig_title, y = 0.973,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left = 0.17, right=0.82, top = 0.93,bottom = 0.12, hspace = 0.01)
            ax = fig.add_subplot(gs[0, 0])

            ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
            ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

            [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
            ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        colorVal = scalarMap.to_rgba(isc)
        ax.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = str(sc_name[isc]))

        if isc == nb_sc - 1:
            # Label of x axis is real time
            most_recent_tle_among_all_sc = np.max(save_last_tle_date) # this is the max of seconds ellapsed among all CYGNSS between the latest tle and the date_ref
            dt_xlabel = most_recent_tle_among_all_sc / nb_ticks_xlabel
            xticks = np.arange(0, most_recent_tle_among_all_sc, dt_xlabel)
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x) for x in xticks]
            for i in range(len(xticks)):
                date_list_str.append( str(date_list[i])[5:10] )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
    #        ax.margins(0,0)
            ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])


            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="CYGNSS", fontsize = fontsize_plot)
            legend.get_title().set_fontsize(str(fontsize_plot))


            # Save/show results
            if save_results == 1:
                fig_save_name = '0109_eccentricity'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            if show_plots == 1:
                plt.show(); plt.show()

if "sma" in sys.argv:
    # Sma plot
    isc = 6
    if isc == 6:
    # for isc in range(nb_sc):
        nb_tle_for_sc = len(ecc[isc])
        nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
        for itle in range(nb_tle_for_sc):
            nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
        save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
        if isc == 6:
            #        fig_title = 'SMA - ' + r'$\mathbf{R_E}$' + ' from TLEs for each CYGNSS VS time'
            fig_title = 'SMA(t) / SMA(t0) from TLEs for each CYGNSS VS time'
            y_label = 'SMA(t) / SMA(t0)'
            x_label = 'Real time'
        x_axis = nb_sec_since_date_ref 
        y_axis = ( np.array(sma[isc]) ) / ( np.array(sma[isc][0]) )
        y_axis_sma = y_axis
        x_axis_sma = x_axis
        
        if isc == 6:
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle(fig_title, y = 0.973,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left = 0.17, right=0.82, top = 0.93,bottom = 0.12, hspace = 0.01)
            ax = fig.add_subplot(gs[0, 0])

            ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
            ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

            [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
            ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        colorVal = scalarMap.to_rgba(isc)
        ax.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = str(sc_name[isc]))

        if isc == 6:#nb_sc - 1:
            # Label of x axis is real time
            most_recent_tle_among_all_sc = np.max(save_last_tle_date) # this is the max of seconds ellapsed among all CYGNSS between the latest tle and the date_ref
            dt_xlabel = most_recent_tle_among_all_sc / nb_ticks_xlabel
            xticks = np.arange(0, most_recent_tle_among_all_sc+1, dt_xlabel)
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x) for x in xticks]
            for i in range(len(xticks)):
                date_list_str.append( str(date_list[i])[5:10] )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
            ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])
            ax.margins(0,0)
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="CYGNSS", fontsize = fontsize_plot)
            legend.get_title().set_fontsize(str(fontsize_plot))

            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.6f'))

            # Save/show results
            if save_results == 1:
                fig_save_name = '0109_sma'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            if show_plots == 1:
                plt.show(); plt.show()
            

if "dsma" in sys.argv:
    # Dsmadday plot (sma(day) - sma(day-1))
#    for isc in range(nb_sc):
    isc = 6
    if isc == 6:
        if ( (sc_name[isc] != "FM03")):
            nb_tle_for_sc = len(ecc[isc])
            nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
            for itle in range(nb_tle_for_sc):
                nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
            save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
            if isc == 6:
                fig_title = 'SMA(day) - SMA(day-1) from TLEs for each CYGNSS VS time'
                y_label = 'DSMADDAY - ' + r'$\mathbf{R_E}$' + ' (km)'
                x_label = 'Real time'

            x_axis = nb_sec_since_date_ref 
            y_axis = np.concatenate((np.array([0]), np.array(dsmadday[isc])), axis = 0)

            if isc == 6:
                fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
                fig.suptitle(fig_title, y = 0.973,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold
                gs = gridspec.GridSpec(1, 1)
                gs.update(left = 0.17, right=0.82, top = 0.93,bottom = 0.12, hspace = 0.01)
                ax = fig.add_subplot(gs[0, 0])

                ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
                ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

                [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
                ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                plt.rc('font', weight='bold') ## make the labels of the ticks in bold

            colorVal = scalarMap.to_rgba(isc)
            ax.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = str(sc_name[isc]))

            if isc == 6:#nb_sc - 1:
                # Label of x axis is real time
                most_recent_tle_among_all_sc = np.max(save_last_tle_date) # this is the max of seconds ellapsed among all CYGNSS between the latest tle and the date_ref
                dt_xlabel = most_recent_tle_among_all_sc / nb_ticks_xlabel
                xticks = np.arange(0, most_recent_tle_among_all_sc, dt_xlabel)
                date_list_str = []
                date_list = [date_ref + timedelta(seconds=x) for x in xticks]
                for i in range(len(xticks)):
                    date_list_str.append( str(date_list[i])[5:10] )
                ax.xaxis.set_ticks(xticks)
                ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])
        #        ax.margins(0,0)
                legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="CYGNSS", fontsize = fontsize_plot)
                legend.get_title().set_fontsize(str(fontsize_plot))


                # Save/show results
                if save_results == 1:
                    fig_save_name = '0109_dsmadday'
                    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                if show_plots == 1:
                    plt.show(); plt.show()


if "alt_p" in sys.argv:
    # Altitude perigee plot
    for isc in range(nb_sc):
        nb_tle_for_sc = len(ecc[isc])
        nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
        for itle in range(nb_tle_for_sc):
            nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
        save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
        if isc == 0:
            fig_title = 'Altitude of perigee from TLEs for each CYGNSS VS time'
            y_label = 'Altitude of perigee (km)'
            x_label = 'Real time'

        x_axis = nb_sec_since_date_ref 
        y_axis = alt_perigee[isc]

        if isc == 0:
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle(fig_title, y = 0.973,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left = 0.17, right=0.82, top = 0.93,bottom = 0.12, hspace = 0.01)
            ax = fig.add_subplot(gs[0, 0])

            ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
            ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

            [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
            ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        colorVal = scalarMap.to_rgba(isc)
        ax.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = str(sc_name[isc]))

        if isc == nb_sc - 1:
            # Label of x axis is real time
            most_recent_tle_among_all_sc = np.max(save_last_tle_date) # this is the max of seconds ellapsed among all CYGNSS between the latest tle and the date_ref
            dt_xlabel = most_recent_tle_among_all_sc / nb_ticks_xlabel
            xticks = np.arange(0, most_recent_tle_among_all_sc, dt_xlabel)
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x) for x in xticks]
            for i in range(len(xticks)):
                date_list_str.append( str(date_list[i])[5:10] )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
            ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])
    #        ax.margins(0,0)
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="CYGNSS", fontsize = fontsize_plot)
            legend.get_title().set_fontsize(str(fontsize_plot))


            # Save/show results
            if save_results == 1:
                fig_save_name = '0109_alt_perigee'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            if show_plots == 1:
                plt.show(); plt.show()
            
if "alt_a" in sys.argv:
    # Altitude apogee plot
    for isc in range(nb_sc):
        nb_tle_for_sc = len(ecc[isc])
        nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
        for itle in range(nb_tle_for_sc):
            nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
        save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
        if isc == 0:
            fig_title = 'Altitude of apogee from TLEs for each CYGNSS VS time'
            y_label = 'Altitude of apogee (km)'
            x_label = 'Real time'

        x_axis = nb_sec_since_date_ref 
        y_axis = alt_apogee[isc]

        if isc == 0:
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle(fig_title, y = 0.973,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold
            gs = gridspec.GridSpec(1, 1)
            gs.update(left = 0.17, right=0.82, top = 0.93,bottom = 0.12, hspace = 0.01)
            ax = fig.add_subplot(gs[0, 0])

            ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
            ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

            [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
            ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        colorVal = scalarMap.to_rgba(isc)
        ax.plot(x_axis, y_axis, linewidth = 2, color = colorVal, label = str(sc_name[isc]))

        if isc == nb_sc - 1:
            # Label of x axis is real time
            most_recent_tle_among_all_sc = np.max(save_last_tle_date) # this is the max of seconds ellapsed among all CYGNSS between the latest tle and the date_ref
            dt_xlabel = most_recent_tle_among_all_sc / nb_ticks_xlabel
            xticks = np.arange(0, most_recent_tle_among_all_sc, dt_xlabel)
            date_list_str = []
            date_list = [date_ref + timedelta(seconds=x) for x in xticks]
            for i in range(len(xticks)):
                date_list_str.append( str(date_list[i])[5:10] )
            ax.xaxis.set_ticks(xticks)
            ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
            ax.set_xlim([ax.get_xlim()[0], most_recent_tle_among_all_sc])
    #        ax.margins(0,0)
            legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="CYGNSS", fontsize = fontsize_plot)
            legend.get_title().set_fontsize(str(fontsize_plot))


            # Save/show results
            if save_results == 1:
                fig_save_name = '0109_alt_apogee'
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
            if show_plots == 1:
                plt.show(); plt.show()
            

# smatemp = np.array(sma[0])
# n = len(smatemp)
# nb_interpo_point = 10
# ninterpo = (n-1)*nb_interpo_point + 1
# smainterpo = np.zeros([ninterpo])
# smainterpo[-1] = smatemp[-1]
# for i in range(ninterpo-1):
#     j = (int)(i / nb_interpo_point)
#     smainterpo[i] = ( smatemp[j+1] - smatemp[j] ) / nb_interpo_point * i +  smatemp[j] 




# FOR CHRIS, perigee and apogee vs time
file_sc = open("./cygnss/radius_apogee_perigee.txt", "w")                
print >> file_sc, "This file shows the radii of perigee and apogee for each CYGNSS. They have been calculated from the TLEs at spacetrack.org, using the eccentricity and the mean motion (#rev/day, that was then converted to semi-major axis). The altitudes of perigee and apogee are taken to be radius - earth mean radius (6378.137 km). You can contact cbv@umich.edu for questions."
for isc in range(nb_sc):
    print >> file_sc, "CYG"+sc_name[isc]
    print >> file_sc,  "time altitude_perigee(km) radius_perigee(km) altitude_apogee(km) radius_apogee(km)"
    nb_tle_for_sc = len(ecc[isc])
    for itle in range(nb_tle_for_sc):
        date_tle = convert_tle_date_to_date( date[isc][itle] )
        print >> file_sc, date_tle, alt_perigee[isc][itle], alt_perigee[isc][itle] + earth_radius, alt_apogee[isc][itle],alt_apogee[isc][itle]+earth_radius
    if isc < nb_sc - 1:
        print  >> file_sc,""
file_sc.close()
