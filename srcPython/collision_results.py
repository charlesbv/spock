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
# This script read the file created by list_read_collision_file_new.py (usually on Big). This file shows the Pc and TCA for a list of runs. This script plots some statistics about this.

import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from datetime import datetime, timedelta



filename_results = 'results_combined_with_results_complete.txt' 
#for inclination angle difference: 'results_combined_with_results_complete.txt'. 
#For eccentricity: results_encounterecc_argper180.txt. 
# for bc: results_encounterbc.txt
#For f107/ap: results_encounterf107.txt
#(results_encounterecc.txt is for collision apogee but at local time 180 deg off compared to runs from other stduies so different temperature so differnt density)
file_results = open(filename_results)
read_file_results = file_results.readlines()
nb_run = len(read_file_results)
pc = np.zeros([nb_run])
tca = []
filename_run = []
alt = np.zeros([nb_run, 2])
inc = np.zeros([nb_run, 2])
arg_per = np.zeros([nb_run, 2])
raan = np.zeros([nb_run, 2])
true_ano = np.zeros([nb_run, 2])
ecc = np.zeros([nb_run, 2])
bc = np.zeros([nb_run])
f107 = np.zeros([nb_run])
quantile = np.zeros([nb_run])
alt_arr = []; inc_arr = []; arg_per_arr = []; raan_arr = []; true_ano_arr = []; ecc_arr = []; bc_arr = [];f107_arr = []; quantile_arr = []
alt0_arr = []; inc0_arr = []; arg_per0_arr = []; raan0_arr = []; true_ano0_arr = []; ecc0_arr = []
alt1_arr = []; inc1_arr = []; arg_per1_arr = []; raan1_arr = []; true_ano1_arr = []; ecc1_arr = []
for irun in range(nb_run):
    filename_run.append( read_file_results[irun].split()[0] )
    pc[irun] =  np.float( read_file_results[irun].split()[1] )
    tca.append( read_file_results[irun].split()[2]  )#datetime.strptime(read_file_results[irun].split()[2] , "%Y-%m-%dT%H:%M:%S.%f") )
    alt[irun,0] = np.float( filename_run[-1].split('alt')[1].split('_')[0].split('-')[0] )
    alt[irun,1] = np.float( filename_run[-1].split('alt')[1].split('_')[0].split('-')[1] )
    inc[irun,0] = np.float( filename_run[-1].split('inc')[1].split('_')[0].split('-')[0] )
    inc[irun,1] = np.float( filename_run[-1].split('inc')[1].split('_')[0].split('-')[1] )
    arg_per[irun,0] = np.float( filename_run[-1].split('arg_per')[1].split('_')[0].split('-')[0] )
    arg_per[irun,1] = np.float( filename_run[-1].split('arg_per')[1].split('_')[0].split('-')[1] )
    raan[irun,0] = np.float( filename_run[-1].split('raan')[1].split('_')[0].split('-')[0] )
    raan[irun,1] = np.float( filename_run[-1].split('raan')[1].split('_')[0].split('-')[1] )
    true_ano[irun,0] = np.float( filename_run[-1].split('true_ano')[1].split('_')[0].split('-')[0] )
    true_ano[irun,1] = np.float( filename_run[-1].split('true_ano')[1].split('_')[0].split('-')[1] )
    ecc[irun,0] = np.float( filename_run[-1].split('ecc')[1].split('_')[0].split('-')[0] )
    ecc[irun,1] = np.float( '-'.join(filename_run[-1].split('ecc')[1].split('_')[0].split('-')[1:])  )
    if 'bc' in filename_run[-1]:
        bc[irun] = np.float( filename_run[-1].split('bc')[1].split('_')[0] )
    f107[irun] = np.float( filename_run[-1].split('f107')[1].split('_')[0] )
    quantile[irun] = np.float( filename_run[-1].split('quantile')[1].split('.')[0].split('_')[0] )


    if ( alt[irun,0] in alt_arr ) == False:
        alt_arr.append(alt[irun,0]) # all altitude, regardeless of sc 0 or 1
    if ( inc[irun,0] in inc_arr ) == False:
        inc_arr.append(inc[irun,0])
    if ( arg_per[irun,0] in arg_per_arr ) == False:
        arg_per_arr.append(arg_per[irun,0])
    if ( raan[irun,0] in raan_arr ) == False:
        raan_arr.append(raan[irun,0])
    if ( true_ano[irun,0] in true_ano_arr ) == False:
        true_ano_arr.append(true_ano[irun,0])
    if ( ecc[irun,0] in ecc_arr ) == False:
        ecc_arr.append(ecc[irun,0])

    if ( alt[irun,1] in alt_arr ) == False:
        alt_arr.append(alt[irun,1]) # all altitude, regardeless of sc 0 or 1
    if ( inc[irun,1] in inc_arr ) == False:
        inc_arr.append(inc[irun,1])
    if ( arg_per[irun,1] in arg_per_arr ) == False:
        arg_per_arr.append(arg_per[irun,1])
    if ( raan[irun,1] in raan_arr ) == False:
        raan_arr.append(raan[irun,1])
    if ( true_ano[irun,1] in true_ano_arr ) == False:
        true_ano_arr.append(true_ano[irun,1])
    if ( ecc[irun,1] in ecc_arr ) == False:
        ecc_arr.append(ecc[irun,1])
    if ( bc[irun] in bc_arr ) == False:
        bc_arr.append(bc[irun])

    if ( f107[irun] in f107_arr ) == False:
        f107_arr.append(f107[irun])
    if ( quantile[irun] in quantile_arr ) == False:
        quantile_arr.append(quantile[irun])


    if ( alt[irun,0] in alt0_arr ) == False:
        alt0_arr.append(alt[irun,0]) #altitudes of only sc 0
    if ( alt[irun,1] in alt1_arr ) == False:
        alt1_arr.append(alt[irun,1]) #altitudes of only sc 1
    if ( inc[irun,0] in inc0_arr ) == False:
        inc0_arr.append(inc[irun,0])
    if ( inc[irun,1] in inc1_arr ) == False:
        inc1_arr.append(inc[irun,1])
    if ( arg_per[irun,0] in arg_per0_arr ) == False:
        arg_per0_arr.append(arg_per[irun,0])
    if ( arg_per[irun,1] in arg_per1_arr ) == False:
        arg_per1_arr.append(arg_per[irun,1])
    if ( raan[irun,0] in raan0_arr ) == False:
        raan0_arr.append(raan[irun,0]) 
    if ( raan[irun,1] in raan1_arr ) == False:
        raan1_arr.append(raan[irun,1]) 
    if ( true_ano[irun,0] in true_ano0_arr ) == False:
        true_ano0_arr.append(true_ano[irun,0])
    if ( true_ano[irun,1] in true_ano1_arr ) == False:
        true_ano1_arr.append(true_ano[irun,1])
    if ( ecc[irun,0] in ecc0_arr ) == False:
        ecc0_arr.append(ecc[irun,0]) 
    if ( ecc[irun,1] in ecc1_arr ) == False:
        ecc1_arr.append(ecc[irun,1]) 

alt_arr = np.sort(np.array(alt_arr))
inc_arr = np.sort(np.array(inc_arr))
arg_per_arr = np.sort(np.array(arg_per_arr))
raan_arr = np.sort(np.array(raan_arr))
true_ano_arr = np.sort(np.array(true_ano_arr))
ecc_arr = np.sort(np.array(ecc_arr))
bc_arr = np.sort(np.array(bc_arr))
f107_arr = np.sort(np.array(f107_arr))
quantile_arr = np.sort(np.array(quantile_arr))

alt0_arr = np.sort(np.array(alt0_arr))
inc0_arr = np.sort(np.array(inc0_arr))
arg_per0_arr = np.sort(np.array(arg_per0_arr))
raan0_arr = np.sort(np.array(raan0_arr))
true_ano0_arr = np.sort(np.array(true_ano0_arr))
ecc0_arr = np.sort(np.array(ecc0_arr))

alt1_arr = np.sort(np.array(alt1_arr))
inc1_arr = np.sort(np.array(inc1_arr))
arg_per1_arr = np.sort(np.array(arg_per1_arr))
raan1_arr = np.sort(np.array(raan1_arr))
true_ano1_arr = np.sort(np.array(true_ano1_arr))
ecc1_arr = np.sort(np.array(ecc1_arr))


nb_alt = len(alt_arr)
nb_inc = len(inc_arr)
nb_arg_per = len(arg_per_arr)
nb_raan = len(raan_arr)
nb_true_ano = len(true_ano_arr)
nb_ecc = len(ecc_arr)
nb_bc = len(bc_arr)
nb_f107 = len(f107_arr)
nb_quantile = len(quantile_arr)

nb_alt0 = len(alt0_arr)
nb_inc0 = len(inc0_arr)
nb_arg_per0 = len(arg_per0_arr)
nb_raan0 = len(raan0_arr)
nb_true_ano0 = len(true_ano0_arr)
nb_ecc0 = len(ecc0_arr)

nb_alt1 = len(alt1_arr)
nb_inc1 = len(inc1_arr)
nb_arg_per1 = len(arg_per1_arr)
nb_raan1 = len(raan1_arr)
nb_true_ano1 = len(true_ano1_arr)
nb_ecc1 = len(ecc1_arr)


nb_run_diff_oe_f107 = nb_alt * nb_inc * nb_arg_per * nb_raan * nb_true_ano * nb_ecc * nb_bc * nb_f107
# gather all runs that have the same oe and 107 together. So they only differ by their quantile
pc_arr = np.zeros([nb_alt0, nb_alt1, nb_inc0, nb_inc1, nb_arg_per0, nb_arg_per1, nb_raan0, nb_raan1, nb_true_ano0, nb_true_ano1, nb_ecc0, nb_ecc1, nb_bc, nb_f107, nb_quantile]) - 1 # if no run for a given oe/bc/f107/quantile then pc is -1
filename_run_arr = np.chararray([nb_alt0, nb_alt1, nb_inc0, nb_inc1, nb_arg_per0, nb_arg_per1, nb_raan0, nb_raan1, nb_true_ano0, nb_true_ano1, nb_ecc0, nb_ecc1, nb_bc, nb_f107, nb_quantile], itemsize = 500)
tca_arr = np.chararray([nb_alt0, nb_alt1, nb_inc0, nb_inc1, nb_arg_per0, nb_arg_per1, nb_raan0, nb_raan1, nb_true_ano0, nb_true_ano1, nb_ecc0, nb_ecc1, nb_bc, nb_f107, nb_quantile], itemsize = 30)
pc_median = np.zeros([nb_alt0, nb_alt1, nb_inc0, nb_inc1, nb_arg_per0, nb_arg_per1, nb_raan0, nb_raan1, nb_true_ano0, nb_true_ano1, nb_ecc0, nb_ecc1, nb_bc, nb_f107]) - 1 # if no run for a given oe/bc/f107/quantile then pc is -1
for irun in range(nb_run):
    irun, nb_run - 1
    ialt0 = np.where( alt0_arr == alt[irun, 0] )[0][0]
    ialt1 = np.where( alt1_arr == alt[irun, 1] )[0][0]
    iinc0 = np.where( inc0_arr == inc[irun, 0] )[0][0]
    iinc1 = np.where( inc1_arr == inc[irun, 1] )[0][0]
    iarg_per0 = np.where( arg_per0_arr == arg_per[irun, 0] )[0][0]
    iarg_per1 = np.where( arg_per1_arr == arg_per[irun, 1] )[0][0]
    iraan0 = np.where( raan0_arr == raan[irun, 0] )[0][0]
    iraan1 = np.where( raan1_arr == raan[irun, 1] )[0][0]
    itrue_ano0 = np.where( true_ano0_arr == true_ano[irun, 0] )[0][0]
    itrue_ano1 = np.where( true_ano1_arr == true_ano[irun, 1] )[0][0]
    iecc0 = np.where( ecc0_arr == ecc[irun, 0] )[0][0]
    iecc1 = np.where( ecc1_arr == ecc[irun, 1] )[0][0]
    ibc = np.where( bc_arr == bc[irun] )[0][0]
    if107 = np.where( f107_arr == f107[irun] )[0][0]
    iquantile = np.where( quantile_arr == quantile[irun] )[0][0]
    
    pc_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, iquantile] = pc[irun]
    filename_run_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, iquantile] = filename_run[irun]
    tca_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, iquantile] = tca[irun]

    if quantile[irun] == 5:
        pc_median[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] = pc[irun]


# PLOT
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3
color_arr = ['b', 'r','cornflowerblue','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']


raise Exception



# Contour of delta Pc for uncertainty VS median f107 vs different altitudes
fig_title = 'Contour map of the relative uncertainty in $\mathrm{P_{c}}$\nas a function of F10.7 and altitude'
y_label = 'F10.7'
x_label = 'Altitude (km)'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.99,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold


origin = 'lower'
#origin = 'upper'


iinc0 = 0
iarg_per0 = 0
iraan0 = 0
itrue_ano0 = 0
iecc0 = 0
iinc1 = 0
iarg_per1 = 0
iraan1 = 0
itrue_ano1 = 0
iecc1 = 0
ibc = 0
# delta = 0.025
# x = y = np.arange(-3.0, 3.01, delta)
# X, Y = np.meshgrid(x, y)
x = alt0_arr
y = f107_arr
X, Y = np.meshgrid(x, y)
Z = np.zeros([nb_f107, nb_alt0])



for if107 in  range(nb_f107):
    for ialt in range(nb_alt0):

        pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
        Z[if107, ialt] = ( pc_max_here - pc_min_here ) / pc_nom * 100 


nr, nc = Z.shape

Z = np.ma.array(Z)

# We are using automatic selection of contour levels;
# this is usually not such a good idea, because they don't
# occur on nice boundaries, but we do it here for purposes
# of illustration.
levels = np.arange(0,100+5,5)
CS1 = ax.contourf(X, Y, Z, levels,
                  #[-1, -0.1, 0, 0.1],
                  #alpha=0.5,
                  cmap = plt.cm.get_cmap("rainbow"),
                  origin=origin,
                 extend='both')

CS2 = ax.contour(X, Y, Z, levels,
                  colors=('k',),
                  linewidths=(3,),
                  origin=origin)

CS1.cmap.set_under('white')
CS1.cmap.set_over('white')

cbar = plt.colorbar(CS1, ax = ax)#, title = '$(\mathrm{P_{c, max} - P_{c, min}})/\mathrm{P_{c, nom}}$')
cbar.ax.set_ylabel('$(\mathrm{P_{c, max} - P_{c, min}})/\mathrm{P_{c, nom}}$', fontsize = fontsize_plot, weight = 'bold')

fig_save_name = 'contour_delta_pc_vs_f107_more_level.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Uncertainty VS median f107 for different altitudes
fig_title = 'Relative uncertainty in $\mathrm{P_{c}}$ VS F10.7 (nominal) at different altitudes'
y_label = '$(\mathrm{P_{c, max} - P_{c, min}})/\mathrm{P_{c, nom}}$'
x_label = 'Nominal F10.7'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

iinc0 = 0
iarg_per0 = 0
iraan0 = 0
itrue_ano0 = 0
iecc0 = 0
iinc1 = 0
iarg_per1 = 0
iraan1 = 0
itrue_ano1 = 0
iecc1 = 0
ibc = 0
max_yaxis = -1
alt_arr_temp = [300, 360, 400, 460, 500]
for ialt_temp in range(len(alt_arr_temp)):#nb_alt): # here sc 0 and sc 1 have the same altitude
    ialt = np.where(alt_arr == alt_arr_temp[ialt_temp])[0][0]
    x_axis = np.zeros([nb_f107])
    y_axis = np.zeros([nb_f107])
    print ialt, nb_alt
    for if107 in range(nb_f107):
        x_axis[if107] = f107_arr[if107]
        pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
        y_axis[if107] = ( pc_max_here - pc_min_here ) / pc_nom * 100
        print f107_arr[if107]
        print pc_min_here, pc_max_here, pc_nom, ( pc_max_here - pc_min_here ) / pc_nom * 100
        print ''
    if np.max(y_axis) > max_yaxis:
        max_yaxis = np.max(y_axis)
    ax.plot(x_axis, y_axis, linewidth = 2, color = color_arr[ialt_temp], label = str(alt_arr[ialt]))
    ax.scatter(x_axis, y_axis, linewidth = 2, color = color_arr[ialt_temp])

ax.set_ylim([0, 105])  # ax.set_ylim([0, max_yaxis+10])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Altitude", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
fig_save_name = 'delta_pc_vs_f107.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception


# Uncertainty VS bc for different altitudes
fig_title = 'Relative uncertainty in $\mathrm{P_{c}}$ VS ballistic coefficient at different altitudes'
y_label = '$(\mathrm{P_{c, max} - P_{c, min}})/\mathrm{P_{c, nom}}$'
x_label = '$B_c$'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

iinc0 = 0
iarg_per0 = 0
iraan0 = 0
itrue_ano0 = 0
iecc0 = 0
iinc1 = 0
iarg_per1 = 0
iraan1 = 0
itrue_ano1 = 0
iecc1 = 0
if107 = 0
max_yaxis = -1
for ialt in range(nb_alt): # here sc 0 and sc 1 have the same altitude
    x_axis = np.zeros([nb_bc])
    y_axis = np.zeros([nb_bc])
    print ialt, nb_alt
    for ibc in range(nb_bc):
        x_axis[ibc] = bc_arr[ibc]
        pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
        y_axis[ibc] = ( pc_max_here - pc_min_here ) / pc_nom * 100
        print bc_arr[ibc]
        print pc_min_here, pc_max_here, pc_nom, ( pc_max_here - pc_min_here ) / pc_nom * 100
        print ''
    if np.max(y_axis) > max_yaxis:
        max_yaxis = np.max(y_axis)
    ax.plot(x_axis, y_axis, linewidth = 2, color = color_arr[ialt], label = str(alt_arr[ialt]))
    ax.scatter(x_axis, y_axis, linewidth = 2, color = color_arr[ialt])

ax.set_ylim([0, 105])  # ax.set_ylim([0, max_yaxis+10])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Altitude", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
fig_save_name = 'delta_pc_vs_bc.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

raise Exception


# Uncertainty VS eccentricity difference for different altitudes
fig_title = 'Relative uncertainty in $\mathrm{P_{c}}$ VS eccentricity of sc 2 at different altitudes'
y_label = '$(\mathrm{P_{c, max} - P_{c, min}})/\mathrm{P_{c, nom}}$'
x_label = 'Eccentricity of spacecraft 2'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

iarg_per0 = 0
iraan0 = 0
itrue_ano0 = 0
iinc0 = 0
iarg_per1 = 0
iraan1 = 0
itrue_ano1 = 0
iinc1 = 0
ibc = 0
if107= 0
max_yaxis = -1
for ialt in range(nb_alt): # here sc 0 and sc 1 have the same altitude
    x_axis = np.zeros([nb_ecc1])
    y_axis = np.zeros([nb_ecc1])
    print alt1_arr[ialt]
    for iecc1 in range(nb_ecc1):
        x_axis[iecc1] = ecc1_arr[iecc1] 
        pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
        pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
        y_axis[iecc1] = ( pc_max_here - pc_min_here ) / pc_nom * 100
        print ecc1_arr[iecc1]
        #            print pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]
        print pc_min_here, pc_max_here, pc_nom, ( pc_max_here - pc_min_here ) / pc_nom * 100
        print ''
    if np.max(y_axis) > max_yaxis:
        max_yaxis = np.max(y_axis)
    ax.semilogx(x_axis, y_axis, linewidth = 2, color = color_arr[ialt], label = str(alt_arr[ialt]))
    ax.scatter(x_axis, y_axis, linewidth = 2, color = color_arr[ialt])

ax.set_ylim([0, 105])  # ax.set_ylim([0, max_yaxis+10])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Altitude", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
fig_save_name = 'delta_pc_vs_eccentricity_argper180.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# alt0 = 400; alt1 = 400
# inc0 = 90; inc1 = 30
# arg_per0 = 0; arg_per1 = 0
# raan0 = 0; raan1 = 0
# true_ano0 = 0; true_ano1 = 0  
# ecc0 = 0; ecc1 = 0
# f10701 = 100;
# quantile01 = 3

# ialt0 = np.where( alt0_arr == alt0 )[0][0]
# ialt1 = np.where( alt1_arr == alt1 )[0][0]
# iinc0 = np.where( inc0_arr == inc0 )[0][0]
# iinc1 = np.where( inc1_arr == inc1 )[0][0]
# iarg_per0 = np.where( arg_per0_arr == arg_per0 )[0][0]
# iarg_per1 = np.where( arg_per1_arr == arg_per1 )[0][0]
# iraan0 = np.where( raan0_arr == raan0 )[0][0]
# iraan1 = np.where( raan1_arr == raan1 )[0][0]
# itrue_ano0 = np.where( true_ano0_arr == true_ano0 )[0][0]
# itrue_ano1 = np.where( true_ano1_arr == true_ano1 )[0][0]
# iecc0 = np.where( ecc0_arr == ecc0 )[0][0]
# iecc1 = np.where( ecc1_arr == ecc1 )[0][0]
# if107 = np.where( f107_arr == f10701 )[0][0]
# iquantile = np.where( quantile_arr == quantile01 )[0][0]

# print alt0, alt1, inc0, inc1, arg_per0, arg_per1, raan0, raan1, true_ano0, true_ano1, ecc0, ecc1, f10701, quantile01
# print pc_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile]
# print filename_run_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile]
# print tca_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile]
# print ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile








# Uncertainty VS inclination difference for different altitudes
fig_title = 'Relative uncertainty in $\mathrm{P_{c}}$ VS absolute difference in inclination at different altitudes'
y_label = '$(\mathrm{P_{c, max} - P_{c, min}})/\mathrm{P_{c, nom}}$'
x_label = 'Absolute difference in inclination ' + u'(\N{DEGREE SIGN})'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

nb_inc_combinations = nb_inc0 * nb_inc1
iarg_per0 = 0
iraan0 = 0
itrue_ano0 = 0
iecc0 = 0
iarg_per1 = 0
iraan1 = 0
itrue_ano1 = 0
iecc1 = 0
ibc = 0
if107= 0
max_yaxis = -1
for ialt in range(nb_alt): # here sc 0 and sc 1 have the same altitude
    x_axis = np.zeros([nb_inc_combinations])
    y_axis = np.zeros([nb_inc_combinations])
    print ialt, nb_alt
    iinc_comb = 0
    for iinc0 in range(nb_inc0):
        for iinc1 in range(nb_inc1):
            x_axis[iinc_comb] = np.abs(inc1_arr[iinc1] - inc0_arr[iinc0])
            pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
            pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
            pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
            y_axis[iinc_comb] = ( pc_max_here - pc_min_here ) / pc_nom * 100
            iinc_comb = iinc_comb + 1
            print inc1_arr[iinc1], x_axis[iinc_comb-1]
#            print pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]
            print pc_min_here, pc_max_here, pc_nom, ( pc_max_here - pc_min_here ) / pc_nom * 100
            print ''
    if np.max(y_axis) > max_yaxis:
        max_yaxis = np.max(y_axis)
    ax.plot(x_axis, y_axis, linewidth = 2, color = color_arr[ialt], label = str(alt_arr[ialt]))
    ax.scatter(x_axis, y_axis, linewidth = 2, color = color_arr[ialt])

ax.set_ylim([0, 105])  # ax.set_ylim([0, max_yaxis+10])

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Altitude", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
fig_save_name = 'delta_pc_vs_difference_inclination_all.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



# Pc VS quantiles at two different altitude for a given inclination of satellite 2
fig_title = ''
y_label = '$P_c$'
x_label = 'Quantile of the density distribution'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in bold

nb_inc_combinations = nb_inc0 * nb_inc1
iarg_per0 = 0
iraan0 = 0
itrue_ano0 = 0
iecc0 = 0
iarg_per1 = 0
iraan1 = 0
itrue_ano1 = 0
iecc1 = 0
ibc = 0
if107= 0
max_yaxis = -1
iinc0 = 0
x_axis = np.arange(1,10)

ialt = 2
iinc1 = 0
y_axis = pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]
pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
delta_pc = ( pc_max_here - pc_min_here ) / pc_nom * 100
ax.plot(x_axis, y_axis, linewidth = 2, color = 'b', label =  str(alt1_arr[ialt]) + ' km | ' + str(np.abs(inc1_arr[iinc1] - inc0_arr[iinc0])) +  u'\N{DEGREE SIGN} | ' + format(delta_pc, ".0f") + '%')
ax.scatter(x_axis, y_axis, linewidth = 2, color = 'b')
iinc1 = -1
pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
delta_pc = ( pc_max_here - pc_min_here ) / pc_nom * 100
y_axis = pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]
ax.plot(x_axis, y_axis, linewidth = 2, color = 'b', label =  str(alt1_arr[ialt]) + ' km | ' + str(np.abs(inc1_arr[iinc1] - inc0_arr[iinc0]))+  u'\N{DEGREE SIGN} | ' + format(delta_pc, ".0f") + '%', linestyle = 'dashed')
ax.scatter(x_axis, y_axis, linewidth = 2, color = 'b')

ialt = -1
iinc1 = 0
pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
delta_pc = ( pc_max_here - pc_min_here ) / pc_nom * 100
y_axis = pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]
ax.plot(x_axis, y_axis, linewidth = 2, color = 'r', label =  str(alt1_arr[ialt]) + ' km | ' + str(np.abs(inc1_arr[iinc1] - inc0_arr[iinc0])) +  u'\N{DEGREE SIGN} | ' + format(delta_pc, ".0f") + '%') 
ax.scatter(x_axis, y_axis, linewidth = 2, color = 'r')
iinc1 = -1
pc_max_here = np.max(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_min_here = np.min(pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]) 
pc_nom = pc_median[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107] 
delta_pc = ( pc_max_here - pc_min_here ) / pc_nom * 100
y_axis = pc_arr[ialt, ialt, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, ibc, if107, :]
ax.plot(x_axis, y_axis, linewidth = 2, color = 'r', label =  str(alt1_arr[ialt]) + ' km | ' + str(np.abs(inc1_arr[iinc1] - inc0_arr[iinc0]))+  u'\N{DEGREE SIGN} | ' + format(delta_pc, ".0f") + '%', linestyle = 'dashed')
ax.scatter(x_axis, y_axis, linewidth = 2, color = 'r')

ax.plot([1,9], [1e-4, 1e-4], linewidth = 2, linestyle = 'dashed', color = 'k')
ax.text(5,1e-4,'Maneuver\nrecommended', color =  'k', horizontalalignment = 'center', verticalalignment = 'center', fontsize = fontsize_plot)

legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="        Altitude | " + r'$\Delta$inc. | ' + r'$\Delta$$P_c$' , fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

fig_save_name = 'two_ex_pc_vs_quantiles.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# alt0 = 400; alt1 = 400
# inc0 = 90; inc1 = 30
# arg_per0 = 0; arg_per1 = 0
# raan0 = 0; raan1 = 0
# true_ano0 = 0; true_ano1 = 0  
# ecc0 = 0; ecc1 = 0
# f10701 = 100;
# quantile01 = 3

# ialt0 = np.where( alt0_arr == alt0 )[0][0]
# ialt1 = np.where( alt1_arr == alt1 )[0][0]
# iinc0 = np.where( inc0_arr == inc0 )[0][0]
# iinc1 = np.where( inc1_arr == inc1 )[0][0]
# iarg_per0 = np.where( arg_per0_arr == arg_per0 )[0][0]
# iarg_per1 = np.where( arg_per1_arr == arg_per1 )[0][0]
# iraan0 = np.where( raan0_arr == raan0 )[0][0]
# iraan1 = np.where( raan1_arr == raan1 )[0][0]
# itrue_ano0 = np.where( true_ano0_arr == true_ano0 )[0][0]
# itrue_ano1 = np.where( true_ano1_arr == true_ano1 )[0][0]
# iecc0 = np.where( ecc0_arr == ecc0 )[0][0]
# iecc1 = np.where( ecc1_arr == ecc1 )[0][0]
# if107 = np.where( f107_arr == f10701 )[0][0]
# iquantile = np.where( quantile_arr == quantile01 )[0][0]

# print alt0, alt1, inc0, inc1, arg_per0, arg_per1, raan0, raan1, true_ano0, true_ano1, ecc0, ecc1, f10701, quantile01
# print pc_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile]
# print filename_run_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile]
# print tca_arr[ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile]
# print ialt0, ialt1, iinc0, iinc1, iarg_per0, iarg_per1, iraan0, iraan1, itrue_ano0, itrue_ano1, iecc0, iecc1, if107, iquantile






