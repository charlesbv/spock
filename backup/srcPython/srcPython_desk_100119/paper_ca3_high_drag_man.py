# This script looks at phase angle and phase rate between sc.
# ASSUMPTIONS:
# - all run have the same dates (start end, time step)

import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython")
from read_output_file import *
from read_input_file import *
import numpy as np
from find_in_read_input_order_variables import *
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import matplotlib.lines as mlines

same_sma = 0 # set this to 1 to plot the same sma runs. 0 otherwise

# runs to plot the phase and phase rate 
run_list =  ['nadir.txt', 'highdrag.txt']


nb_run = len(run_list)
var_to_read = ['phase_angle',  'nb_seconds_since_start', 'sma']
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_radius        = 6378.137; # mean equatorial radius (km)
for irun in range(nb_run):
    isc = 0
    input_filename = run_list[irun]
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
    dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    nb_steps = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
    nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')]; 
    if irun == 0: # see assumptions in header
        phase_angle = np.zeros([nb_run, nb_steps])
        phase_rate = np.zeros([nb_run, nb_steps])
        sma = np.zeros([nb_run, nb_steps])

    var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
    if irun == 0:
        date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        nb_seconds_since_start = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
    phase_angle[irun, :] = var_out[find_in_read_input_order_variables(var_out_order, 'phase_angle')]
    sma[irun, :] = var_out[find_in_read_input_order_variables(var_out_order, 'sma')] 
    # calculate phase rate from SMA (keplerian equations) because if from phase_angle problems every orbit
    mean_motion = np.sqrt( earth_mu / sma[irun, :]**3 ) # in s^-1
    phase_rate[irun, :] = mean_motion * 180. / np.pi # in deg/s




# PLOTS
## Parameters for the figure
height_fig = 15.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 28 
ratio_fig_size = 4./3


fig_title = ''#Phase angle difference as a function of time'

x_label = 'Time (days)'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
plt.rc('font', weight='bold') ## make the labels of the ticks in normal
gs = gridspec.GridSpec(3, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1)


## SMA of all runs
y_label = '$\Delta$SMA (m)'
ax_sma = fig.add_subplot(gs[0, 0])
ax_sma.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#ax_sma.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_sma.spines.itervalues()] # change the width of the frame of the figure
ax_sma.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in normal

x_axis = nb_seconds_since_start / 3600. / 24

ax_sma.plot(x_axis, (sma[1, :] - sma[0, :])*1000, linewidth = 2, color = 'k')

ax_sma.xaxis.set_ticklabels([], fontsize = fontsize_plot)#, rotation='vertical')
# ax_sma.set_xlim([0,10])
# ax_sma.set_ylim([0,.15])
ax_sma.margins(0,0)
# legend = ax_sma.legend(loc='upper right', bbox_to_anchor=(1, 1), numpoints = 1,  title="", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))
ax_sma.set_ylim([ax_sma.get_ylim()[0]-np.abs(0.1*ax_sma.get_ylim()[0]), ax_sma.get_ylim()[1]+np.abs(0.1*ax_sma.get_ylim()[0])])


## Phase angle difference 
y_label = '$\Phi$ (' + u'\N{DEGREE SIGN}' + ')'
ax_angle = fig.add_subplot(gs[1, 0])
ax_angle.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
#ax_angle.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_angle.spines.itervalues()] # change the width of the frame of the figure
ax_angle.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in normal

x_axis = nb_seconds_since_start / 3600. / 24


ax_angle.plot(x_axis, phase_angle[1, :] - phase_angle[0, :], linewidth = 2, color = 'k')

ax_angle.xaxis.set_ticklabels([], fontsize = fontsize_plot)#, rotation='vertical')

# if same_sma == 1:
#     ax_angle.set_ylim([9.8,15]) # same sma
# else:
#     ax_angle.set_ylim([0,30]) # different sma
ax_angle.margins(0,0)
# legend = ax_angle.legend(loc='lower right', bbox_to_anchor=(1, 0), numpoints = 1,  title="", fontsize = fontsize_plot,handles=[red_line,blue_line])
# legend.get_title().set_fontsize(str(fontsize_plot))
ax_angle.set_ylim([ax_angle.get_ylim()[0]-0.1*ax_angle.get_ylim()[1], ax_angle.get_ylim()[1]*1.1])


## Phase rate difference with 1st run in run_list
y_label = '$\Delta \dot \Phi$ (' + u'\N{DEGREE SIGN}' + '/day)'
ax_rate = fig.add_subplot(gs[2, 0])
ax_rate.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
ax_rate.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax_rate.spines.itervalues()] # change the width of the frame of the figure
ax_rate.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='bold') ## make the labels of the ticks in normal 

x_axis = nb_seconds_since_start / 3600. / 24

ax_rate.plot(x_axis, ( phase_rate[1, :] - phase_rate[0, :] ) * 3600. * 24, linewidth = 2, color = 'k', label = '500 km')


# if same_sma == 1:
#     ax_rate.set_ylim([-0.1,.7]) # same sma
# else:
#     ax_rate.set_ylim([-0.1,1.25]) # different sma
ax_rate.margins(0,0)
# legend = ax_rate.legend(loc='upper right', bbox_to_anchor=(1, 1), numpoints = 1,  title="", fontsize = fontsize_plot)
# legend.get_title().set_fontsize(str(fontsize_plot))

ax_rate.set_ylim([ax_rate.get_ylim()[0]-0.1*ax_rate.get_ylim()[1], ax_rate.get_ylim()[1]*1.1])

fig_save_name = 'dragman.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


