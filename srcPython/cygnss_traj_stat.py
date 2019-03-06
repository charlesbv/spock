# This script was made to answer Scott's email on Dec 4, 2018.
# The script:
# 1- reads the position and velocity (ECEF r/v) of a given FM between date_start_str and date_stop_str
# 2- converts ECEF r/v to ECI r/v
# 3- converts ECI r/v to osculating orbital elements
# 4- for each orbital element, compute distributions
# INPUTS:
# - cygfm: which FM to look at
# - date_start: start date of the simulation
# - date_stop: stop date of the simulation

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
cygfm = 2
date_start_str = '2018-08-18' # YYYY-MM-DD
date_stop_str = '2018-08-20' # YYYY-MM-DD
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
from datetime import datetime, timedelta
import os
from cygnss_read_netcdf_and_convert_to_eci_and_oe import *
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches

def plot_dist(var_name, unit, x, y):
    y_label = 'Percentage'
    x_label = var_name + unit
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_title('Distribution of ' + var_name + ' from ' + date_start_str + ' to ' + date_stop_str  + ' (FM0' + cygfm_str+ ')', weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    binsize_actual = x[1] - x[0]
    ax.bar(x, y, binsize_actual)
    if var_name.lower() == 'sma':
        ax.plot([sma_thr_min, sma_thr_min], [np.min(y), np.max(y)], linewidth = 2, linestyle = 'dashed', color = 'r')
    ax.set_xlim([np.min(x), np.max(x)])
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
    fig.set_figheight(height_fig)
    fig.set_figwidth(height_fig*ratio_fig_size)
    fig_save_name = 'dist_' + var_name.replace(' ','_').lower() + '_' + date_start_str.replace('-','').replace(':','') + '_to_' + date_stop_str.replace('-','').replace(':','') + '_FM0' + cygfm_str + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')

def plot_var(var_name, unit, y, date_plot_start, date_plot_stop):
    # Plot the sma over a day
    date_plot_start_date =  datetime.strptime(date_plot_start, '%Y-%m-%dT%H:%M:%S')
    date_plot_stop_date =  datetime.strptime(date_plot_stop, '%Y-%m-%dT%H:%M:%S')
    istart = np.where(date_conc == date_plot_start_date)[0][0]
    istop = np.where(date_conc == date_plot_stop_date)[0][0] + 1
    y_label = var_name + unit
    x_label = 'Real time'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.17)
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ntime = istop - istart 
    nb_seconds_since_start_day = np.zeros([ntime])
    for itime in range(ntime):
        nb_seconds_since_start_day[itime] = (date_conc[istart + itime] - date_conc[istart]).total_seconds()
    ax.scatter(nb_seconds_since_start_day, y[istart:istop], linewidth = 2,color = 'b')
    ax.plot(nb_seconds_since_start_day, y[istart:istop], linewidth = 2,color = 'b')
    #in red are the ouliers
    ax.scatter(nb_seconds_since_start_day[sma_out_i - istart], y[sma_out_i], linewidth = 2,color = 'r', zorder = 5)
    label = ['Valid','Anomaly']
    color = ['blue', 'red']
    handles_arr = []
    for icat in range(len(label)):
        handles_arr.append(mpatches.Patch(color=color[icat], label=label[icat]))

    ax.legend( loc='center left',  bbox_to_anchor=(0, 0.93), fontsize = fontsize_plot, handles=handles_arr
               , frameon=False)#, columnspacing = -4.45, title =)    
    dt_xlabel =  5*60 # dt for ticks on x axis (in seconds)
    xticks = np.arange(0, nb_seconds_since_start_day[-1], dt_xlabel)
    date_list_str = []
    date_list = [date_plot_start_date + timedelta(seconds=x-xticks[0]) for x in xticks]
    for i in range(len(xticks)):
        date_list_str.append( str(date_list[i])[11:16] )
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot, rotation='vertical')
    ax.set_title(var_name + ' from ' + date_plot_start + ' to ' + date_plot_stop  + ' (FM0' + cygfm_str+ ')', weight = 'normal', fontsize  = fontsize_plot)
    ax.margins(0,0);
    ax.set_ylim(np.min(y[istart:istop]*0.9), np.max(y[istart:istop]*1.1))
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
    fig.set_figheight(height_fig)
    fig.set_figwidth(height_fig*ratio_fig_size)
    fig_save_name = var_name.lower() + '_' + date_plot_start.replace('-','').replace(':','') + '_to_' + date_plot_stop.replace('-','').replace(':','') + '_FM0' + cygfm_str + '.pdf' 
    fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')

    
cygfm_str = str(cygfm)
netcdf_dir = '/Users/cbv/cygnss/netcdf'
if netcdf_dir[-1] != '/':
    netcdf_dir = netcdf_dir + '/'

date_start = datetime.strptime(date_start_str, '%Y-%m-%d')
date_stop = datetime.strptime(date_stop_str, '%Y-%m-%d')
nb_day = (int)((date_stop - date_start).total_seconds()/3600./24)
sma = []; eccentricity = []; inclination = []; long_an = []; w = []; f = []; period = []; phase_angle = []
date_flight_rounded_date = []
for iday in range(nb_day):
    print 'day ' + str(iday) + ' out of ' + str(nb_day) + ' day(s)'
    date = date_start + timedelta(days = iday)
    yy = datetime.strftime(date, '%Y')
    doy = datetime.strftime(date, '%j')
    # Read netcdf file to get the ECEF r/v
    netcdf_path = netcdf_dir + yy + '/' + doy + '/'
    filename = [filename for filename in os.listdir(netcdf_path) if filename.startswith('cyg0' + cygfm_str)][0]
    netcdf_filename = netcdf_path + filename
    if iday == 0:
        load_spice_override = 1
    else:
        load_spice_override = 0
    date_flight_rounded_temp, lon_cyg_temp, lat_cyg_temp, lon_spec_temp, lat_spec_temp, fom_temp, gps_temp,\
        x_cyg_temp, y_cyg_temp, z_cyg_temp, vx_cyg_temp, vy_cyg_temp, vz_cyg_temp,date_flight_rounded_date_temp,\
        r_eci_cyg_temp, v_eci_cyg_temp, \
        sma_temp, eccentricity_temp, inclination_temp, long_an_temp, w_temp, \
        f_temp, period_temp, phase_angle_temp \
        = cygnss_read_netcdf_and_convert_to_eci_and_oe(netcdf_filename, load_spice_override)
    sma.append(sma_temp); eccentricity.append(eccentricity_temp); inclination.append(inclination_temp);
    long_an.append(long_an_temp); w.append(w_temp); f.append(f_temp); period.append(period_temp);
    phase_angle.append(phase_angle_temp)
    date_flight_rounded_date.append(date_flight_rounded_date_temp)

# Concatenate all days
sma_conc = []; eccentricity_conc = []; inclination_conc = []; long_an_conc = [];
w_conc = []; f_conc = []; period_conc = []; phase_angle_conc = []
date_conc = []
for iday in range(nb_day):
    sma_conc = sma_conc + sma[iday]; eccentricity_conc = eccentricity_conc + eccentricity[iday]
    inclination_conc = inclination_conc + inclination[iday]
    long_an_conc = long_an_conc + long_an[iday]; w_conc = w_conc + w[iday]
    f_conc = f_conc + f[iday]; period_conc = period_conc + period[iday]
    phase_angle_conc = phase_angle_conc + phase_angle[iday]
    date_conc = date_conc + date_flight_rounded_date[iday]
sma_conc = np.array(sma_conc); eccentricity_conc = np.array(eccentricity_conc); inclination_conc = np.array(inclination_conc); long_an_conc = np.array(long_an_conc);
w_conc = np.array(w_conc); f_conc = np.array(f_conc); period_conc = np.array(period_conc); phase_angle_conc = np.array(phase_angle_conc)
date_conc = np.array(date_conc)

    
    
sma_quant = np.zeros([5])
sma_quant[0] = np.percentile(sma_conc, 10); sma_quant[1] = np.percentile(sma_conc, 25)
sma_quant[2] = np.median(sma_conc); sma_quant[3] = np.percentile(sma_conc, 75);
sma_quant[4] = np.percentile(sma_conc, 90) 

sma_thr_min = sma_quant[0] - 100*(sma_quant[4] - sma_quant[0])
sma_thr_max = sma_quant[4] + 100*(sma_quant[4] - sma_quant[0])

# Compute distributions of the orbital elements for the entire period [day_start, day_stop]
height_fig = 15.  # the width is calculated as height_fig * 4/3. 
fontsize_plot = 25
ratio_fig_size = 4./3
fig_title = ''#VLLH distribution specular points for a particular ispec at a particular time
y_label = 'Percentage'

# Dates below are not for the distributions but for the plot that show the variables (e.g., sma) as a function of time
date_plot_start = '2018-08-18T16:00:00'
date_plot_stop = '2018-08-18T17:55:00'


## SMA
min_bin = np.min(sma_conc)#6898.
max_bin = np.max(sma_conc)#6906.
nbin = 10
bin_size = (max_bin - min_bin)/nbin
dist_sma_data = np.histogram(sma_conc, bins = np.arange(min_bin, max_bin + bin_size, bin_size))
bin_array_temp = dist_sma_data[1]
bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
dist_sma = dist_sma_data[0] * 100. / len(sma_conc)
var_name = 'SMA'
unit = ' (km)'
x = bin_array
y = dist_sma
plot_dist(var_name, unit, x, y)


sma_out_i = np.where((sma_conc < sma_thr_min) | (sma_conc > sma_thr_max))[0]

y = sma_conc
plot_var(var_name, unit, y, date_plot_start, date_plot_stop)

## eccentricity
dist_eccentricity_data = np.histogram(eccentricity_conc, bins = 10)
bin_array_temp = dist_eccentricity_data[1]
bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
dist_eccentricity = dist_eccentricity_data[0] * 100. / len(eccentricity_conc)
var_name = 'Eccentricity'
unit = ''
x = bin_array
y = dist_eccentricity
plot_dist(var_name, unit, x, y)
y = eccentricity_conc
plot_var(var_name, unit, y, date_plot_start, date_plot_stop)

## inclination
dist_inclination_data = np.histogram(inclination_conc, bins = 10)
bin_array_temp = dist_inclination_data[1]
bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
dist_inclination = dist_inclination_data[0] * 100. / len(inclination_conc)
var_name = 'Inclination'
unit =  u' (\N{DEGREE SIGN})'
x = bin_array
y = dist_inclination
plot_dist(var_name, unit, x, y)
y = inclination_conc
plot_var(var_name, unit, y, date_plot_start, date_plot_stop)


# ## long_an
# dist_long_an_data = np.histogram(long_an_conc, bins = 10)
# bin_array_temp = dist_long_an_data[1]
# bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
# dist_long_an = dist_long_an_data[0] * 100. / len(long_an_conc)
# var_name = 'RAAN'
# unit =  u' (\N{DEGREE SIGN})'
# x = bin_array
# y = dist_long_an
# plot_dist(var_name, unit, x, y)

# ## w
# dist_w_data = np.histogram(w_conc, bins = 10)
# bin_array_temp = dist_w_data[1]
# bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
# dist_w = dist_w_data[0] * 100. / len(w_conc)
# var_name = 'Argument of perigee'
# unit =  u' (\N{DEGREE SIGN})'
# x = bin_array
# y = dist_w
# plot_dist(var_name, unit, x, y)

# ## period
# dist_period_data = np.histogram(period_conc, bins = 10)
# bin_array_temp = dist_period_data[1]
# bin_array = ( bin_array_temp[:-1] + np.roll(bin_array_temp,-1)[:-1] ) /2.
# dist_period = dist_period_data[0] * 100. / len(period_conc)
# var_name = 'Period'
# unit =  ' (seconds)'
# x = bin_array
# y = dist_period
# plot_dist(var_name, unit, x, y)

