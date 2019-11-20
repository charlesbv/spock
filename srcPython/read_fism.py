import scipy
from scipy import io
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
import numpy as np
import os
from datetime import datetime, timedelta


filename_list_not_sorted = []
date_file = []
dir = 'flare_data/2017'
for file in os.listdir(dir):
    if file.endswith(".sav"):
        filename_list_not_sorted.append(os.path.join(dir, file))
        date_file.append(file.split('_')[2].split('_')[0]) # FISM_60sec_2017251_v01_03
date_file = np.array(date_file)
sort_date_file = np.argsort(date_file)
filename_list = np.array(filename_list_not_sorted)[sort_date_file]
nfile = len(filename_list)
seconds_since_start_all = []
flux_at_wavelength_all = []
flux_integrated_all = []
wavelength_chosen = 154

for ifile in range(nfile): # !!!!!! uncomment
#ifile = 0 # !!!!!! remove and loop over all files

    filename= filename_list[ifile]
    dict = scipy.io.readsav(filename, idict=None, python_dict=False, uncompressed_file_name=None, verbose=False)

    date_str = dict('date')
    date = datetime.strptime(str(date_str), "%Y%j")
    utc = dict('utc')
    if ifile == 0:
        date_ref = date
    seconds_ref = (date - date_ref).total_seconds()
    seconds_since_start_all = seconds_since_start_all + list( utc + seconds_ref )

    flux = dict('fism_pred')
    wavelength = dict('fism_wv')
    iwave = np.where(wavelength == wavelength_chosen + 0.5)[0][0]
    flux_at_wavelength = list(flux[iwave, :])
    flux_at_wavelength_all = flux_at_wavelength_all + flux_at_wavelength

    ntime = flux.shape[1]
    flux_integrated = []
    for itime in range(ntime):
        where_not_nan_and_not_0 = np.where((np.isnan(flux[:, itime]) == False) & (flux[:, itime] != 0))[0]
        flux_integrated.append(np.sum(flux[where_not_nan_and_not_0, itime]))
    flux_integrated_all = flux_integrated_all + flux_integrated

# end of loop
seconds_since_start_all = np.array(seconds_since_start_all)
flux_at_wavelength_all = np.array(flux_at_wavelength_all)
flux_integrated_all = np.array(flux_integrated_all)
# FIGURE

height_fig = 11
ratio_fig_size = 4./3
fontsize_plot = 25

# Integrated flux
fig_title = ''#local time of perigee vs delta eccentricity
ax_title = 'FISM 1-minute resolution integrated flux'
x_label = 'Real time'
y_label = 'W/m$^2$'

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])


ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_title(ax_title, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in bold

x_axis = np.arange(0,1440)
stop_index = -12*60 
ax.plot(seconds_since_start_all[:stop_index]/3600., flux_integrated_all[:stop_index], linewidth = 2, color = 'k')

dt = 24 # in hour
date_arr = np.array([date_ref + timedelta(hours=i) for i in np.arange(0, seconds_since_start_all[-1]/3600.,dt)])
nticks = len(date_arr)
xticks = []
for itick in range(nticks):
    xticks.append((date_arr[itick] - date_ref).total_seconds()/3600)
date_list_str = []
date_list = [date_ref + timedelta(hours=x) for x in xticks]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[5:10])# + "\n" + str(date_list[i])[11:16] )

ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0)
fig_save_name = 'fig/fism_integrated.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
