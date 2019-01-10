# plot the lon lat from the files in toshare

import numpy as np
import matplotlib
import pickle
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import pickle
import fileinput
import time
from datetime import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import subprocess

from datetime import datetime, timedelta

filename = 'toshare/20180808T000000_to_20180808T235959.txt'

filespec = open(filename)
read_filespec = filespec.readlines()
n = len(read_filespec)
lon = np.zeros([n])
lat = np.zeros([n])
date = []
for i in range(n):
    date_temp  = read_filespec[i].split()[0]
    date_temp_date = datetime.strptime(date_temp, "%Y-%m-%dT%H:%M:%S")
    date.append(date_temp_date)
    lon[i] = np.float( read_filespec[i].split()[1] )
    lat[i] = np.float( read_filespec[i].split()[2] )

date = np.array(date)
date_start = date[0]
n_first_hours_val = 10 # in min
n_first_hours_val_sec = n_first_hours_val * 60.
index_n_first_hours = np.where(date < date[0] + timedelta(seconds = n_first_hours_val_sec))[0]
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 29
ratio_fig_size = 4./3
width_fig = 25

dindex = 10
icount = -1

min_lat_range = 0.
max_lat_range = 25.
min_lon_range = 115.
max_lon_range = 140.

for iii in np.arange(0,n,10):#np.arange(index_n_first_hours[0], index_n_first_hours[-1], dindex):
    icount = icount + 1
    print date[iii]
    fig_title = ''#Probability per $P_C$ bin'
    fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in normal
    gs = gridspec.GridSpec(1, 2)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)

    ax = fig.add_subplot(gs[0, 0])

    y_label = 'Latitude (degree)'
    x_label = 'Longitude (degree)'

    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='normal') ## make the labels of the ticks in normal   
    ax.text(0.5, 0.01, str(date[iii]), fontsize = fontsize_plot, transform = ax.transAxes, horizontalalignment = 'center' )
    #ax.scatter(lon[index_n_first_hours], lat[index_n_first_hours],  color = 'b', s = 10)
    ax.scatter(lon[iii], lat[iii],  color = 'b', s = 10)
    ax.set_xlim([min_lon_range, max_lon_range])
    ax.set_ylim([min_lat_range, max_lat_range])
    fig_save_name = 'ani/' + filename.split('/')[-1].replace(".txt", '_lon_lat_' + str(icount) + '.png')
    fig.set_figheight(height_fig)
    fig.set_figwidth(width_fig)

    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

os.system('ffmpeg -y -r 10 -i ani/' + filename.split('/')[-1].replace(".txt", '_lon_lat_%d.png') + ' -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p ani/' +  filename.split('/')[-1].replace(".txt", '_lon_lat.mp4'))
