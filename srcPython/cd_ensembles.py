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

import numpy as np
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import sys
import fileinput
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
plt.ion()

mu, sigma = 2.2, 0.2 # mean and standard deviation
nb_ensembles = 5000

s = np.random.normal(mu, sigma, nb_ensembles)



# Parameters of figures
height_fig = 9.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
when_plot_in_hour = 7*24 #5 + 50/60. # can be overwritten for each plot by uncommenting the same line in each plot section
index_when_plot = 0 # (int) (when_plot_in_hour * 3600L / dt)
step_std = 1./60 # step in hours to calculate the standard deviation
hour_time_step_xticks = 24. # time step of ticks when plotting a function as a function of time
ratio_fig_size = 4./3

fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle('', fontsize = (int)(fontsize_plot*1.1), weight = 'bold')
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left= 0.09, right=0.935, top = 0.93,bottom = 0.12)
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title('', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
ax1.set_ylabel('# of ensembles (out of ' + str(nb_ensembles)+')', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Cd', fontsize = fontsize_plot, weight = 'bold')
ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

n, bins, patches = ax1.hist(s , nb_ensembles / 10,  histtype='stepfilled', alpha = 0.7, color = 'b',label = '')

        #ax1.hist(s, nb_ensembles / 2, normed=True)

ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

plt.show(); plt.show()
fig_save_name = 'cd_ensemble.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + fig_save_name)
