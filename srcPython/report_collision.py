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
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from seconds_since_epoch_start_to_utc import *
import sys
from read_input_file import *
from get_prop_dir import *
from find_in_read_input_order_variables import *
from datetime import datetime, timedelta
from read_collision_file import *

show_plot = 1

if show_plot == 1:
    plt.ion()

# arguments:
# 1: name of output report file (written by this script)
name_mission = 'other'
# Below set up the list of runs
list_run = ['devel_1126_ok_quartile_f107_1_quartile_ap_1.txt',
             '1126_ok_quartile_f107_2_quartile_ap_2.txt',
             '1126_ok_quartile_f107_3_quartile_ap_3.txt',
             '1126_ok_quartile_f107_4_quartile_ap_4.txt',
             '1126_ok_quartile_f107_5_quartile_ap_5.txt',
             '1126_ok_quartile_f107_6_quartile_ap_6.txt',
             '1126_ok_quartile_f107_7_quartile_ap_7.txt',
             '1126_ok_quartile_f107_8_quartile_ap_8.txt',
             '1126_ok_quartile_f107_9_quartile_ap_9.txt']


report_filename = get_prop_dir(1) + "srcPython/collision/" + sys.argv[1]
report_file = open(report_filename, "w+")
run_dir = "run_collision/"
input_filename_list = []
nb_runs = len(list_run)
tca_array = [] 
cpc_final_array = []
tca_all = []
cpc_final_all = []
cpc_array = []
cpc_all = []
for irun in range(nb_runs):
    input_filename = get_prop_dir(1) + run_dir + 'input/main_input/' + list_run[irun]
    input_filename_list.append( input_filename.split('/')[-1]) 
    param_in, param_in_name = read_input_file(input_filename)
    output_file_path_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_path_list')]
    date_start = param_in[find_in_read_input_order_variables(param_in_name, 'date_start')]

    output_run_dir_name = output_file_path_list[0].split('/')[-3]
    collision_filename = get_prop_dir(1) + run_dir + 'output/' + output_run_dir_name + '/' + output_run_dir_name + '_collision.txt'
    date_collision, nb_collisions_each_dt, cpc, cpc_final, tca = read_collision_file( collision_filename )
    dt_coll = ( date_collision[0][1] - date_collision[0][0] ).total_seconds()
    nb_ca = cpc.shape[0]

    for ica in range(nb_ca):
        tca_all.append([datetime.strptime(tca[ica], "%Y-%m-%dT%H:%M:%S.%f"), irun])
        cpc_final_all.append([cpc_final[ica], irun])
        cpc_all.append([cpc[ica], [date_collision[ica][0], date_collision[ica][-1]]])
    if irun == 0:
        min_tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
        min_tca_seconds = (datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f") - date_start).total_seconds()
    if (datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f") - date_start).total_seconds() < min_tca_seconds:
        min_tca_seconds = (datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f") - date_start).total_seconds()
        min_tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
    if irun == 0:
        max_tca = datetime.strptime(tca[-1], "%Y-%m-%dT%H:%M:%S.%f")
        max_tca_seconds = (datetime.strptime(tca[-1], "%Y-%m-%dT%H:%M:%S.%f") - date_start).total_seconds()
    if (datetime.strptime(tca[-1], "%Y-%m-%dT%H:%M:%S.%f") - date_start).total_seconds() > max_tca_seconds:
        max_tca_seconds = (datetime.strptime(tca[-1], "%Y-%m-%dT%H:%M:%S.%f") - date_start).total_seconds()
        max_tca = datetime.strptime(tca[-1], "%Y-%m-%dT%H:%M:%S.%f")
    #tca_array.append([datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f") for t in tca])
    tca_array.append(tca)
    cpc_final_array.append(cpc_final)
    cpc_array.append(cpc)
    print >> report_file, "#" + input_filename_list[-1]
    for ica in range(nb_ca):
        print >> report_file, tca[ica].replace("T"," "), '{0:20f}'.format(cpc_final[ica]).replace(" ", "")
    print >> report_file, ""


tca_array = np.array(tca_array)
cpc_final_array = np.array(cpc_final_array)
cpc_array = np.array(cpc_array)

report_file.close()


tca_all_sorted = sorted(tca_all)
nb_ca_total = len(tca_all)
index_tca_all_sorted = [i[0] for i in sorted(enumerate(tca_all), key=lambda x:x[1])]
cpc_final_all_sorted_on_tca_all_sorted = []
cpc_all_sorted_on_tca_all_sorted = []
for i in range(nb_ca_total):
    cpc_final_all_sorted_on_tca_all_sorted.append(cpc_final_all[index_tca_all_sorted[i]])
    cpc_all_sorted_on_tca_all_sorted.append(cpc_all[index_tca_all_sorted[i]])

min_tca_list = []
max_tca_list = []
min_tca_list.append(tca_all_sorted[0][0])
tca_by_group = []
tca_by_group_sublist = []
cpc_final_by_group = []
cpc_final_by_group_sublist = []
cpc_by_group = []
cpc_by_group_sublist = []
for ica in range(nb_ca_total-1):
    tca_by_group_sublist.append(tca_all_sorted[ica])
    cpc_final_by_group_sublist.append(cpc_final_all_sorted_on_tca_all_sorted[ica])
    cpc_by_group_sublist.append(cpc_all_sorted_on_tca_all_sorted[ica])
    if (tca_all_sorted[ica+1][0] - tca_all_sorted[ica][0]).total_seconds() > 300:
        min_tca_list.append(tca_all_sorted[ica+1][0])
        max_tca_list.append(tca_all_sorted[ica][0])
        tca_by_group.append(tca_by_group_sublist)
        tca_by_group_sublist = []
        cpc_final_by_group.append(cpc_final_by_group_sublist)
        cpc_final_by_group_sublist = []
        cpc_by_group.append(cpc_by_group_sublist)
        cpc_by_group_sublist = []

if len(max_tca_list) > 0:
    if max_tca_list[-1] != tca_all_sorted[-1][0]:
        max_tca_list.append(tca_all_sorted[-1][0])
        tca_by_group_sublist.append(tca_all_sorted[-1])
        tca_by_group.append(tca_by_group_sublist)
        cpc_final_by_group_sublist.append(cpc_final_all_sorted_on_tca_all_sorted[-1])
        cpc_final_by_group.append(cpc_final_by_group_sublist)
        cpc_by_group_sublist.append(cpc_all_sorted_on_tca_all_sorted[-1])
        cpc_by_group.append(cpc_by_group_sublist)
else:
    max_tca_list.append(tca_all_sorted[-1][0])
    tca_by_group_sublist.append(tca_all_sorted[-1])
    tca_by_group.append(tca_by_group_sublist)
    cpc_final_by_group_sublist.append(cpc_final_all_sorted_on_tca_all_sorted[-1])
    cpc_final_by_group.append(cpc_final_by_group_sublist)
    cpc_by_group_sublist.append(cpc_all_sorted_on_tca_all_sorted[-1])
    cpc_by_group.append(cpc_by_group_sublist)
    
nb_group_of_ca = len(min_tca_list)


# Plot
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9

# GENERATE DISCTINCT COLORS
NCURVES = nb_runs
np.random.seed(101)
curves = [np.random.random(20) for i in range(NCURVES)]
values = range(NCURVES)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

# END OF GENERATE DISCTINCT COLORS


# import colorsys
# N = nb_runs
# HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
# RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
# color_arr = ['k','b','r','g','m','y']
# marker_arr = ['o', 'v', '^', '<', '>', 's', '*','x','D']

error_pc_sampling = 0.0

if cpc_final_by_group[0][0][0] < 0:
    problem_in_cpc_final = 1 # for some reason, sometimes the cpc final written in the report collision file by SpOCK is messed up. In that case, we need to read the nb of collision per step and sum all of them


# PLOT CPC FINAL BY GROUP OF CA
min_pc_threshold_to_plot = 1*10**(-5) # any group of ca with a pc smaller than this threshold won't generate a figure. This is done because pc under these thrsholds are proabbly not accurate because the number of ensembles run for this simu are not enough
for igroup in range(nb_group_of_ca):
    nb_ca_for_this_group = len(cpc_final_by_group[igroup])
    if problem_in_cpc_final == 1:
        cpc_one_particular_group = np.array(cpc_by_group[igroup]) # this is kind of messed up but it's ok. The only things to know are that cpc_one_particular_group[:,0][ica] is the cpc of ica for this group, and cpc_one_particular_group[:,1][ica] is the start and end of span for this ica
        cpc_final_one_particular_group = np.zeros([nb_ca_for_this_group, 2])
        for ica in range(nb_ca_for_this_group):
            cpc_final_one_particular_group[ica, 0] = cpc_one_particular_group[:,0][ica][-1]
            cpc_final_one_particular_group[ica, 1] = cpc_final_by_group[igroup][ica][1]
    else:
        cpc_final_one_particular_group = np.array(cpc_final_by_group[igroup])

        
    tca_one_particular_group = np.array(tca_by_group[igroup])

    if min(cpc_final_one_particular_group[:,0]) > min_pc_threshold_to_plot:

        fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
        gs = gridspec.GridSpec(1, 1)
        gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
        fig.suptitle('', fontsize = 22)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        ax1 = fig.add_subplot(gs[0, 0])
        #ax1 = fig.add_subplot(111)
        ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
        ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
        #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

        nb_seconds_for_this_group = (max_tca_list[igroup] - min_tca_list[igroup]).total_seconds()

        date_ref = tca_one_particular_group[0, 0] # date_start
        seconds_from_date_ref_to_tca_one_particular_group = np.zeros([nb_ca_for_this_group])
        for ica in range(nb_ca_for_this_group):
            colorVal = scalarMap.to_rgba(values[tca_one_particular_group[ica, 1]])
            seconds_from_date_ref_to_tca_one_particular_group[ica] = ( tca_one_particular_group[ica, 0] - date_ref ).total_seconds()
            ax1.plot([seconds_from_date_ref_to_tca_one_particular_group[ica], seconds_from_date_ref_to_tca_one_particular_group[ica]], [cpc_final_one_particular_group[ica,0]*(1-error_pc_sampling/2.), cpc_final_one_particular_group[ica,0]*(1+error_pc_sampling/2.)], marker = 's', color = colorVal, label = input_filename_list[tca_one_particular_group[ica, 1]].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%" , markersize = 10,markeredgecolor = 'none',linestyle="None" ) #label = 'F10.7 = ' + input_filename_list[tca_one_particular_group[ica, 1]].split('f107_')[1].split('_ap')[0] + ' | Ap = ' + input_filename_list[tca_one_particular_group[ica, 1]].split('ap_')[1].split('.')[0])

        ax1.plot([min(seconds_from_date_ref_to_tca_one_particular_group),max(seconds_from_date_ref_to_tca_one_particular_group)],[0.0001, 0.0001], color = 'r', linewidth = 4, linestyle = 'dashed')
        ax1.text(max(seconds_from_date_ref_to_tca_one_particular_group),0.0001,'Maneuver\nrecommended', color = 'r', horizontalalignment = 'right', fontsize = fontsize_plot, verticalalignment = 'center')
        
    # lines = []
    # for idx in range(len(curves)):
    #     line = curves[idx]
    #     colorVal = scalarMap.to_rgba(values[idx])
    #     colorText = (
    #         'color: (%4.2f,%4.2f,%4.2f)'%(colorVal[0],colorVal[1],colorVal[2])
    #         )
    #     retLine, = ax.plot(line,
    #                        color=colorVal,
    #                        label=colorText)

    #    ax1.plot(seconds_from_date_ref_to_tca_one_particular_group, cpc_final_one_particular_group[:,0], color = 'k', linewidth = 2)



#
        ax1.set_xlim([ax1.get_xlim()[0] - ax1.get_xlim()[1]*0.05, ax1.get_xlim()[1]*1.05])
        ax1.set_ylim([min(cpc_final_one_particular_group[:,0])*0.99, max(cpc_final_one_particular_group[:,0])*1.01])
        ax1.margins(1,1) # autoscale both axes(fist value is for the x axis, second value for the y axis)
        legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quartile of\nF10.7/Ap dist.", fontsize = fontsize_plot)
        legend.get_title().set_fontsize(str(fontsize_plot))
        ax1.set_title('Probability of collision VS TCA for different solar activity conditions', weight = 'bold', fontsize = fontsize_plot , y = 1.005)  # Probability of collision VS TCA (group of CA #'+ str(igroup + 1) +')
        ax1.set_xlabel('Seconds since ' + str(date_ref), fontsize = fontsize_plot, weight = 'bold')
        ax1.set_ylabel('Pc', fontsize = fontsize_plot, weight = 'bold')


        fig_save_name = './collision/' + report_filename.split('/')[-1] + '_cpc_final_ca_group_' + str(igroup+1)+ '.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
        os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./" + name_mission)

#label = 'F10.7 = ' + input_filename_list[tca_one_particular_group[ica, 1]].split('f107_')[1].split('_ap')[0] + ' | Ap = ' + input_filename_list[tca_one_particular_group[ica, 1]].split('ap_')[1].split('.')[0]
#.split('f107')[1].split('.')[0]





#OLD
#no_error_1_4m_50000ens_quartile_f107_5_quartile_ap_5.txt no_error_1_4m_50000ens_quartile_f107_1_quartile_ap_1.txt no_error_1_4m_50000ens_quartile_f107_2_quartile_ap_2.txt no_error_1_4m_50000ens_quartile_f107_3_quartile_ap_3.txt no_error_1_4m_50000ens_quartile_f107_6_quartile_ap_6.txt no_error_1_4m_50000ens_quartile_f107_4_quartile_ap_4.txt no_error_1_4m_50000ens_quartile_f107_8_quartile_ap_8.txt no_error_1_4m_50000ens_quartile_f107_7_quartile_ap_7.txt no_error_1_4m_50000ens_quartile_f107_9_quartile_ap_9.txt

#ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_1_quartile_ap_1.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_2_quartile_ap_2.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_3_quartile_ap_3.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_4_quartile_ap_4.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_5_quartile_ap_5.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_6_quartile_ap_6.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_7_quartile_ap_7.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_8_quartile_ap_8.txt ca_10000m_gravity2_1_3m_1126_sc_100kg_55000ens_quartile_f107_9_quartile_ap_9.txt 

# ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_-0_25_sigma_ap_-0_25_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_-0_5_sigma_ap_-0_5_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_-0_75_sigma_ap_-0_75_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_-1_5_sigma_ap_-1_5_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_-1_sigma_ap_-1_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_-2_sigma_ap_-2_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_0_25_sigma_ap_0_25_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_0_5_sigma_ap_0_5_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_0_75_sigma_ap_0_75_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_0_sigma_ap_0_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_1_5_sigma_ap_1_5_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_1_sigma_ap_1_sigma.txt ca_5000m_gravity2_1_3m_1126_sc_100kg_sigma_f107_2_sigma_ap_2_sigma.txt  



# # PLOT CPC VS TIME BY GROUP OF CA
# min_pc_threshold_to_plot = 1*10**(-5) # any group of ca with a pc smaller than this threshold won't generate a figure. This is done because pc under these thrsholds are proabbly not accurate because the number of ensembles run for this simu are not enough
# for igroup in range(nb_group_of_ca):
#     cpc_one_particular_group = np.array(cpc_by_group[igroup]) # this is kind of messed up but it's ok. The only things to know are that cpc_one_particular_group[:,0][ica] is the cpc of ica for this group, and cpc_one_particular_group[:,1][ica] is the start and end of span for this ica
#     tca_one_particular_group = np.array(tca_by_group[igroup])

#     cpc_final_one_particular_group = np.array(cpc_final_by_group[igroup])
#     if min(cpc_final_one_particular_group[:,0]) > min_pc_threshold_to_plot:
#         fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
#         gs = gridspec.GridSpec(1, 1)
#         gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06) 
#         fig.suptitle('', fontsize = 22)
#         plt.rc('font', weight='bold') ## make the labels of the ticks in bold
#         ax1 = fig.add_subplot(gs[0, 0])
#         #ax1 = fig.add_subplot(111)
#         ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
#         ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#         #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
#         [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

#         nb_seconds_for_this_group = (max_tca_list[igroup] - min_tca_list[igroup]).total_seconds()

#         date_ref = tca_one_particular_group[0, 0] # date_start
#         nb_ca_for_this_group = cpc_final_one_particular_group.shape[0]
#         seconds_from_date_ref_to_tca_one_particular_group = np.zeros([nb_ca_for_this_group])
#         seconds_from_date_ref_to_start_span_one_particular_group = np.zeros([nb_ca_for_this_group])
#         nb_seconds_from_start_span_to_end_span_one_particular_group = np.zeros([nb_ca_for_this_group])
#         for ica in range(nb_ca_for_this_group):
#             nb_time_step_in_ca = len(cpc_one_particular_group[:,0][ica])
#             colorVal = scalarMap.to_rgba(values[tca_one_particular_group[ica, 1]])
#             seconds_from_date_ref_to_tca_one_particular_group[ica] = ( tca_one_particular_group[ica, 0] - date_ref ).total_seconds()
#             seconds_from_date_ref_to_start_span_one_particular_group[ica] = ( cpc_one_particular_group[:,1][ica][0] - date_ref ).total_seconds()
#             nb_seconds_from_start_span_to_end_span_one_particular_group[ica] = ( cpc_one_particular_group[:,1][ica][1] - cpc_one_particular_group[:,1][ica][0]  ).total_seconds()
#             xaxis = np.arange(seconds_from_date_ref_to_start_span_one_particular_group[ica] + dt_coll, seconds_from_date_ref_to_start_span_one_particular_group[ica]  + nb_seconds_from_start_span_to_end_span_one_particular_group[ica] + 2*dt_coll, dt_coll) # the + dt_coll is so that the number of collision for an interval (of length dt_coll) corresponds to theend of this interval and not the beginning of this is interval. This makes the unpertrubed TCA appear in the correct interval. If you don't understand this sentence, don't worry about it.
#             ax1.plot(xaxis, cpc_one_particular_group[:,0][ica], marker = '_', linewidth = 4, color = colorVal, label = input_filename_list[tca_one_particular_group[ica, 1]].split('f107_')[1].split('_ap')[0].replace("_sigma", "") + ' | '+ input_filename_list[tca_one_particular_group[ica, 1]].split('ap_')[1].split('.')[0].replace("_sigma", ""))
# # different label:           ax1.plot(xaxis, cpc_one_particular_group[:,0][ica], marker = '_', linewidth = 4, color = colorVal, label = 'F10.7 = ' + input_filename_list[tca_one_particular_group[ica, 1]].split('f107_')[1].split('_ap')[0] + ' | Ap = ' + input_filename_list[tca_one_particular_group[ica, 1]].split('ap_')[1].split('.')[0])



#         ax1.plot([min(seconds_from_date_ref_to_tca_one_particular_group),max(seconds_from_date_ref_to_tca_one_particular_group)],[0.0001, 0.0001], color = 'r', linewidth = 4, linestyle = 'dashed')
#     # lines = []
#     # for idx in range(len(curves)):
#     #     line = curves[idx]
#     #     colorVal = scalarMap.to_rgba(values[idx])
#     #     colorText = (
#     #         'color: (%4.2f,%4.2f,%4.2f)'%(colorVal[0],colorVal[1],colorVal[2])
#     #         )
#     #     retLine, = ax.plot(line,
#     #                        color=colorVal,
#     #                        label=colorText)

#     #    ax1.plot(seconds_from_date_ref_to_tca_one_particular_group, cpc_one_particular_group[:,0], color = 'k', linewidth = 2)


#         ax1.set_ylim([ax1.get_ylim()[0], ax1.get_ylim()[1]* (1 + ( max(cpc_final_one_particular_group[:,0])  -  min(cpc_final_one_particular_group[:,0] ) ) / 10. ) ])
#         ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)b
#         ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))#ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol = 3)
#         ax1.set_title('Cumulative probability of collision VS time', weight = 'bold', fontsize = fontsize_plot , y = 1.005)  # Probability of collision VS TCA (group of CA #'+ str(igroup + 1) +')
#         ax1.set_xlabel('Seconds since ' + str(date_ref) + '', fontsize = fontsize_plot, weight = 'bold')
#         ax1.set_ylabel('Cumul. probability', fontsize = fontsize_plot, weight = 'bold')


#         fig_save_name = './collision/' + report_filename.split('/')[-1] + '_cpc_vs_time_ca_group_' + str(igroup+1)+ '.pdf'
#         fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#         os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")

    

