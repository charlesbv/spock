from get_name_mission import *
from cadre_read_last_tle import *
from get_prop_dir import *
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
import os
import subprocess
from read_collision_file import *


## NOTE 0: to run this script, you first need to run distance_ensemble_to_main_sc.py (with first_time = 1). This will create the data that plot_ensembles will then use to makes plots
## NOTE 1: to use this script, the only 3 parameters you have to set are:
## - if yes or no you want to save the plots (save_results = 1 if yes, 0 otherwise)
## - if yes or no you want to show the plots (show_plots = 1 if yes, 0 otherwise)
## - the name of the propagator main input file 
## - the path of the folder where you want to store the results (pickle, image, video): called 'path_folder_results'. In this folder, there must be the following subfolders: 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', and 'other'. In each of these subfolders, there must be the 2 subsubfolders: 'result', and 'pickle'. In the subsubfolder 'result', there must be the 2 subsubsubfolders: 'image', and 'video'.
## NOTE 2: this can be run if only ONE MAIN satellite was run (with ensembles of course)
## NOTE 3: the results (pickle, image, video) are stored in a subfolder of path_folder_results. The subfolder is either 'CYGNSS', 'CADRE', 'AERIE', 'SCION', 'QB50', or 'other'. This code figures out which folder the simulation corresponds to by reading the name of the output file chosen by the user in the main input file of the propagator (third line of section #SPACECRAFTS): it tries to find 'cygnss', 'cadre', 'aerie', 'scion', or 'qb50' in the name of the output file. If it does not find it, then the results here will be stored in path_folder_results/other/
## NOTE 4: to run this script, and any python script in the propagator, you need to be one subfolder deep from the main folder where the propagator runs are made. So if path_to_propagator/PropSim is the folder where the propagator runs are made, then the python scripts must be run for example in path_to_propagator/PropSim/subfolder_where_python_scipts_are_run


# SET THE PARAMETERS BELOW:
## Save or not the plots
save_results = 1

## Show or not the plots
show_plots = 1
## path of the folder where you want to store the results (pickle, image, video)
path_folder_results = '/raid3/Armada/Charles/python/' #get_prop_dir(2) + 'output/python_propagator/'

# end of SET THE PARAMETERS BELOW:

if show_plots == 1:
    plt.ion()
## main input file (argument in command line)
nb_run = len(sys.argv) - 2
colors = ['b','r','k', 'g', 'm', 'y','c']
tca_list = []
for irun in range(nb_run):
    main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2+irun]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'



    # read input file
    input_variables, order_input_variables = read_input_file(main_input_file_name)
    date_start = input_variables[0]; date_stop = input_variables[1]; dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6]; satellite_to_plot = input_variables[7]; nb_ensembles_coe = input_variables[12]; nb_ensembles_attitude = input_variables[13]; nb_ensembles_cd = input_variables[11]; ensemble_to_plot_temp = input_variables[14]; 
    nb_ensembles_density = input_variables[17]
    collision_filename = input_variables[18]
    n = nb_steps        

    ## Nb of ensembles
    nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd, nb_ensembles_density]
    nb_ensembles = np.max(nb_ensembles_array)
    for i in range(len(nb_ensembles_array)):
        if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
            nb_ensembles = nb_ensembles_array[i]

    name_mission = get_name_mission(satellite_to_plot[0])
    root_save_fig_name = path_folder_results + name_mission + '/result/image/' + satellite_to_plot[0][:-5]+ "_"

    # Read collision file
    date_collision, nb_collisions_each_dt, cpc, cpc_final, tca = read_collision_file( collision_filename )
    nb_ca = cpc.shape[0]
    nb_steps_collisions = cpc.shape[1]
    dt_collision = ( date_collision[0,1] - date_collision[0,0] ).total_seconds() 
    tca_list.append(tca)
    # Plot
    dt = dt_collision
    nb_steps = nb_steps_collisions
    ## Set up plot parameters 
    height_fig = 9.  # the width is calculated as height_fig * 4/3.
    fontsize_plot = 20 
    hour_time_step_xticks = 1. # time step of ticks when plotting a function as a function of time

    ## Make plots
    ratio_fig_size = 4./3
    tca = tca[0]### ONLY ONE TCA FOR THIS SCRIPT
    if irun == 0:
        earliest_tca = tca 
    delta_tca_with_earliest_tca = ( datetime.strptime(tca, "%Y-%m-%dT%H:%M:%S.%f")- datetime.strptime(earliest_tca, "%Y-%m-%dT%H:%M:%S.%f") ).total_seconds()
    if delta_tca_with_earliest_tca < 0:
        earliest_tca = tca

    tca_date = datetime.strptime(tca, "%Y-%m-%dT%H:%M:%S.%f")
    nb_seconds_span = dt * len(date_collision[0])  # LEO: nb_seconds_span/60 = 46.3 min
    nb_seconds_from_date_start_to_tca = ( tca_date - date_start ).total_seconds()
    nb_seconds_from_date_start_to_span_start = nb_seconds_from_date_start_to_tca - nb_seconds_span / 2
    nb_seconds_from_date_start_to_span_end = nb_seconds_from_date_start_to_tca + nb_seconds_span / 2

    nb_seconds_to_remove_from_span_at_start = 22 * 60. # can't be negative
    nb_seconds_to_remove_from_span_at_end = nb_seconds_to_remove_from_span_at_start

    nb_seconds_from_date_start_to_xaxis_start = nb_seconds_from_date_start_to_span_start +  nb_seconds_to_remove_from_span_at_start
    nb_seconds_from_date_start_to_xaxis_end = nb_seconds_from_date_start_to_span_end - nb_seconds_to_remove_from_span_at_end

    x_axis = np.arange(nb_seconds_from_date_start_to_xaxis_start, nb_seconds_from_date_start_to_xaxis_end , dt)


#    nb_seconds_from_start_span_to_tca = nb_seconds_span / 2  + delta_tca_with_earliest_tca # to plot the vertical line representing the tca, needs to take into account that the xaxis starts with respect to the EARLIEST tca



    ### Cumlative proability of collision for all close approaches
    fig_title = ''
    y_label = 'Cumul. probability'
    factor_on_y = 1.
    ### Plot with these parameters
    nb_pairs_ca = (int) ( round( nb_ca / 2. ) )
    for ipair in range(nb_pairs_ca): # plot two close approaches per figure ##!!!!!!!!!!!!actually workds only if one TCA
        if irun == 0:
            fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
            fig.suptitle(fig_title, y = 0.973,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
            plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        for ica in range(ipair*2,(ipair+1)*2):
            if ica < nb_ca:
                start_index = (int)( nb_seconds_to_remove_from_span_at_start / dt)
                end_index = start_index + len(x_axis)
                y_axis = cpc[ica, start_index:end_index]
                if irun == 0:
                    if ica == nb_steps_collisions - 1:
                        gs = gridspec.GridSpec(2-np.mod(nb_ca,2), 1)
                    else:
                        gs = gridspec.GridSpec(1 , 1) # !!!!!! used to be gridspec.GridSpec(2 , 1)
                    gs.update(left = 0.11, right=0.89, top = 0.96,bottom = 0.12, hspace = 0.45)
                    ax = fig.add_subplot(gs[np.mod(ica,2), 0])
                    ax.set_title('', weight = 'bold', fontsize  = fontsize_plot)
                    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)

                ax.plot(x_axis, y_axis * factor_on_y, linewidth = 2, color = colors[irun])
                ax.plot([nb_seconds_from_date_start_to_tca, nb_seconds_from_date_start_to_tca],[0,1], linewidth = 3, color = colors[irun],linestyle = 'dotted')
                if irun == 0:
                    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
                    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
                    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

                    # hour_time_step_xticks_in_seconds = hour_time_step_xticks * 3600 
                    # xticks = np.arange(0, nb_seconds_span, hour_time_step_xticks_in_seconds)
                    # date_list_str = []
                    # nb_hours_simu = nb_steps * dt/ 3600.
                    # date_list = [date_collision[ica,0] + timedelta(hours=x) for x in np.arange(0, nb_hours_simu+1, hour_time_step_xticks)]
                    # for i in range(len(xticks)):
                    #     if hour_time_step_xticks < 12:
                    #         if i == 0:
                    #             date_list_str.append("h+" +  str(int(xticks[i] / 3600.)) + "\n(TCA - T/2)") 
                    #         else:
                    #             date_list_str.append("+" + str(int(xticks[i]  / 3600.)))
                    #     else:
                    #         date_list_str.append( str(date_list[i])[5:10] + "\n(day + " + str(int(xticks[i] * step_plot / 24.)) + ")")
                    # ax.xaxis.set_ticks(xticks)
                    # ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot)#, rotation='vertical')
                    ax.xaxis.set_ticks([]) 
                    ax.set_xlabel('Real time', weight = 'bold', fontsize  = fontsize_plot) # the span here is the span of the run that is the first called in the line command
                if irun == 0:
                    ax.set_ylim([0, max(cpc_final)*1.1])
                if ((irun == 1) & (max(cpc_final) > ax.get_ylim()[1])):
                    ax.set_ylim([0, max(cpc_final)*1.1])
                ax.margins(0,1) #ax.margins(0,1)
                
                if irun == 0:
                    min_xaxis_min_of_all_run = nb_seconds_from_date_start_to_xaxis_start
                    max_xaxis_max_of_all_run = nb_seconds_from_date_start_to_xaxis_end
                    min_xaxis_max_of_all_run = nb_seconds_from_date_start_to_xaxis_start
                else:
                    if nb_seconds_from_date_start_to_xaxis_start < min_xaxis_min_of_all_run:
                        min_xaxis_min_of_all_run = nb_seconds_from_date_start_to_xaxis_start
                    if nb_seconds_from_date_start_to_xaxis_end > max_xaxis_max_of_all_run:
                        max_xaxis_max_of_all_run = nb_seconds_from_date_start_to_xaxis_end 
                    if  nb_seconds_from_date_start_to_xaxis_start > min_xaxis_max_of_all_run:
                        min_xaxis_max_of_all_run = nb_seconds_from_date_start_to_xaxis_start
                if irun == nb_run - 1:
                    ax.text(min_xaxis_min_of_all_run, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0])/ 20., str(date_start + timedelta( seconds = min_xaxis_min_of_all_run) )[11:-3], fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')
                    ax.text(max_xaxis_max_of_all_run - dt, ax.get_ylim()[0] - (ax.get_ylim()[1] - ax.get_ylim()[0])/ 20., str(date_start + timedelta( seconds = max_xaxis_max_of_all_run) )[11:-3], fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')
                ax.text(min_xaxis_max_of_all_run, ax.get_ylim()[1] - (ax.get_ylim()[1] - ax.get_ylim()[0])/ 20. * (irun+1), tca, color = colors[irun], fontsize = fontsize_plot, weight = 'bold')
                # if irun == 0:
                #     max_nb_seconds_from_date_start_to_xaxis_start = nb_seconds_from_date_start_to_xaxis_start
                #     min_nb_seconds_from_date_start_to_xaxis_start = nb_seconds_from_date_start_to_xaxis_start
                # if nb_seconds_from_date_start_to_xaxis_start > max_nb_seconds_from_date_start_to_xaxis_start:
                #     max_nb_seconds_from_date_start_to_xaxis_start = nb_seconds_from_date_start_to_xaxis_start
                #     date_end_xaxis = datetime.strptime(tca, "%Y-%m-%dT%H:%M:%S.%f") + timedelta(seconds = max_nb_seconds_from_date_start_to_xaxis_start)
                #     ax.text(max_nb_seconds_from_date_start_to_xaxis_start, ax.get_ylim()[0]- (ax.get_ylim()[1] - ax.get_ylim()[0])/ 10., str(date_end_xaxis)[11:-3], fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')

                # if nb_seconds_from_date_start_to_xaxis_start < min_nb_seconds_from_date_start_to_xaxis_start:
                #     min_nb_seconds_from_date_start_to_xaxis_start = nb_seconds_from_date_start_to_xaxis_start
                #     date_start_xaxis = datetime.strptime(tca, "%Y-%m-%dT%H:%M:%S.%f") - timedelta(seconds = min_nb_seconds_from_date_start_to_xaxis_start)
                #     ax.text(min_nb_seconds_from_date_start_to_xaxis_start, ax.get_ylim()[0]- (ax.get_ylim()[1] - ax.get_ylim()[0])/ 10., str(date_start_xaxis)[11:-3], fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')
                # 

        if save_results == 1:
            if irun == nb_run - 1:
                fig_save_name = 'cumulative_probability_ca_' + str(ipair*2) + '_and_' + str(ipair*2+1)
                fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
                fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
                os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission)

if show_plots == 1:
    plt.show(); plt.show()


# # Plot nb_collisions_each_dt
# ## Parameters of figure
# height_fig = 9.  # the width is calculated as height_fig * 4/3.
# fontsize_plot = 20 
# ## Plot the figure
# x_axis_nb_collisions_each_dt = np.zeros([nb_time_steps_collisions])
# date_start_collision = date_collision[0]
# for idt in range(nb_time_steps_collisions):
#     x_axis_nb_collisions_each_dt[idt] = ( date_collision[idt] - date_start_collision ).total_seconds() / 3600.

# ratio_fig_size = 4./3
# fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# gs = gridspec.GridSpec(1, 1)
# gs.update(left= 0.10, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

# ax2 = fig.add_subplot(gs[0, 0])
# ax2.set_title('Probability of Collision VS Time Since Epoch', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
# ax2.set_ylabel('Probability', fontsize = fontsize_plot, weight = 'bold')
# ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.plot(x_axis_nb_collisions_each_dt, nb_collisions_each_dt / nb_ensembles**2, 'k', linewidth = 2)
# #ax2.get_xaxis().get_major_formatter().set_useOffset(False)
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# ax2.set_xlabel('Time since Epoch (hours)', fontsize = fontsize_plot, weight = 'bold')

# ax2.legend(fontsize = fontsize_plot, loc = 3)

# if save_results == 1:
#     fig_save_name = 'probability_collision_at_each_dt'
#     fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#     os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

# if show_plots == 1:
#     plt.show(); plt.show()

# # Plot cpc
# ## Parameters of figure
# height_fig = 9.  # the width is calculated as height_fig * 4/3.
# fontsize_plot = 20 
# ## Plot the figure
# x_axis_cpc = np.zeros([nb_time_steps_collisions])
# date_start_collision = date_collision[0]
# for idt in range(nb_time_steps_collisions):
#     x_axis_cpc[idt] = ( date_collision[idt] - date_start_collision ).total_seconds() / 3600.

# ratio_fig_size = 4./3
# fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
# fig.suptitle('', y = 0.975,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# gs = gridspec.GridSpec(1, 1)
# gs.update(left= 0.10, right=0.97, top = 0.93,bottom = 0.08, hspace = 0.01)

# ax2 = fig.add_subplot(gs[0, 0])
# ax2.set_title('Cumulative Probability of Collision VS Time Since Epoch', weight = 'bold', fontsize = fontsize_plot,  y = 1.008) 
# ax2.set_ylabel('Cumulative Probability', fontsize = fontsize_plot, weight = 'bold')
# ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) # set a bunch of parameters for the ticks and their labels
# [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
# ax2.plot(x_axis_cpc, cpc, 'k', linewidth = 2)
# #ax2.get_xaxis().get_major_formatter().set_useOffset(False)
# ax2.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)
# ax2.set_xlabel('Time since Epoch (hours)', fontsize = fontsize_plot, weight = 'bold')

# ax2.legend(fontsize = fontsize_plot, loc = 3)

# if save_results == 1:
#     fig_save_name = 'cumulative_probability_collision'
#     fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
#     os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' + name_subfolder_save)

# if show_plots == 1:
#     plt.show(); plt.show()


