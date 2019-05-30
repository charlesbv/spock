# PLOT PRNS VS TIME
import matplotlib.colors as colors
import matplotlib.patches as mpatches
plot_selected_prn = 1 # if this is set to 1 then 2 PRNs selected for the overpass will also be reported at the bottom of the plot
color_gain = ['grey', 'blue', 'limegreen', 'red']
label_gain = ['0', '2-5', '6-10', '11-15']
handles_arr = []
for icat in range(len(label_gain)):
    handles_arr.append(mpatches.Patch(color=color_gain[icat], label=label_gain[icat]))

marker_ant = [10, 11]
dant_netcdf = [0, 0.15]
dant_spock = [0, -0.15]

# figure out the number of different prn (SpOCK and on-board) during that time intervals
prn_list = []
for itime in range(itime_start, itime_stop):
    for ispec in range(4):
        if ( prn_selected[itime][ispec] in prn_list ) == False:
            prn_list.append(prn_selected[itime][ispec])

prn_list = np.array(prn_list)
nprn = len(prn_list)
prn_list_sort = prn_list[np.argsort(prn_list)]

height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25      
ratio_fig_size = 4./3
fig_title = ''#Accuracy VS RCG
y_label = 'PRN'
x_label = 'Real time'
ax_title =  date_entire[itime_start].replace('T', ' ') + ' to ' + date_entire[itime_stop][11:19] + ' UTC - sat-bop'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                                                                      
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])
ax.set_title(ax_title, weight = 'normal', fontsize  = fontsize_plot)
ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure                                       
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7,bottom=1,  left=1, right=1)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold                           
gps_satbop_value = []
for itime in range(itime_start, itime_stop):
    gain_spock_now = fom_selected[itime]
    gain_sort_index = np.argsort(-gain_spock_now) # descending order
    for ispec in range(4):
        prn_spock = prn_selected[itime][ispec]
        if which_ant_selected[itime][ispec] == 2: # starboard
            iant_spock = 0
        elif which_ant_selected[itime][ispec] == 3: # starboard
            iant_spock = 1
        prn_spock_value = np.where(prn_list_sort == prn_spock)[0][0]


        xaxis = (nb_seconds_since_start[itime] - nb_seconds_since_start[itime_start])/60.
        if fom_selected[itime][ispec] == 0: # !!!!!!!!!!! if change gain limmits then need to change also label_gain
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[0], s = 20)   
        elif ((fom_selected[itime][ispec] >= 1) & (fom_selected[itime][ispec] <= 5)):
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[1], s = 20)   
        elif ((fom_selected[itime][ispec] >= 6) & (fom_selected[itime][ispec] <= 10)):
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[2], s = 20)   
        elif ((fom_selected[itime][ispec] >= 11) & (fom_selected[itime][ispec] <= 15)):
            ax.scatter(xaxis, prn_spock_value + dant_spock[iant_spock] + 0.95,  marker = '.', color = color_gain[3], s = 20)   
        else:
            print "***! Error: the gain is not between 0 and 15. !***"; sys.exit()

tmax = (nb_seconds_since_start[itime_stop] - nb_seconds_since_start[itime_start])/60./2
tstar = tmax + 130/60.
ax.plot([tmax, tmax], [-0.6, -0.6 + 0.05], linewidth = 2, color = 'k')
ax.plot([tstar, tstar], [-0.6, -0.6 + 0.05], linewidth = 2, color = 'k')

dt_xlabel =  1. # min
xticks = np.arange(0, inter_dur_sec/60.+1, dt_xlabel)
date_list_str = []
date_list = [date_datetime_entire[itime_start] + timedelta(seconds=(x-xticks[0])*60.) for x in xticks]
for i in range(len(xticks)):
    date_list_str.append( str(date_list[i])[11:19] )
ax.xaxis.set_ticks(xticks)
ax.xaxis.set_ticklabels(date_list_str, fontsize = fontsize_plot, rotation=60, horizontalalignment = 'center')

ax.yaxis.set_ticks(np.arange(1, nprn+1))
ax.yaxis.set_ticklabels(prn_list_sort, fontsize = fontsize_plot)#, rotation='vertical')
ax.margins(0,0)
ax.set_ylim([-0.6, nprn+0.5])
ax.text(tmax, ax.get_ylim()[0], ' Tmax', rotation = 90, fontsize = fontsize_plot, horizontalalignment  = 'center', verticalalignment = 'bottom')
ax.text(tstar, ax.get_ylim()[0], ' Tmax+130s', rotation = 90, fontsize = fontsize_plot, horizontalalignment  = 'center', verticalalignment = 'bottom')
legend = ax.legend( loc='center left',  bbox_to_anchor=(1, 0.5), fontsize = fontsize_plot, handles=handles_arr, ncol=1, frameon=False, title = 'PRN gain')
legend.get_title().set_fontsize(str(fontsize_plot)) 

fig_save_name = '/Users/cbv/satbop_' + satbop_output_dir[:-1] + '_itimeStart' + str(itime_start) + '_itimeStop' + str(itime_stop) + '_prn_select.pdf'
#'testfm0' + str(cygfm)+ '.pdf'#+str(itime_in) + '_score_3d_binom.pdf'#time_diagram_prn_spock_onboard_iday' + str(idate) + '_itimeStart' + str(itime_start) + '_itimeStop' + +str(itime_stop) + '.pdf'
fig.savefig(fig_save_name, facecolor=fig  .get_facecolor(), edgecolor='none', bbox_inches='tight')
print fig_save_name
