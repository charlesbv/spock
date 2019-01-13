# This script is copy paste of kalman_6state.py adpated to 9 state. kalman.py is too old to be used (but has 9 elt in the state too).
import sys
sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/kalman/spock_development_new_structure_kalman_dev/srcPython")
from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from read_input_file import *
from read_output_file import *
from read_odtk import *

# How to run this script:
# python kalman.py arg
# arg can be:
# - 3d then plot 3d orbit
# - r then plot position on three graphs on same page
# - v then plot speed on three graphs on same page
# - ad
# - rho
# - Pdiag
# - dX
# - y
# - drag_sigma

# Assumptions:
# - the time in the measurement file must be the same as in the kalman file


# Read measurement file
## Read SpOCK main input file to figure out the name of the measurement file
input_filename = sys.argv[1]
var_in, var_in_order = read_input_file(input_filename)
filename_kalman_meas_in = var_in[find_in_read_input_order_variables(var_in_order, 'filename_kalman_meas_in')];
filename_kalman_meas_out = var_in[find_in_read_input_order_variables(var_in_order, 'filename_kalman_meas_out')];
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')];
mass = var_in[find_in_read_input_order_variables(var_in_order, 'mass')];
area = var_in[find_in_read_input_order_variables(var_in_order, 'area')];
cd = var_in[find_in_read_input_order_variables(var_in_order, 'cdd')];
date_start = var_in[find_in_read_input_order_variables(var_in_order, 'date_start')];
date_stop = var_in[find_in_read_input_order_variables(var_in_order, 'date_stop')];



date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S.%f")
date_stop_str = datetime.strftime(date_stop, "%Y-%m-%dT%H:%M:%S.%f")

# here just to read the cd and area from "_density" file -> no KF estimatino of Cd for this version of the KF that we use here
isc = 0
output_file_path_list_kalm = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list_kalm = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
var_to_read = ["cd" ,"tot_area_drag"]
var_out, var_out_order = read_output_file( output_file_path_list_kalm[isc] + output_file_name_list_kalm[isc], var_to_read )
cd_kalm = var_out[find_in_read_input_order_variables(var_out_order, 'cd')]
tot_area_drag_kalm = var_out[find_in_read_input_order_variables(var_out_order, 'tot_area_drag')]

# "spock/true/true1/noise_true1.txt" #"/Users/cbv/kalman/cbv/spock/netcdf/output/cyg07.ddmi.s20170609-000000-e20170609-235959.l1.power-brcs.sand004.txt"
file_meas = open(output_file_path_list[isc] + filename_kalman_meas_out[isc]) #filename_kalman_meas)
read_file_meas = file_meas.readlines()
nb_header = 0
# while (read_file_meas[nb_header].split()[0] != '#START'):
#     nb_header = nb_header + 1
# nb_header = nb_header + 1
n_meas = len(read_file_meas) - nb_header
r_meas = np.zeros([n_meas,3])
v_meas = np.zeros([n_meas,3])
date_meas_str = []
date_meas = []
nb_seconds_since_start_meas = []
n_meas_until_end_epoch = 0
itime = 0

while ( itime < n_meas ) :
    date_meas_str.append(read_file_meas[itime+nb_header].split()[0])
    if datetime.strptime(date_meas_str[-1], "%Y-%m-%dT%H:%M:%S.%f" ) > date_stop:
        break
    date_meas.append( datetime.strptime(date_meas_str[-1], "%Y-%m-%dT%H:%M:%S.%f" ) )
    nb_seconds_since_start_meas.append( (date_meas[-1] - date_meas[0]).total_seconds() )
    r_meas[itime,0] = read_file_meas[itime+nb_header].split()[1]
    r_meas[itime,1] = read_file_meas[itime+nb_header].split()[2]
    r_meas[itime,2] = read_file_meas[itime+nb_header].split()[3]
    v_meas[itime,0] = read_file_meas[itime+nb_header].split()[4]
    v_meas[itime,1] = read_file_meas[itime+nb_header].split()[5]
    v_meas[itime,2] = read_file_meas[itime+nb_header].split()[6]

    n_meas_until_end_epoch = n_meas_until_end_epoch + 1
    itime = itime + 1
n_meas = n_meas_until_end_epoch
date_meas_str = np.array(date_meas_str)
file_meas.close()

# Read Kalman file

#### Read SpOCK main input file to figure out stuff to then read the output
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
dt_kalm = var_in[find_in_read_input_order_variables(var_in_order, 'dt')]; # !!!be carfeul: dt_kalm is read from SpOCK main input file. hwover, the time step can be different (smaller) than this value if the observation time does not fall at an exact time step. On other wrods, the actual time step is not constant: it can be different from dt_kalm at each observation


## Read kalman filter output file
isc = 0 # !!!!!!! all sc
filename_kalm = output_file_path_list[isc] + "kalman_" + output_file_name_list[isc]
file_kalm = open(filename_kalm)
read_file_kalm = file_kalm.readlines()
nb_header = 0
# while (read_file_kalm[nb_header].split()[0] != '#START'):
#     nb_header = nb_header + 1
# nb_header = nb_header + 1
n_kalm = len(read_file_kalm) - nb_header
if n_kalm != n_meas:
    print "***! The number of measurements is different from the number of points in the Kalman filter output. !***";

r_kalm = np.zeros([n_kalm,3]) # ECI
v_kalm = np.zeros([n_kalm,3]) # ECI
a_kalm = np.zeros([n_kalm,3]) # 
a_kalm_total_drag_from_c = np.zeros([n_kalm,3]) # 
a_kalm_total_drag_from_c_mag = np.zeros([n_kalm]) # 
a_no_kalm_drag_from_c = np.zeros([n_kalm,3]) # 
a_no_kalm_drag_from_c_mag = np.zeros([n_kalm]) # 
Pdiag = np.zeros([n_kalm,7])
dX = np.zeros([n_kalm,7])
y = np.zeros([n_kalm,3])
sy = np.zeros([n_kalm,3])
y_pf = np.zeros([n_kalm,3])
ad_kalm = np.zeros([n_kalm])
rho_kalm = np.zeros([n_kalm])
rho_kalm_from_c = np.zeros([n_kalm])
drag_sigma_kalm = np.zeros([n_kalm])
tau = np.zeros([n_kalm])
sum_cd_a_cos = np.zeros([n_kalm])
is_obs = np.zeros([n_kalm]) # 1 if this time corresponds to a measurement, 0 if prediction only (the state is not estimated but only predicted (estimated is when we use the equations of the Kalman filter to use the measrument to get the estimate))
y_pf = np.zeros([n_kalm,6])
#cd_kalm = np.zeros([n_kalm])


date_kalm = []
date = []
nb_seconds_since_start_kalm = []
for itime in range(n_kalm):
    date_kalm.append(datetime.strptime( read_file_kalm[itime+nb_header].split()[0], "%Y-%m-%dT%H:%M:%S.%f" ))
    nb_seconds_since_start_kalm.append( (date_kalm[-1] - date_kalm[0]).total_seconds() )
    r_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[1]
    r_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[2]
    r_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[3]
    v_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[4]
    v_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[5]
    v_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[6]
    a_no_kalm_drag_from_c[itime,0] = read_file_kalm[itime+nb_header].split()[7] #only prediction 
    a_no_kalm_drag_from_c[itime,1] = read_file_kalm[itime+nb_header].split()[8]
    a_no_kalm_drag_from_c[itime,2] = read_file_kalm[itime+nb_header].split()[9]
    a_no_kalm_drag_from_c_mag[itime] = np.linalg.norm( a_no_kalm_drag_from_c[itime,:] )

    tau[itime] = read_file_kalm[itime+nb_header].split()[10]

    a_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[11] # only unknonwn acceleratoin
    a_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[12]
    a_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[13]

    rho_kalm_from_c[itime] = read_file_kalm[itime+nb_header].split()[14] 
    sum_cd_a_cos[itime] = read_file_kalm[itime+nb_header].split()[15] 
    if (read_file_kalm[itime+nb_header].split()[19] == 'obs'):
        is_obs[itime] = 1
        y_pf[itime,0] = read_file_kalm[itime+nb_header].split()[20]
        y_pf[itime,1] = read_file_kalm[itime+nb_header].split()[21]
        y_pf[itime,2] = read_file_kalm[itime+nb_header].split()[22]
        y_pf[itime,3] = read_file_kalm[itime+nb_header].split()[23]
        y_pf[itime,4] = read_file_kalm[itime+nb_header].split()[24]
        y_pf[itime,5] = read_file_kalm[itime+nb_header].split()[25]

    a_kalm_total_drag_from_c[itime, :] = a_kalm[itime,:] + a_no_kalm_drag_from_c[itime,:] # prediction + unknown = total drag
    a_kalm_total_drag_from_c_mag[itime] =  np.linalg.norm( a_kalm_total_drag_from_c[itime, :] )

where_obs = np.where(is_obs == 1)[0]
date_kalm = np.array(date_kalm)

# OLD (below it's when c code used to output much more things. I used this colde below for most of the time. To use it again, you need to uncomment the // OLD block in the c code (function kalman_write_out))
#     Pdiag[itime,0] = read_file_kalm[itime+nb_header].split()[7]
#     Pdiag[itime,1] = read_file_kalm[itime+nb_header].split()[8]
#     Pdiag[itime,2] = read_file_kalm[itime+nb_header].split()[9]
#     Pdiag[itime,3] = read_file_kalm[itime+nb_header].split()[10]
#     Pdiag[itime,4] = read_file_kalm[itime+nb_header].split()[11]
#     Pdiag[itime,5] = read_file_kalm[itime+nb_header].split()[12]
#     dX[itime,0] = read_file_kalm[itime+nb_header].split()[13]
#     dX[itime,1] = read_file_kalm[itime+nb_header].split()[14]
#     dX[itime,2] = read_file_kalm[itime+nb_header].split()[15]
#     dX[itime,3] = read_file_kalm[itime+nb_header].split()[16]
#     dX[itime,4] = read_file_kalm[itime+nb_header].split()[17]
#     dX[itime,5] = read_file_kalm[itime+nb_header].split()[18]
#     y[itime,0] = read_file_kalm[itime+nb_header].split()[19]
#     y[itime,1] = read_file_kalm[itime+nb_header].split()[20]
#     y[itime,2] = read_file_kalm[itime+nb_header].split()[21]
#     sy[itime,0] = read_file_kalm[itime+nb_header].split()[22]
#     sy[itime,1] = read_file_kalm[itime+nb_header].split()[23]
#     sy[itime,2] = read_file_kalm[itime+nb_header].split()[24]
#     y_pf[itime,0] = read_file_kalm[itime+nb_header].split()[25]
#     y_pf[itime,1] = read_file_kalm[itime+nb_header].split()[26]
#     y_pf[itime,2] = read_file_kalm[itime+nb_header].split()[27]
#     a_kalm[itime,0] = read_file_kalm[itime+nb_header].split()[28] # only unknonwn acceleratoin
#     a_kalm[itime,1] = read_file_kalm[itime+nb_header].split()[29]
#     a_kalm[itime,2] = read_file_kalm[itime+nb_header].split()[30]
#     ad_kalm[itime] = read_file_kalm[itime+nb_header].split()[31]
#     rho_kalm_from_c[itime] = read_file_kalm[itime+nb_header].split()[32] 
#     dX[itime,6] = read_file_kalm[itime+nb_header].split()[33]
#     Pdiag[itime,6] = read_file_kalm[itime+nb_header].split()[34]
#     drag_sigma_kalm[itime] = read_file_kalm[itime+nb_header].split()[35]
#     tau[itime] = read_file_kalm[itime+nb_header].split()[36]
#     a_no_kalm_drag_from_c[itime,0] = read_file_kalm[itime+nb_header].split()[37] #only prediction 
#     a_no_kalm_drag_from_c[itime,1] = read_file_kalm[itime+nb_header].split()[38]
#     a_no_kalm_drag_from_c[itime,2] = read_file_kalm[itime+nb_header].split()[39]
#     a_no_kalm_drag_from_c_mag[itime] = np.linalg.norm( a_no_kalm_drag_from_c[itime,:] )

#     a_kalm_total_drag_from_c[itime, :] = a_kalm[itime,:] + a_no_kalm_drag_from_c[itime,:] # prediction + unknown = total drag
#     a_kalm_total_drag_from_c_mag[itime] =  np.linalg.norm( a_kalm_total_drag_from_c[itime, :] )

# end of OLD

#    cd_kalm[itime] = read_file_kalm[itime+nb_header].split()[37]

file_kalm.close()

# Results with propagation only 
var_to_read = ["position","velocity", "acceleration_lvlh", "acceleration_lvlh_drag", "acceleration_eci_drag", "density", "nb_seconds_since_start"]
var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
v_spock = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
a_spock = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration_eci_drag')]
date_spock = var_out[find_in_read_input_order_variables(var_out_order, 'date')] # this date is different from the one from the kalman filter sicne the kalman filter also includes the observation times, which can potentially not be included in date_spock if the times of the observations do not fall exactly at the time step of SpOCK
rho_spock = var_out[find_in_read_input_order_variables(var_out_order, 'density')] 
nb_seconds_since_start_spock = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')] 
a_spock_mag_here = np.linalg.norm(a_spock, axis = 1)

# test if times are the same in measurement and kalman file
# if date_kalm != date_meas_str:
#     print "***! Times in the measurement file are different from times in the Kalman file. !***"

# Time axis (in hours if more than 3 hours, otherwise in minutes)
nb_seconds_since_start_kalm = np.array(nb_seconds_since_start_kalm)
if ( date_kalm[-1] - date_kalm[0] ).total_seconds() > 3*3600:
    x_unit = "hour"
    x_axis = nb_seconds_since_start_kalm/3600.
    x_axis_spock = nb_seconds_since_start_spock/3600.
else:
    x_unit = "min"
    x_axis = nb_seconds_since_start_kalm/60.
    x_axis_spock = nb_seconds_since_start_spock/60.




# Time axis_meas: always the same as measurement file because the time step is variable sometimes in the netcdf files for CYGNSS
nb_seconds_since_start_meas = np.array(nb_seconds_since_start_meas)
x_axis_meas = nb_seconds_since_start_meas / 3600.
index_plot_meas = 60 # !!!!!! careful: the time step in the measurement file is not always constant (particularly of the netcdf files). So make sure you know what you're doing here. In doubt, put 1




# What to plot
plot_3d = 0
plot_rv = 0
plot_v = 0
plot_Pdiag = 0
plot_dX = 0
plot_y = 0
plot_pred = 0
plot_rho = 0
if '3d' in sys.argv:
    plot_3d = 1
if 'rv' in sys.argv:
    plot_rv = 1
if 'Pdiag' in sys.argv:
    plot_Pdiag = 1
if 'dX' in sys.argv:
    plot_dX = 1
if 'y' in sys.argv:
    plot_y = 1
if 'rho' in sys.argv:
    plot_rho = 1
if 'adrag' in sys.argv:
    plot_adrag = 1
if 'y_pf' in sys.argv:
    plot_y_pf = 1


compare_out = ""
for i in range(2,len(sys.argv)):
    if '.txt' in sys.argv[i]:
        compare_out = sys.argv[i]
if 'pred' in sys.argv:
    plot_pred = 1
    
## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 20 
ratio_fig_size = 4./3

dt_plot = 60. # in seconds
if dt_plot <= dt_kalm: 
    index_plot = 1
else:
    index_plot = (int)(dt_plot / dt_kalm)




# Plot 3d orbits
if plot_3d == 1:
    plt.ioff()
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig_title = '3D orbits'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold

    ax = fig.gca(projection='3d')
    ax.plot(r_kalm[::index_plot,0], r_kalm[::index_plot,1], r_kalm[::index_plot,2], label='Kalman', color = 'r', linewidth = 4)
    ax.scatter(r_meas[::index_plot,0], r_meas[::index_plot,1], r_meas[::index_plot,2], label='Measurement', color = 'b', linewidth = 6)
    x_label = "X (km)"
    y_label = "Y (km)"
    z_label = "Z (km)"
    ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot, labelpad = 35)
    ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot, labelpad = 35)
    ax.set_zlabel(z_label, weight = 'normal', fontsize  = fontsize_plot, labelpad = 35)
    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)

    legend = ax.legend(loc = 'center', bbox_to_anchor=(0.5, 1),  fontsize = fontsize_plot)
    fig.savefig(output_file_name_list[isc].replace(".txt", "") + "_3d.pdf")
    plt.show()

# Plot position
if plot_rv == 1: # !!!!!! assumes easuremnt and kalman otuput file start at the same date
    fig_title = ''
    x_label = 'Time (h)'

    n_meas_plot = len(x_axis_meas[::index_plot_meas])
 
    date_meas_plot = []
    date_kalm_plot = []
    date_kalm_same_time_as_meas = date_kalm[where_obs]

    r_error_mag = np.zeros([n_meas_plot])
    
    r_kalm_same_time_as_meas = np.zeros([n_meas, 3])
    r_kalm_same_time_as_meas[:, 0] = r_kalm[where_obs, 0]
    r_kalm_same_time_as_meas[:, 1] = r_kalm[where_obs, 1]
    r_kalm_same_time_as_meas[:, 2] = r_kalm[where_obs, 2]

#     r_error_spock_mag = np.zeros([n_meas_plot])
#     r_spock_same_time_as_meas = np.zeros([n_meas, 3])
#     r_spock_same_time_as_meas[:, 0] = r_spock[where_obs, 0]
#     r_spock_same_time_as_meas[:, 1] = r_spock[where_obs, 1]
#     r_spock_same_time_as_meas[:, 2] = r_spock[where_obs, 2]

    v_error_mag = np.zeros([n_meas_plot])
    v_kalm_same_time_as_meas = np.zeros([n_meas, 3])
    v_kalm_same_time_as_meas[:, 0] = v_kalm[where_obs, 0]
    v_kalm_same_time_as_meas[:, 1] = v_kalm[where_obs, 1]
    v_kalm_same_time_as_meas[:, 2] = v_kalm[where_obs, 2]

#     v_error_spock_mag = np.zeros([n_meas_plot])
#     v_spock_same_time_as_meas = np.zeros([n_meas, 3])
#     v_spock_same_time_as_meas[:, 0] = v_spock[where_obs, 0]
#     v_spock_same_time_as_meas[:, 1] = v_spock[where_obs, 1]
#     v_spock_same_time_as_meas[:, 2] = v_spock[where_obs, 2]

    # Read ODTK output (Kyle Nave)
    filename_odtk = '/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/eclipse/cyg_data/F7_20170824_145117_STKdefpred_v001.e'
    r_odtk, v_odtk, date_odtk, date_odtk_raw = read_odtk(filename_odtk, date_start_str, date_stop_str)
    #cygfm_to_ccsds = ['F7','F9','2B','2C','2F','36','37','49']
    ## date_odtk and date_i


    for itime in range(n_meas_plot):
        date_meas_plot.append(date_meas[itime*index_plot_meas])
        date_kalm_plot.append(date_kalm_same_time_as_meas[itime*index_plot_meas])
        r_error_diffce = ( r_meas[itime*index_plot_meas,:] - r_kalm_same_time_as_meas[itime*index_plot_meas, :] )
        r_error_mag[itime] = np.linalg.norm( r_error_diffce )
#         r_error_diffce_spock = ( r_meas[itime*index_plot_meas,:] - r_spock_same_time_as_meas[itime*index_plot_meas, :] )
#         r_error_spock_mag[itime] = np.linalg.norm( r_error_diffce_spock )

        v_error_diffce = ( v_meas[itime*index_plot_meas,:] - v_kalm_same_time_as_meas[itime*index_plot_meas, :] )
        v_error_mag[itime] = np.linalg.norm( v_error_diffce )
#         v_error_diffce_spock = ( v_meas[itime*index_plot_meas,:] - v_spock_same_time_as_meas[itime*index_plot_meas, :] )
#         v_error_spock_mag[itime] = np.linalg.norm( v_error_diffce_spock )


    # Plot error VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    ## r
    y_label = '$|\mathrm{\Delta r}|$ (m)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('$|\mathrm{r_{kalm} - r_{true}}|$', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    
    ax1.plot(x_axis_meas[::index_plot_meas], (r_error_mag)*1000., linewidth = 2,label='Kalman', color = 'b')
#    ax1.plot(x_axis_meas[::index_plot_meas], (r_error_spock_mag)*1000., linewidth = 2,label='Prediction', color = 'r')
    ax1.margins(0,0)

    ## v
    y_label = '$|\mathrm{\Delta v}|$ (m/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('$|\mathrm{v_{kalm} - v_{true}}|$', weight = 'normal', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax2.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#    ax2.xaxis.set_ticklabels("")
    
    ax2.plot(x_axis_meas[::index_plot_meas], (v_error_mag)*1000., linewidth = 2,label='Kalman', color = 'b')
    #ax2.plot(x_axis_meas[::index_plot_meas], (v_error_spock_mag)*1000., linewidth = 2,label='Prediction', color = 'r')
    ax2.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_error_r_v.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Plot Pdiag
if plot_Pdiag == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    ## Position uncertainty
    y_label = 'sigma position (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Position uncertainty', fontsize = fontsize_plot, weight = 'normal')
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], Pdiag[::index_plot,0], linewidth = 2,label='X', color = 'b')
    ax1.plot(x_axis[::index_plot], Pdiag[::index_plot,1], linewidth = 2,label='Y', color = 'r')
    ax1.plot(x_axis[::index_plot], Pdiag[::index_plot,2], linewidth = 2,label='Z', color = 'k')
    legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Velocity uncertainty
    y_label = 'sigma velocity (km/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Velocity uncertainty', fontsize = fontsize_plot, weight = 'normal')
    ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], Pdiag[::index_plot,3], linewidth = 2,label='Vx', color = 'b')
    ax2.plot(x_axis[::index_plot], Pdiag[::index_plot,4], linewidth = 2,label='Vy', color = 'r')
    ax2.plot(x_axis[::index_plot], Pdiag[::index_plot,5], linewidth = 2,label='Vz', color = 'k')
    legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax2.margins(0,0)

    ## ? uncertainty
    y_label = '?'
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title('? uncertainty', fontsize = fontsize_plot, weight = 'normal')
    ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax3.plot(x_axis[::index_plot], Pdiag[::index_plot,6], linewidth = 2,label='?', color = 'b')
    ax3.margins(0,0)    

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_Pdiag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# Plot dX
if plot_dX == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    ## Position 
    y_label = 'dr (km)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Position', fontsize = fontsize_plot, weight = 'normal')
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], dX[::index_plot,0], linewidth = 2,label='dX', color = 'b')
    ax1.plot(x_axis[::index_plot], dX[::index_plot,1], linewidth = 2,label='dY', color = 'r')
    ax1.plot(x_axis[::index_plot], dX[::index_plot,2], linewidth = 2,label='dZ', color = 'k')
    legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.margins(0,0)

    ## Velocity 
    y_label = 'dv (km/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Velocity', fontsize = fontsize_plot, weight = 'normal')
    ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], dX[::index_plot,3], linewidth = 2,label='dVx', color = 'b')
    ax2.plot(x_axis[::index_plot], dX[::index_plot,4], linewidth = 2,label='dVy', color = 'r')
    ax2.plot(x_axis[::index_plot], dX[::index_plot,5], linewidth = 2,label='dVz', color = 'k')
    legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax2.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_dX.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# # Plot y, sy, y_pf
# if plot_y == 1:
#     fig_title = ''
#     if x_unit == "hour":
#         x_label = 'Time (hr)'
#     else:
#         x_label = 'Time (min)'

#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

#     fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(3, 1)
#     gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
#     ## Position uncertainty
#     y_label = 'y (km)'
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('y', fontsize = fontsize_plot, weight = 'normal')
#     ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax1.xaxis.set_ticklabels("")
#     ax1.plot(x_axis[::index_plot], y[::index_plot,0], linewidth = 2,label='Y[0]', color = 'b')
#     ax1.plot(x_axis[::index_plot], y[::index_plot,1], linewidth = 2,label='Y[1]', color = 'r')
#     ax1.plot(x_axis[::index_plot], y[::index_plot,2], linewidth = 2,label='Y[2]', color = 'k')

#     legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
#     ax1.margins(0,0)

#     ## Velocity uncertainty
#     y_label = 'sy (?)'
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('sy', fontsize = fontsize_plot, weight = 'normal')
#     ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax2.xaxis.set_ticklabels("")
#     ax2.plot(x_axis[::index_plot], sy[::index_plot,0], linewidth = 2,label='sy[0]', color = 'b')
#     ax2.plot(x_axis[::index_plot], sy[::index_plot,1], linewidth = 2,label='sy[1]', color = 'r')
#     ax2.plot(x_axis[::index_plot], sy[::index_plot,2], linewidth = 2,label='sz[2]', color = 'k')
#     legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)
#     ax2.margins(0,0)

#     ## Velocity uncertainty
#     y_label = 'y_pf (km)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('y_pf', fontsize = fontsize_plot, weight = 'normal')
#     ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax3.plot(x_axis[::index_plot], y_pf[::index_plot,0], linewidth = 2,label='y_pf[0]', color = 'b')
#     ax3.plot(x_axis[::index_plot], y_pf[::index_plot,1], linewidth = 2,label='y_pf[1]', color = 'r')
#     ax3.plot(x_axis[::index_plot], y_pf[::index_plot,2], linewidth = 2,label='y_pf[2]', color = 'k')
#     legend = ax3.legend(loc = 'upper right',  fontsize = fontsize_plot)
#     ax3.margins(0,0)
    
#     fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_y_sy_ypf.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



# Plot rho_kalm_from_c
if plot_rho == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    

    y_label = 'Mass density (kg/m$^3$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Mass density as a function of time', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    #    ax1.plot(x_axis[::index_plot], rho_kalm[::index_plot], linewidth = 2,label='Kalman', color = 'b')
    ax1.plot(x_axis[::index_plot], rho_kalm_from_c[::index_plot], linewidth = 2,label='SpOCK', color = 'b')
    ax1.plot(x_axis_spock[::index_plot], rho_spock[::index_plot], linewidth = 2,label='NRLMSIS00e', color = 'r')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)

    ax1.set_ylim([min(rho_spock[::index_plot])*0.9, max(rho_spock[::index_plot]) *1.7]) 
    ax1.margins(0,0)
    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_rho.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



# Compare drag acceleration estimated by Kalman (= a_no_kalm + a_gauss_markov) to the drag acceleration without a Kalman filter (from the main output file (which is written only if the kalman filter is not used))
if plot_adrag == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    y_label = 'a (m/s$^2$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Drag acceleration as a function of time', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], a_kalm_total_drag_from_c_mag[::index_plot]*1000, linewidth = 2,label='Estimate (MSIS + correction)', color = 'b')
    ax1.plot(x_axis_spock[::index_plot], a_spock_mag_here[::index_plot]*1000, linewidth = 2,label='Prediction (MSIS)', color = 'r')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    legend = ax1.legend(bbox_to_anchor=(0, 1),loc = 'upper left',  fontsize = fontsize_plot, ncol = 2)

    ax1.set_ylim([min(a_spock_mag_here[::index_plot]*1000)*0.9, max(a_spock_mag_here[::index_plot]*1000) *1.2]) 
    ax1.margins(0,0)
    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_adrag_mag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



if plot_y_pf == 1:
    fig_title = ''
    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'

    y_pf_pos = np.linalg.norm(y_pf[:,:3], axis = 1)
    y_pf_vel = np.linalg.norm(y_pf[:,3:6], axis = 1)

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(2, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.27)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold

    ax1 = fig.add_subplot(gs[0, 0])
    y_label = 'Difference in position (m)'
    ax1.set_title('Difference in position and velocity estimated - measurements', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[np.where(y_pf_pos>0)[0]], y_pf_pos[np.where(y_pf_pos>0)[0]]*1000., linewidth = 2,label='Estimate', color = 'b')
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)

    #ax1.set_ylim([min(y_pf_pos[np.where(y_pf_pos>0)[0]]*1000)*0.9, max(y_pf_pos[np.where(y_pf_pos>0)[0]]*1000) *1.7]) 
    ax1.set_ylim([0,80])
    ax1.margins(0,0)


    ax2 = fig.add_subplot(gs[1, 0])
    y_label = 'Difference in velocity (m/s)'
    ax2.set_title('', weight = 'normal', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax2.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[np.where(y_pf_vel>0)[0]][::index_plot], y_pf_vel[np.where(y_pf_vel>0)[0]][::index_plot]*1000., linewidth = 2,label='Estimate', color = 'b')
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #legend = ax2.legend(loc = 'upper right',  fontsize = fontsize_plot)

    ax2.set_ylim([min(y_pf_vel[np.where(y_pf_vel>0)[0]]*1000)*0.9, max(y_pf_vel[np.where(y_pf_vel>0)[0]]*1000) *1.7]) 
    ax2.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_y_pf.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  





if compare_out != "":
    var_in, var_in_order = read_input_file(compare_out)
    output_file_path_list_comp = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list_comp = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
    var_to_read = ["position","velocity", "acceleration_lvlh", "acceleration_lvlh_drag", "acceleration_eci_drag","density", "cd" ,"tot_area_drag"]
    var_out, var_out_order = read_output_file( output_file_path_list_comp[isc] + output_file_name_list_comp[isc], var_to_read )
    r_comp = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
    v_comp = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
    a_comp = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration_eci_drag')]
    cd_comp = var_out[find_in_read_input_order_variables(var_out_order, 'cd')]
    tot_area_drag_comp = var_out[find_in_read_input_order_variables(var_out_order, 'tot_area_drag')]
    rho_comp = var_out[find_in_read_input_order_variables(var_out_order, 'density')] 
    n_comp = len(a_comp)
    ad_comp = np.zeros([n_comp])
    for itime in range(n_comp):
        ad_comp[itime] = np.linalg.norm(a_comp[itime,:])

    if x_unit == "hour":
        x_label = 'Time (hr)'
    else:
        x_label = 'Time (min)'
    fig_title = ''
    # fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    # fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    # gs = gridspec.GridSpec(3, 1)
    # gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    # ## X
    # y_label = 'X (m)'
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Position - Kalman (b), True (r)', weight = 'normal', fontsize = fontsize_plot)
    # ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax1.xaxis.set_ticklabels("")
    # ax1.scatter(x_axis[::index_plot], ( r_kalm[::index_plot,0] - r_comp[::index_plot,0] ) * 1000, linewidth = 2,label='Kalman', color = 'b')
    # ax1.plot(x_axis[::index_plot], ( r_comp[::index_plot,0] - r_comp[::index_plot,0] ) * 1000, linewidth = 4,label='True', color = 'r')
    # #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    # ax1.margins(0,0)

    # ## Y
    # y_label = 'Y (m)'
    # ax2 = fig.add_subplot(gs[1, 0])
    # ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    # ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax2.xaxis.set_ticklabels("")
    # ax2.scatter(x_axis[::index_plot], ( r_kalm[::index_plot,1] - r_comp[::index_plot,1] ) * 1000, linewidth = 2,label='Kalman', color = 'b')
    # ax2.plot(x_axis[::index_plot], ( r_comp[::index_plot,1] - r_comp[::index_plot,1] ) * 1000, linewidth = 4,label='True', color = 'r')
    # ax2.margins(0,0)

    # ## X
    # y_label = 'Z (m)'
    # ax3 = fig.add_subplot(gs[2, 0])
    # ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    # ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    # ax3.scatter(x_axis[::index_plot], ( r_kalm[::index_plot,2] - r_comp[::index_plot,2] ) * 1000, linewidth = 2,label='Kalman', color = 'b')
    # ax3.plot(x_axis[::index_plot], ( r_comp[::index_plot,2] - r_comp[::index_plot,2] ) * 1000, linewidth = 4,label='True', color = 'r')
    # ax3.margins(0,0)    

    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_true_vs_kalm_r.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    
    # if x_unit == "hour":
    #     x_label = 'Time (hr)'
    # else:
    #     x_label = 'Time (min)'
    # fig_title = ''
    # fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    # fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    # gs = gridspec.GridSpec(3, 1)
    # gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.1, wspace = 0.35)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    # ## X
    # y_label = 'Vx (m/s)'
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Velocity - Kalman (b), True (r)', weight = 'normal', fontsize = fontsize_plot)
    # ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax1.xaxis.set_ticklabels("")
    # ax1.scatter(x_axis[::index_plot], (v_kalm[::index_plot,0] - v_comp[::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
    # ax1.plot(x_axis[::index_plot], (v_comp[::index_plot,0] - v_comp[::index_plot,0])*1000., linewidth = 4,label='True', color = 'r')
    # #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    # ax1.margins(0,0)

    # ## Y
    # y_label = 'Vy (m/s)'
    # ax2 = fig.add_subplot(gs[1, 0])
    # ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    # ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # ax2.xaxis.set_ticklabels("")
    # ax2.scatter(x_axis[::index_plot], (v_kalm[::index_plot,1] - v_comp[::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
    # ax2.plot(x_axis[::index_plot], (v_comp[::index_plot,1] - v_comp[::index_plot,1])*1000., linewidth = 4,label='True', color = 'r')
    # ax2.margins(0,0)

    # ## X
    # y_label = 'Vz (m/s)'
    # ax3 = fig.add_subplot(gs[2, 0])
    # ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
    # ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    # ax3.scatter(x_axis[::index_plot], (v_kalm[::index_plot,2] - v_comp[::index_plot,2])*1000., linewidth = 2,label='Kalman', color = 'b')
    # ax3.plot(x_axis[::index_plot], (v_comp[::index_plot,2] - v_comp[::index_plot,2])*1000., linewidth = 4,label='True', color = 'r')
    # ax3.margins(0,0)    

    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_true_vs_kalm_v.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Measurements vs true
    r_error_meas =  np.zeros([n_kalm, 3])
    r_error_meas_mag = np.zeros([n_kalm])
    v_error_meas =  np.zeros([n_kalm, 3])
    v_error_meas_mag = np.zeros([n_kalm])    
    # Kalman vs true
    r_error = np.zeros([n_kalm, 3])
    v_error = np.zeros([n_kalm, 3])
    a_error = np.zeros([n_kalm, 3])
    r_error_mag = np.zeros([n_kalm])
    v_error_mag = np.zeros([n_kalm])
    a_error_mag = np.zeros([n_kalm])
    a_kalm_mag = np.zeros([n_kalm])
    # propagation only VS true
    r_spock_error = np.zeros([n_kalm, 3])
    v_spock_error = np.zeros([n_kalm, 3])
    a_spock_error = np.zeros([n_kalm, 3])
    r_spock_error_mag = np.zeros([n_kalm])
    v_spock_error_mag = np.zeros([n_kalm])
    a_spock_error_mag = np.zeros([n_kalm])
    a_spock_mag = np.zeros([n_kalm])
    # True
    a_comp_mag = np.zeros([n_kalm])
    v_comp_mag = np.zeros([n_kalm])
    v_kalm_mag = np.zeros([n_kalm])
    omega_rot = 0.00007292158553; # rad/s
    v_kalm_rel_to_atmo = np.zeros([n_kalm, 3])
    v_kalm_rel_to_atmo_mag = np.zeros([n_kalm])
    v_comp_rel_to_atmo = np.zeros([n_comp, 3])
    v_comp_rel_to_atmo_mag = np.zeros([n_comp])

    for itime in range(n_kalm):
        # Measurements vs true
        r_error_meas[itime, :] = r_meas[itime,:] - r_comp[itime,:]
        v_error_meas[itime, :] = v_meas[itime,:] - v_comp[itime,:]
        r_error_meas_mag[itime] = np.linalg.norm(r_error_meas[itime, :])
        v_error_meas_mag[itime] = np.linalg.norm(v_error_meas[itime, :])

        # Kalman vs true
        r_error[itime, :] = r_kalm[itime,:] - r_comp[itime,:]
        v_error[itime, :] = v_kalm[itime,:] - v_comp[itime,:]
        a_error[itime, :] = a_kalm[itime,:]  - a_comp[itime,:] # !!!!!!!!! should it be: a_error[itime, :] = a_kalm[itime,:]  - a_comp[itime,:] or:  a_error[itime, :] = a_spock[itime,:] + a_kalm[itime,:]  - a_comp[itime,:]
        r_error_mag[itime] = np.linalg.norm(r_error[itime, :])
        v_error_mag[itime] = np.linalg.norm(v_error[itime, :])
        a_error_mag[itime] = np.linalg.norm(a_error[itime, :])
        v_kalm_mag[itime] =  np.linalg.norm(v_kalm[itime, :])
        
        # propagation only VS true
        r_spock_error[itime, :] = r_spock[itime,:] - r_comp[itime,:]
        v_spock_error[itime, :] = v_spock[itime,:] - v_comp[itime,:]
        a_spock_error[itime, :] = a_spock[itime,:] - a_comp[itime,:]
        a_spock_mag[itime] = np.linalg.norm(a_spock[itime,:])
        r_spock_error_mag[itime] = np.linalg.norm(r_spock_error[itime, :])
        v_spock_error_mag[itime] = np.linalg.norm(v_spock_error[itime, :])
        a_spock_error_mag[itime] = np.linalg.norm(a_spock_error[itime, :])
        
        # True
        a_comp_mag[itime] = np.linalg.norm(a_comp[itime, :])
        v_comp_mag[itime] = np.linalg.norm(v_comp[itime, :])

        a_kalm_mag[itime] = np.linalg.norm(a_kalm[itime,:] + a_spock[itime,:])

    
        v_kalm_rel_to_atmo[itime, 0] = v_kalm[itime, 0] + omega_rot * r_kalm[itime, 1];
        v_kalm_rel_to_atmo[itime, 1] = v_kalm[itime, 1] - omega_rot * r_kalm[itime, 0];
        v_kalm_rel_to_atmo[itime, 2] = v_kalm[itime, 2];
        v_kalm_rel_to_atmo_mag[itime] = np.linalg.norm(v_kalm_rel_to_atmo[itime,:])

        v_comp_rel_to_atmo[itime, 0] = v_comp[itime, 0] + omega_rot * r_comp[itime, 1];
        v_comp_rel_to_atmo[itime, 1] = v_comp[itime, 1] - omega_rot * r_comp[itime, 0];
        v_comp_rel_to_atmo[itime, 2] = v_comp[itime, 2];
        v_comp_rel_to_atmo_mag[itime] = np.linalg.norm(v_comp_rel_to_atmo[itime,:])
        
        rho_kalm[itime] =  a_kalm_mag[itime] / ( 1./2 * cd_kalm[itime] * tot_area_drag_kalm[itime] / mass * v_kalm_rel_to_atmo_mag[itime]**2 ) / ( 1000**3 ) # in kg/m^3 !!!!!!! one surface only

    # Measurements vs true
    r_rmse_meas = np.sqrt( np.mean(r_error_meas_mag**2) )
    v_rmse_meas = np.sqrt( np.mean(v_error_meas_mag**2) )

    # Kalman vs true
    r_rmse = np.sqrt( np.mean(r_error_mag**2) )
    v_rmse = np.sqrt( np.mean(v_error_mag**2) )
    a_rmse = np.sqrt( np.mean(a_error_mag**2) )
    a_nrmse = np.sqrt( np.mean(a_error_mag**2) / np.mean(a_comp_mag**2))
    ax_rmse = np.sqrt( np.mean(a_error[::index_plot,0]**2) )
    ay_rmse = np.sqrt( np.mean(a_error[::index_plot,1]**2) )
    az_rmse = np.sqrt( np.mean(a_error[::index_plot,2]**2) )
    ax_nrmse = np.sqrt( np.mean(a_error[::index_plot,0]**2) / np.mean(a_comp[::index_plot,0]**2) )
    ay_nrmse = np.sqrt( np.mean(a_error[::index_plot,1]**2) / np.mean(a_comp[::index_plot,1]**2) )
    az_nrmse = np.sqrt( np.mean(a_error[::index_plot,2]**2) / np.mean(a_comp[::index_plot,2]**2) )

    # Propagation only VS true
    r_spock_rmse = np.sqrt( np.mean(r_spock_error_mag**2) )
    v_spock_rmse = np.sqrt( np.mean(v_spock_error_mag**2) )
    a_spock_rmse = np.sqrt( np.mean(a_spock_error_mag**2) )
    a_spock_nrmse = np.sqrt( np.mean(a_spock_error_mag**2) / np.mean(a_comp_mag**2))
    ax_spock_rmse = np.sqrt( np.mean(a_spock_error[::index_plot,0]**2) )
    ay_spock_rmse = np.sqrt( np.mean(a_spock_error[::index_plot,1]**2) )
    az_spock_rmse = np.sqrt( np.mean(a_spock_error[::index_plot,2]**2) )
    ax_spock_nrmse = np.sqrt( np.mean(a_spock_error[::index_plot,0]**2) / np.mean(a_comp[::index_plot,0]**2) )
    ay_spock_nrmse = np.sqrt( np.mean(a_spock_error[::index_plot,1]**2) / np.mean(a_comp[::index_plot,1]**2) )
    az_spock_nrmse = np.sqrt( np.mean(a_spock_error[::index_plot,2]**2) / np.mean(a_comp[::index_plot,2]**2) )


    # Plot rho VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Mass density (kg/m$^3$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('Mass density as a function of time', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
#    ax1.plot(x_axis[::index_plot], rho_kalm[::index_plot], linewidth = 2,label='Kalman', color = 'b')
    ax1.plot(x_axis[::index_plot], rho_kalm_from_c[::index_plot], linewidth = 2,label='Kalman', color = 'b')
    ax1.plot(x_axis[::index_plot], rho_comp[::index_plot], linewidth = 2,label='True', color = 'r')
#    ax1.plot(x_axis[::index_plot], rho_spock[::index_plot], linewidth = 2,label='True', color = 'k')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.set_ylim([min(rho_comp[::index_plot])*0.9, max(rho_comp[::index_plot]) *1.1])
    ax1.margins(0,0)
    # print rho_comp
    # print rho_kalm
    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_rho.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # # Plot ad VS time
    # fig_title = ''
    # fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    # fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    # gs = gridspec.GridSpec(4, 1)
    # gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    # ## X
    # y_label = '$|\mathrm{a_{LVLH}}|$ (m/s$^2$)'
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Magnitude of LVLH acceleration as a function of time', weight = 'normal', fontsize = fontsize_plot)
    # ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # #ax1.xaxis.set_ticklabels("")
    # ax1.scatter(x_axis[::index_plot], ad_kalm[::index_plot]*1000, linewidth = 2,label='Kalman', color = 'b')
    # ax1.scatter(x_axis[::index_plot], a_comp[::index_plot,0]*1000, linewidth = 2,label='True', color = 'r')
    # #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    # ax1.set_ylim([min([min(ad_kalm), min(ad_comp)]), max([max(ad_kalm), max(ad_comp)])])
    # ax1.margins(0,0)

    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_ad.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # #raise Exception

    # Plot error VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(4, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    ## X
    y_label = 'Error (m)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('$|\mathrm{r_{kalm} - r_{true}}|$', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[::index_plot], (r_error_mag[::index_plot])*1000., linewidth = 2,label='Kalman', color = 'b')
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
    ax1.text( max(x_axis[::index_plot]), max((r_error_mag[::index_plot])*1000.) - ( max((r_error_mag[::index_plot])*1000.) - min((r_error_mag[::index_plot])*1000.) )/7, "$\mathrm{RMSE_{kalm}}$ = " + '{:.2f}'.format(r_rmse*1000) + " m", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right')
    ax1.text( max(x_axis[::index_plot]), max((r_error_mag[::index_plot])*1000.) - 2* ( max((r_error_mag[::index_plot])*1000.) - min((r_error_mag[::index_plot])*1000.) )/7, "$\mathrm{RMSE_{meas}}$ = " + '{:.2f}'.format(r_rmse_meas*1000) + " m", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right')
    ax1.text( max(x_axis[::index_plot]), max((r_error_mag[::index_plot])*1000.) - 3* ( max((r_error_mag[::index_plot])*1000.) - min((r_error_mag[::index_plot])*1000.) )/7, "$\mathrm{RMSE_{prop}}$ = " + '{:.2f}'.format(r_spock_rmse*1000) + " m", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right')
    ax1.margins(0,0)

    ## Y
    y_label = 'Error (cm/s)'
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('$|\mathrm{v_{kalm} - v_{true}}|$', weight = 'normal', fontsize = fontsize_plot)
    ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax2.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
    ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax2.xaxis.set_ticklabels("")
    ax2.plot(x_axis[::index_plot], (v_error_mag[::index_plot])*100000., linewidth = 2,label='Kalman', color = 'b')
    ax2.text( max(x_axis[::index_plot]), max((v_error_mag[::index_plot])*100000.) - ( max((v_error_mag[::index_plot])*100000.) - min((v_error_mag[::index_plot])*100000.) )/7, "$\mathrm{RMSE_{kalm}}$ = " + '{:.2f}'.format(v_rmse*100000) + " cm/s", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right')
    ax2.text( max(x_axis[::index_plot]), max((v_error_mag[::index_plot])*100000.) - 2* ( max((v_error_mag[::index_plot])*100000.) - min((v_error_mag[::index_plot])*100000.) )/7, "$\mathrm{RMSE_{meas}}$ = " + '{:.2f}'.format(v_rmse_meas*100000) + " cm/s", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right')
    ax2.text( max(x_axis[::index_plot]), max((v_error_mag[::index_plot])*100000.) - 3* ( max((v_error_mag[::index_plot])*100000.) - min((v_error_mag[::index_plot])*100000.) )/7, "$\mathrm{RMSE_{prop}}$ = " + '{:.2f}'.format(v_spock_rmse*100000) + " cm/s", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right')
    ax2.margins(0,0)
    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_error_mag.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

    # ax2.text( max(x_axis[::index_plot]),  max((v_error_mag)*100000.) - ( max((v_error_mag)*100000.) - min((v_error_mag)*100000.) )/7, "$\mathrm{RMSE_{kalm}}$ = " + '{:.2f}'.format(v_rmse*100000) + " cm/s", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')
    # ax2.text( max(x_axis[::index_plot]),  max((v_error_mag)*100000.) - 2 * ( max((v_error_mag)*100000.) - min((v_error_mag)*100000.) )/7, "$\mathrm{RMSE_{meas}}$ = " + '{:.2f}'.format(v_rmse_meas*100000) + " cm/s", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')
    # ax2.text( max(x_axis[::index_plot]),  max((v_error_mag)*100000.) - 3* ( max((v_error_mag)*100000.) - min((v_error_mag)*100000.) )/7, "$\mathrm{RMSE_{prop}}$ = " + '{:.2f}'.format(v_spock_rmse*100000) + " cm/s", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')


#         ## Z
#     y_label = 'Error (m/s$^2$)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('$|\mathrm{a_{drag, kalm} - a_{drag, true}}|$', weight = 'normal', fontsize = fontsize_plot)
#     ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     # ax3.xaxis.set_ticklabels("")
#     ax3.plot(x_axis[::index_plot], (a_error_mag[::index_plot]*1000.), linewidth = 2,label='Kalman', color = 'b')
#     #print a_error_mag
#     ax3.text( max(x_axis[::index_plot]),  max((a_error_mag*1000.)) - ( max((a_error_mag*1000.)) - min((a_error_mag*1000.)) )/7, "$\mathrm{NRMSE_{kalm}}$ = " + '{:.2f}'.format(a_nrmse) , fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')
#     ax3.text( max(x_axis[::index_plot]),  max((a_error_mag*1000.)) - 2*( max((a_error_mag*1000.)) - min((a_error_mag*1000.)) )/7, "$\mathrm{NRMSE_{prop}}$ = " + '{:.2f}'.format(a_spock_nrmse), fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')
# #    ax3.text( max(x_axis[::index_plot]),  max((a_error_mag*1000.)) - 3 * ( max((a_error_mag*1000.)) - min((a_error_mag*1000.)) )/15, "NRMSE = " + '{:.3f}'.format(a_nrmse), fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')
#     ax3.text( max(x_axis[::index_plot]),  max((a_error_mag*1000.)) - 3 * ( max((a_error_mag*1000.)) - min((a_error_mag*1000.)) )/7, "$\overline{|\mathrm{a_{drag, true}}|}$ = " + '{:.1e}'.format(np.mean(a_comp_mag)*1000.) + " m/s$^2$", fontsize = fontsize_plot, weight = 'normal', color = 'k', horizontalalignment = 'right', verticalalignment = 'top')
#     ax3.margins(0,0)                                                                                                                  
#    ax3.set_ylim([min(a_error_mag*1000.[::index_plot]), max(a_error_mag*1000.[::index_plot])])
    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_error_mag.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    
#     # Plot error acceleration components VS time
#     fig_title = ''
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

#     fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(4, 1)
#     gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.8, wspace = 0.35)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
#     ## X
#     y_label = 'Error (m/s$^2$)'
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('$\mathrm{a_{x,kalm} - a_{x,true}}$', weight = 'normal', fontsize = fontsize_plot)
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax1.xaxis.set_ticklabels("")
#     ax1.plot(x_axis[::index_plot], (a_error[::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax1.margins(0,0)
#     ax1.annotate("RMSE = " + '{:.1e}'.format(ax_rmse*1000) + " m/s$^2$", (0,-0.35), (0,0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal')
#     ax1.annotate("NRMSE = " + '{:.1f}'.format(ax_nrmse), (0.5,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal', ha='center')
#     ax1.annotate("$\overline{|\mathrm{a_{x,true}}|}$  = " + '{:.1e}'.format(np.mean(np.abs(a_comp[::index_plot,0]))*1000.) + " m/s$^2$", (1,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal', ha = 'right')

#         ## X
#     y_label = 'Error (m/s$^2$)'
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('$\mathrm{a_{y,kalm} - a_{y,true}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
# #    ax2.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax2.xaxis.set_ticklabels("")
#     ax2.plot(x_axis[::index_plot], (a_error[::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax2.margins(0,0)
#     ax2.annotate("RMSE = " + '{:.1e}'.format(ay_rmse*1000) + " m/s$^2$", (0,-0.35), (0,0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal')
#     ax2.annotate("NRMSE = " + '{:.1f}'.format(ay_nrmse), (0.5,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal', ha='center')
#     ax2.annotate("$\overline{|\mathrm{a_{y,true}}|}$  = " + '{:.1e}'.format(np.mean(np.abs(a_comp[::index_plot,1]))*1000.) + " m/s$^2$", (1,-0.35), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal', ha = 'right')

#         ## X
#     y_label = 'Error (m/s$^2$)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('$\mathrm{a_{z,kalm} - a_{z,true}}$', weight = 'normal', fontsize = fontsize_plot)
# #    ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot, labelpad = 30)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     #ax3.xaxis.set_ticklabels("")
#     ax3.plot(x_axis[::index_plot], (a_error[::index_plot,2])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax3.margins(0,0)
#     # ax3.tick_params(direction='out', pad=15)
#     ax3.annotate("RMSE = " + '{:.1e}'.format(az_rmse*1000) + " m/s$^2$", (0,-0.60), (0,0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal')
#     ax3.annotate("NRMSE = " + '{:.1f}'.format(az_nrmse), (0.5,-0.60), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal', ha='center')
#     ax3.annotate("$\overline{|\mathrm{a_{z,true}}|}$  = " + '{:.1e}'.format(np.mean(np.abs(a_comp[::index_plot,2]))*1000.) + " m/s$^2$", (1,-0.60), (0, 0), xycoords='axes fraction', textcoords='offset points', va='bottom', fontsize = fontsize_plot, weight = 'normal', ha = 'right')

#     fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_error_acc_components.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


# #        LVLV ACCELERATION COMPONENTS OF KALMAN
#     # Plot error VS time
#     fig_title = ''
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

#     fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(3, 1)
#     gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
#     ## X
#     y_label = '$\mathrm{a_{x,kalm}}$ (m/s$^2$)'
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('$\mathrm{a_{x,kalm}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax1.xaxis.set_ticklabels("")
#     ax1.plot(x_axis[1::index_plot], (a_kalm[1::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
#     #legend =o ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
#     ax1.margins(0,0)
#     ax1.set_ylim([min(a_kalm[2::index_plot,0])*1000., max(a_kalm[2::index_plot,0])*1000.])

#     ## Y
#     y_label = '$\mathrm{a_{y,kalm}}$ (m/s$^2$)'
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('$\mathrm{a_{y,kalm}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax2.xaxis.set_ticklabels("")
#     ax2.plot(x_axis[1::index_plot], (a_kalm[1::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax2.margins(0,0)
#     ax2.set_ylim([min(a_kalm[2::index_plot,1])*1000., max(a_kalm[2::index_plot,1])*1000.])

#         ## Z
#     y_label = '$\mathrm{a_{z,kalm}}$ (m/s$^2$)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('$\mathrm{a_{z,kalm}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     # ax3.xaxis.set_ticklabels("")
#     ax3.plot(x_axis[1::index_plot], (a_kalm[1::index_plot,2])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax3.margins(0,0)
#     ax3.set_ylim([min(a_kalm[2::index_plot,2])*1000., max(a_kalm[2::index_plot,2])*1000.])

#     fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration_kalman.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

# #    LVLV ACCELERATION COMPONENTS OF COMP
#     # Plot error VS time
#     fig_title = ''
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

#     fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(3, 1)
#     gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
#     ## X
#     y_label = '$\mathrm{a_{x,true}}$ (m/s$^2$)'
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('$\mathrm{a_{x,true}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax1.xaxis.set_ticklabels("")
#     ax1.plot(x_axis[1::index_plot], (a_comp[1::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
#     #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
#     ax1.margins(0,0)

#     ## Y
#     y_label = '$\mathrm{a_{y,true}}$ (m/s$^2$)'
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('$\mathrm{a_{y,true}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax2.xaxis.set_ticklabels("")
#     ax2.plot(x_axis[1::index_plot], (a_comp[1::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax2.margins(0,0)


#         ## Z
#     y_label = '$\mathrm{a_{z,true}}$ (m/s$^2$)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('$\mathrm{a_{z,true}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     # ax3.xaxis.set_ticklabels("")
#     ax3.plot(x_axis[1::index_plot], (a_comp[1::index_plot,  2])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax3.margins(0,0)

#     fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration_' + output_file_name_list_comp[isc].replace(".txt", ".pdf")
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


#     itime = 0
#     mag_delta_r_temp = r_kalm[itime,:] - r_comp[itime,:]
#     mag_delta_v_temp = v_kalm[itime,:] - v_comp[itime,:]
#     print "Magnitude of delta r at time " + str(itime) + ": ", np.linalg.norm(mag_delta_r_temp) * 1000, "m"
#     print "Magnitude of delta v at time " + str(itime) + ": ", np.linalg.norm(mag_delta_v_temp)*1000., "m/s"


# #        LVLV ACCELERATION COMPONENTS OF SpOCK
#     # Plot error VS time
#     fig_title = ''
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

#     fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(3, 1)
#     gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
#     ## X
#     y_label = '$\mathrm{a_{x,spock}}$ (m/s$^2$)'
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('$\mathrm{a_{x,spock}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax1.xaxis.set_ticklabels("")
#     ax1.plot(x_axis[1::index_plot], (a_spock[1::index_plot,0])*1000., linewidth = 2,label='Spock', color = 'b')
#     #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
#     ax1.margins(0,0)

#     ## Y
#     y_label = '$\mathrm{a_{y,spock}}$ (m/s$^2$)'
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('$\mathrm{a_{y,spock}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax2.xaxis.set_ticklabels("")
#     ax2.plot(x_axis[1::index_plot], (a_spock[1::index_plot,1])*1000., linewidth = 2,label='Spock', color = 'b')
#     ax2.margins(0,0)

#         ## Z
#     y_label = '$\mathrm{a_{z,spock}}$ (m/s$^2$)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('$\mathrm{a_{z,spock}}$', weight = 'normal', fontsize = fontsize_plot)
#     ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     # ax3.xaxis.set_ticklabels("")
#     ax3.plot(x_axis[1::index_plot], (a_spock[1::index_plot,2])*1000., linewidth = 2,label='Spock', color = 'b')
#     ax3.margins(0,0)

#     fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration_propagation_only.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  




# #        LVLV ACCELERATION COMPONENTS OF KALMAN AND TRUE
#     # Plot error VS time
#     fig_title = ''
#     fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
#     fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
#     gs = gridspec.GridSpec(3, 1)
#     gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
#     plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
#     ## X
#     y_label = '$\mathrm{a_{x,kalm}}$ (m/s$^2$)'
#     ax1 = fig.add_subplot(gs[0, 0])
#     ax1.set_title('$\mathrm{a_{x,kalm}}$ (b), $\mathrm{a_{x,true}}$ (r), $\mathrm{a_{x,spock}}$ (k)', weight = 'normal', fontsize = fontsize_plot)
#     ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
#     ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     #ax1.xaxis.set_ticklabels("")
#     ax1.plot(x_axis[1::index_plot], (a_kalm[1::index_plot,0]+a_spock[1::index_plot,0])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax1.plot(x_axis[1::index_plot], (a_comp[1::index_plot,0])*1000., linewidth = 2,label='True', color = 'r')
#     ax1.plot(x_axis[1::index_plot], (a_spock[1::index_plot,0])*1000., linewidth = 2,label='Propagation', color = 'k')
#     #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
# #    ax1.set_ylim([min(a_kalm[::index_plot,0][2*60:]+a_spock[::index_plot,0][2*60:])*1000, max(a_kalm[::index_plot,0][2*60:]+a_spock[::index_plot,0][2*60:])*1000])
#     ax1.set_ylim([min([min(a_spock[2*3600::index_plot,0]*1000.), min((a_comp[2*3600::index_plot,0])*1000.), min((a_kalm[2*3600::index_plot,0]+a_spock[2*3600::index_plot,0])*1000.)]), max([max(a_spock[2*3600::index_plot,0]*1000.), max((a_comp[2*3600::index_plot,0])*1000.), max((a_kalm[2*3600::index_plot,0]+a_spock[2*3600::index_plot,0])*1000.)])])
#     ax1.margins(0,0)

#     # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration_kalman_and_true_and_propagation_only.pdf'
#     # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  

#     ## Y
#     y_label = '$\mathrm{a_{y,kalm}}$ (m/s$^2$)'
#     ax2 = fig.add_subplot(gs[1, 0])
#     ax2.set_title('$\mathrm{a_{y,kalm}}$ (b), $\mathrm{a_{y,true}}$ (r), $\mathrm{a_{y,spock}}$ (k)', weight = 'normal', fontsize = fontsize_plot)
#     ax2.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax2.spines.itervalues()] # change the width of the frame of the figure
#     ax2.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     ax2.xaxis.set_ticklabels("")
#     ax2.plot(x_axis[1::index_plot], (a_kalm[1::index_plot,1]+a_spock[1::index_plot,1])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax2.plot(x_axis[1::index_plot], (a_comp[1::index_plot,1])*1000., linewidth = 2,label='True', color = 'r')
#     ax2.plot(x_axis[1::index_plot], (a_spock[1::index_plot,1])*1000., linewidth = 2,label='Propagation', color = 'k')
#     ax2.set_ylim([min([min(a_spock[2*3600::index_plot,1]*1000.), min((a_comp[2*3600::index_plot,1])*1000.), min((a_kalm[2*3600::index_plot,1]+a_spock[2*3600::index_plot,1])*1000.)]), max([max(a_spock[2*3600::index_plot,1]*1000.), max((a_comp[2*3600::index_plot,1])*1000.), max((a_kalm[2*3600::index_plot,1]+a_spock[2*3600::index_plot,1])*1000.)])])
#     ax2.margins(0,0)


#         ## Z
#     y_label = '$\mathrm{a_{z}}$ (m/s$^2$)'
#     ax3 = fig.add_subplot(gs[2, 0])
#     ax3.set_title('$\mathrm{a_{z,kalm}}$ (b), $\mathrm{a_{z,true}}$ (r), $\mathrm{a_{z,spock}}$ (k)', weight = 'normal', fontsize = fontsize_plot)
#     ax3.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
#     ax3.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
#     [i.set_linewidth(2) for i in ax3.spines.itervalues()] # change the width of the frame of the figure
#     ax3.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
#     # ax3.xaxis.set_ticklabels("")
#     ax3.plot(x_axis[1::index_plot], (a_kalm[1::index_plot,2]+a_spock[1::index_plot,2])*1000., linewidth = 2,label='Kalman', color = 'b')
#     ax3.plot(x_axis[1::index_plot], (a_comp[1::index_plot,2])*1000., linewidth = 2,label='True', color = 'r')
#     ax3.plot(x_axis[1::index_plot], (a_spock[1::index_plot,2])*1000., linewidth = 2,label='Propagation', color = 'k')
#     ax3.set_ylim([min([min(a_spock[2*3600::index_plot,2]*1000.), min((a_comp[2*3600::index_plot,2])*1000.), min((a_kalm[2*3600::index_plot,2]+a_spock[2*3600::index_plot,2])*1000.)]), \
#                   max([max(a_spock[2*3600::index_plot,2]*1000.), max((a_comp[2*3600::index_plot,2])*1000.), max((a_kalm[2*3600::index_plot,2]+a_spock[2*3600::index_plot,2])*1000.)])])
#     ax3.margins(0,0)

#     fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_acceleration_kalman_and_true_and_propagation.pdf'
#     fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


#     MAGNITUDE OF ACCELERATION KALMAN TRUE PROPAGATION
    # Plot error VS time
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    
    ## X
    y_label = '$\mathrm{a}$ (m/s$^2$)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('$\mathrm{a_{kalm}}$ (b), $\mathrm{a_{true}}$ (r), $\mathrm{a_{spock}}$ (k)', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    #ax1.plot(x_axis[1::index_plot], (a_kalm_mag[1::index_plot])*1000., linewidth = 2,label='Kalman', color = 'b')
    ax1.plot(x_axis[1::index_plot], (a_kalm_total_drag_from_c_mag[1::index_plot])*1000., linewidth = 2,label='Kalman', color = 'b')
    ax1.plot(x_axis[1::index_plot], (a_comp_mag[1::index_plot])*1000., linewidth = 2,label='True', color = 'r')
    ax1.plot(x_axis[1::index_plot], (a_spock_mag[1::index_plot])*1000., linewidth = 2,label='Propagation', color = 'k')
    
    #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
#    ax1.set_ylim([min(a_kalm[::index_plot,0][2*60:]+a_spock[::index_plot,0][2*60:])*1000, max(a_kalm[::index_plot,0][2*60:]+a_spock[::index_plot,0][2*60:])*1000])
#     ax1.set_ylim([min([min(a_spock_mag[2*3600::index_plot]*1000.), min((a_comp_mag[2*3600::index_plot])*1000.), min(a_kalm_mag[2*3600::index_plot]*1000.)]),\
#                       max([max(a_spock_mag[2*3600::index_plot]*1000.), max((a_comp_mag[2*3600::index_plot])*1000.), max(a_kalm_mag[2*3600::index_plot]*1000.)])])

    ax1.set_ylim([min((a_comp_mag[1::index_plot])*1000.) * 0.9, max((a_comp_mag[1::index_plot])*1000.) * 1.1])
    ax1.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_magnitude_acceleration_kalman_and_true_and_propagation_only.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    ## TAU OF GAUSS MARKOV POCESS
    fig_title = ''
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(3, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    plt.rc('font', weight='normal') ## make the labels of the ticks in bold

    y_label = 'tau (s)'
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_title('tau in Gauss-Markov process', weight = 'normal', fontsize = fontsize_plot)
    ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    #ax1.xaxis.set_ticklabels("")
    ax1.plot(x_axis[1::index_plot], (tau[1::index_plot]), linewidth = 2,label='Kalman', color = 'b')
    ax1.margins(0,0)

    fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_tau.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    #     ## Cd
    # fig_title = ''
    # fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    # fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold
    # gs = gridspec.GridSpec(3, 1)
    # gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
    # plt.rc('font', weight='normal') ## make the labels of the ticks in bold

    # y_label = 'Cd'
    # ax1 = fig.add_subplot(gs[0, 0])
    # ax1.set_title('Cd', weight = 'normal', fontsize = fontsize_plot)
    # ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
    # ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
    # [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
    # ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
    # #ax1.xaxis.set_ticklabels("")
    # ax1.plot(x_axis[1::index_plot], (cd_kalm[1::index_plot]), linewidth = 2,label='Kalman', color = 'b')
    # ax1.margins(0,0)

    # fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_cd.pdf'
    # fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  



    if plot_pred == 1: # this is the plot same as Kalman but using only the propagation, not the Kalman filter (equivalent to K = 0)
        var_to_read = ["position","velocity", "acceleration_lvlh", "density"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        r_pred = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
        v_pred = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]
        a_pred = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration_lvlh')]
        rho_pred = var_out[find_in_read_input_order_variables(var_out_order, 'density')] * 10**9 # * 10**9 to convert from kg/m^3 to kg/km^3

        # Plot rho VS time
        fig_title = ''
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

        fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
        plt.rc('font', weight='normal') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
        plt.rc('font', weight='normal') ## make the labels of the ticks in bold

        ## X
        y_label = 'Mass density (kg/m$^3$)'
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Mass density as a function of time', weight = 'normal', fontsize = fontsize_plot)
        ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
        ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
        #ax1.xaxis.set_ticklabels("")
        ax1.scatter(x_axis[::index_plot], rho_pred[::index_plot], linewidth = 2,label='Prediction', color = 'b')
        ax1.scatter(x_axis[::index_plot], rho_comp[::index_plot], linewidth = 2,label='True', color = 'k')
        #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
        ax1.set_ylim([min([min(rho_pred), min(rho_comp)]), max([max(rho_pred), max(rho_comp)])])
        ax1.margins(0,0)

        fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_rho_pred.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


        # Plot ad VS time
        fig_title = ''
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

        fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
        plt.rc('font', weight='normal') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(4, 1)
        gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.3, wspace = 0.35)
        plt.rc('font', weight='normal') ## make the labels of the ticks in bold

        ## X
        y_label = '$|\mathrm{a_{LVLH}}|$ (m/s$^2$)'
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.set_title('Magnitude of LVLH acceleration as a function of time', weight = 'normal', fontsize = fontsize_plot)
        ax1.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
        ax1.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)
        [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure
        ax1.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7)
        #ax1.xaxis.set_ticklabels("")
        ax1.scatter(x_axis[::index_plot], a_pred[::index_plot,0]*1000, linewidth = 2,label='Predictio', color = 'b')
        ax1.scatter(x_axis[::index_plot], a_comp[::index_plot,0]*1000, linewidth = 2,label='True', color = 'k')
        #legend = ax1.legend(loc = 'upper right',  fontsize = fontsize_plot)
        ax1.set_ylim([min([min(a_pred[::index_plot,0]*1000), min(a_comp[::index_plot,0]*1000)]), max([max(a_pred[::index_plot,0]*1000), max(a_comp[::index_plot,0]*1000)])])
        ax1.margins(0,0)

        fig_save_name = output_file_name_list[isc].replace(".txt", "") + '_ad_pred.pdf'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


