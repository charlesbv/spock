# this script compares the ECEF position between SpOCK and the netcdf data, as well as the VLH position of SpOCK with respect to netcdf
# it does so over multiple SpOCK's 1-day simulations between date_start_analysis and date_stop_analysis for the FMs included in fm_analysis.
# Each simlation uses 7-day old CYGNSS TLE so that the position error reflects errors for propagation over 7 to 8 days
# ASSUMPTIONS:
# - see section PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
# - the netcdf files must previously be downloaded and located in netcdf_dir_root

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
date_start_analysis = '2017-03-18'
date_stop_analysis = '2019-01-03'
fm_analysis = [1,2,3,4,5,6,7,8]
netcdf_dir_root = '/Users/cbv/cygnss/netcdf' #'/Volumes/Seagate_Expansion_Drive/Backups.backupdb/srblap2017-0085/2019-04-14-150952/Macintosh HD/Users/cbv/cygnss/netcdf'#'/Users/cbv/cygnss/netcdf'
sgp4 = 1
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

from datetime import datetime, timedelta
import os
import ipdb
import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
from read_input_file import *
from read_output_file import *
from find_in_read_input_order_variables import *
from cygnss_read_netcdf_and_convert_to_eci_and_oe import *
from ecef_to_lvlh import *
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
if netcdf_dir_root[-1] != '/':
    netcdf_dir_root = netcdf_dir_root + '/'

nfm = len(fm_analysis)

date_start_analysis_date = datetime.strptime(date_start_analysis, "%Y-%m-%d")
date_stop_analysis_date = datetime.strptime(date_stop_analysis, "%Y-%m-%d")
nday = (date_stop_analysis_date - date_start_analysis_date).days + 1
day_jump = 7
date_array = np.array([date_start_analysis_date + timedelta(days=i) for i in np.arange(0,nday,day_jump)])
nday = len(date_array)


along_rms = np.zeros([nday, nfm])-1; cross_rms = np.zeros([nday, nfm])-1; radial_rms = np.zeros([nday, nfm])-1; mag_rms = np.zeros([nday, nfm])-1; npoint = np.zeros([nday, nfm])-1;
load_spice = 1
for iday in range(nday):
    date_start_date = date_array[iday]
    date_stop_date = date_start_date + timedelta(days = 1)
    date_start = datetime.strftime(date_start_date, "%Y-%m-%dT%H:%M:%S")
    date_stop = datetime.strftime(date_stop_date, "%Y-%m-%dT%H:%M:%S")
    print iday, nday-1, date_start
    # run SpOCK
    if sgp4 == 1:
        os.system("spock_cygnss_spec_parallel_using_one_week_old_tle.py " + date_start + " " + date_stop + " 1")
    else:
        os.system("spock_cygnss_spec_parallel_using_one_week_old_tle.py " + date_start + " " + date_stop + " 0")
        
    # Read SpOCK's main input file
    if sgp4 == 1:
        input_filename = 'spock_spec_start_' + date_start.replace(":", "_") + '_end_' + date_stop.replace(":", "_") + '_sgp4_using_one_week_old_tle.txt'
    else:
        input_filename = 'spock_spec_start_' + date_start.replace(":", "_") + '_end_' + date_stop.replace(":", "_") + '_classic_using_one_week_old_tle.txt'
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
    dt = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
    nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];

    # for each FM, read SpOCK r/v
    fm_to_sc = [5, 4, 2, 1, 8, 6, 7, 3]
    isc_count = -1
    for ifm in range(nfm):
        isc_count = isc_count + 1
        fm = fm_analysis[ifm]
        isc = fm_to_sc.index(fm)
        var_to_read = ["position_ecef"]
        var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        if isc_count == 0: 
            date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
            date_rounded_datetime_spock = var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime_round_sec')]
            nstep_spock = len(date)
        r_spock = np.zeros([nstep_spock, 3])
        r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position_ecef')]

        # for each fm, read netcdf r/v
        yy = date_start[:4] + '/'
        doy = datetime.strftime(date_start_date, "%Y-%j").split('-')[1] + '/'
        netcdf_dir = netcdf_dir_root + yy + doy.zfill(3)
        netcdf_filename_raw = [x for x in os.listdir(netcdf_dir) if x.endswith(".nc") and x.startswith("cyg0" + np.str(fm))]
        if len(netcdf_filename_raw) != 0:
            netcdf_filename_raw = netcdf_filename_raw[0]
            netcdf_filename = netcdf_dir + netcdf_filename_raw

            date_flight_rounded, lon_cyg, lat_cyg, lon_spec, lat_spec, fom, gps,\
            x_cyg, y_cyg, z_cyg, vx_cyg, vy_cyg, vz_cyg, date_rounded_datetime_netcdf,\
            r_eci_cyg, v_eci_cyg, sma, eccentricity, inclination, long_an, w, \
            f, period, phase_angle = cygnss_read_netcdf_and_convert_to_eci_and_oe(netcdf_filename, load_spice)
            load_spice = 0 # from now on load_spice = 0. We want load_spice to be 1 for the first call of cygnss_read_netcdf_and_convert_to_eci_and_oe only
            nstep_netcdf = len(date_flight_rounded)
            date_rounded_datetime_netcdf = np.array(date_rounded_datetime_netcdf)
            r_netcdf = np.array([x_cyg, y_cyg, z_cyg]).transpose()/1000.
            v_netcdf = np.array([vx_cyg, vy_cyg, vz_cyg]).transpose()/1000.
            # Compare the r between SpOCK and netcdf
            # netcdf output is every second but not necessarily SpOCK output so make sure the r is cmopared at the same date between SpOCK and netcdf
            inetcdf_prev = 0
            along_error = []; cross_error = []; radial_error = []; mag_error = []
            for ispock in range(nstep_spock):
                inetcdf = np.where(date_rounded_datetime_netcdf[inetcdf_prev:] == date_rounded_datetime_spock[ispock])[0]
                if len(inetcdf) != 0: # the date is present in SpOCK and in netcdf
                    inetcdf = inetcdf[0]  + inetcdf_prev
                    rdiff = r_spock[ispock, :] - r_netcdf[inetcdf, :]
                    rdiff_lvlh = ecef_to_lvlh(r_netcdf[inetcdf, :], v_netcdf[inetcdf, :], rdiff)
                    mag_error.append( np.linalg.norm(rdiff) )
                    along_error.append(rdiff_lvlh[0])
                    cross_error.append(rdiff_lvlh[1])
                    radial_error.append(rdiff_lvlh[2])

                    inetcdf_prev = inetcdf
            along_error = np.array(along_error); cross_error = np.array(cross_error); radial_error = np.array(radial_error); mag_error = np.array(mag_error)
            along_rms[iday, ifm] = np.sqrt(np.mean(along_error**2)); cross_rms[iday, ifm] = np.sqrt(np.mean(cross_error**2)); radial_rms[iday, ifm] = np.sqrt(np.mean(radial_error**2)); mag_rms[iday, ifm] = np.sqrt(np.mean(mag_error**2))
            npoint[iday, ifm] = len(along_error)
            print ifm, (int)(along_rms[iday, ifm]), (int)(cross_rms[iday, ifm]), (int)(radial_rms[iday, ifm]), (int)(mag_rms[iday, ifm]), (int)(npoint[iday, ifm])

## Parameters for the figure
height_fig = 11.  # the width is calculated as height_fig * 4/3.
fontsize_plot = 25
ratio_fig_size = 4./3
fig_title = 'FM 1-week along-track, cross-track and radial RMS'
y_label = 'RMS'
x_label = 'Sample'
fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'normal',)
plt.rc('font', weight='normal') ## make the labels of the ticks in bold
gs = gridspec.GridSpec(1, 1)
gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
ax = fig.add_subplot(gs[0, 0])

ax.set_ylabel(y_label, weight = 'normal', fontsize  = fontsize_plot)
ax.set_xlabel(x_label, weight = 'normal', fontsize  = fontsize_plot)

[i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
plt.rc('font', weight='normal') ## make the labels of the ticks in bold

ax.plot(np.max(along_rms, axis = 1), linewidth = 2, color = 'blue', label = 'Along-track')
ax.plot(np.max(cross_rms, axis = 1), linewidth = 2, color = 'red', label = 'Cross-track')
ax.plot(np.max(radial_rms, axis = 1), linewidth = 2, color = 'black', label  = 'Radial')
ax.plot(np.max(mag_rms, axis = 1), linewidth = 2, color = 'limegreen', label = 'Mag')
ax.margins(0,0)
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1, fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))
if sgp4 == 1:
    fig_save_name = 'spock_sgp4_vs_netcdf.pdf'
else:
    fig_save_name = 'spock_classic_vs_netcdf.pdf'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
