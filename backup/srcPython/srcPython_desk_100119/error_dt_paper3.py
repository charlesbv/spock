import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import numpy as np
import sys
from read_input_file import *
from read_output_file import *
from get_prop_dir import *

plt.ion()

##########IMPORTANT: FOR THE 600 S TIME STEP, GO TO THE OUTPUT FILE AND REMOVE ALL LINES EXCEPT THE FIRST AND LAST ONE
# # File with the results of the difference in positions and velocities 
# file_results = open("./paper3_stk_results/results_error_dt.txt", "w+")
# print >> file_results, "name_run max_diff_pos max_diff_vel (distances are in kilometers)\n"

# Runs
##########IMPORTANT: FOR THE 600 S TIME STEP, GO TO THE OUTPUT FILE AND REMOVE ALL LINES EXCEPT THE FIRST AND LAST ONE
input_filename_array = ['runge_kutta_dt_1s.txt', 'runge_kutta_dt_6s.txt', 'runge_kutta_dt_10s.txt', 'runge_kutta_dt_60s.txt', 'runge_kutta_dt_100s.txt', 'runge_kutta_dt_600s.txt']
dt_arr = [1,5,10,30,60,120,300,600]

# # Runs
# for idt in range(len(dt_arr)):
#     ## Write input file
#     os.chdir("../run_paper3/")
#     input_filename_propagator = 'runge_kutta_dt_' + str(dt_arr[idt] ) + '_s'
#     os.system("cp input/main_input/runge_kutta_dt_1s.txt " + "input/main_input/" + input_filename_propagator + ".txt")
#     line_number = 0
#     for line in fileinput.input("input/main_input/" + input_filename_propagator + ".txt", inplace = True):
#         line_number = line_number + 1
#         if ( ( line_number != 4 ) & ( line_number != 9 ) ):
#             print "%s" % (line),
#         elif ( line_number == 4 ) :
#             string_here = str(dt_arr[idt]) + ", 86400"
#             print "%s\n" % string_here,
#         elif ( line_number == 9 ):
#             string_here = input_filename_propagator + ".txt"
#             print "%s\n" % string_here,
#     ## Run propagator
#     os.system("/usr/local/bin/mpirun -np 1 run_moat " + input_filename_propagator + ".txt")
# os.chdir("../srcPython")

max_r_diff_mag = np.zeros([len(dt_arr)-1]) # -1 because we don't care about teh difference between dt  = 1s and dt = 1s
max_v_diff_mag = np.zeros([len(dt_arr)-1])
##### BEGINNING OF LOOP
for irun in range(len(dt_arr)):
    # Read position and velocity from our propagator
    input_filename_prop = input_filename_propagator = 'runge_kutta_dt_' + str(dt_arr[irun] ) + '_s' + '.txt'
    input_filename_prop_complete = get_prop_dir(1) + "run_paper3/input/main_input/"  + input_filename_prop
    input_variables, order_input_variables = read_input_file(input_filename_prop_complete)
    output_propagator_path = input_variables[6][0]; output_propgator_file = input_variables[7][0]
    to_output = ["position", "velocity"]
    out, out_var = read_output_file(output_propagator_path + output_propgator_file, to_output)
    if irun == 0: # dt = 1s
        r_eci_dt_1s = out[1]; v_eci_dt_1s = out[2]
    else:
        # Difference between positions and velocities
        r_eci_prop = out[1]; v_eci_prop = out[2]
        r_diff = r_eci_dt_1s - r_eci_prop
        v_diff = v_eci_dt_1s - v_eci_prop
##########IMPORTANT: FOR THE 600 S TIME STEP, GO TO THE OUTPUT FILE AND REMOVE ALL LINES EXCEPT THE FIRST AND LAST ONE
        r_diff_mag = np.zeros([len(out[0])])
        v_diff_mag = np.zeros([len(out[0])])
        for i in range(len(out[0])):
            r_diff_mag[i] = np.linalg.norm(r_diff[i]) # in km
            v_diff_mag[i] = np.linalg.norm(v_diff[i]) # in km/s

        max_r_diff_mag[irun-1] = np.max(r_diff_mag)
        max_v_diff_mag[irun-1] = np.max(v_diff_mag)
        # Write results in file
#         print >> file_results, input_filename_prop.split('.')[0] + ' ' + '{0:.3f}'.format(np.max(r_diff_mag)) + ' ' + '{0:.3f}'.format(np.max(v_diff_mag))

#     ##### END OF LOOP

# # Close file with results
# file_results.close()



# PLOT
width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.05, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
ax1.tick_params(axis='x', which='minor',  width = 2, pad = 6, labelsize=fontsize_plot, size = 7)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

ax1.loglog(dt_arr[1:],max_r_diff_mag, 'k', marker = 'o', linewidth = 2, linestyle = '-', markersize = 10)

ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

ax1.set_title('Error as a function of the integration step size', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
ax1.legend(ncol = 2, fontsize = fontsize_plot)
ax1.set_xlabel('Integration step size (s)', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Error (km)', fontsize = fontsize_plot, weight = 'bold')

fig_save_name = './paper3_stk_results/error_dt.eps'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")
