import subprocess
import fileinput
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import numpy as np
import sys
from read_input_file import *
from read_output_file import *
from get_prop_dir import *

plt.ion()

max_order = 40
step_order = 2
order_array = np.arange(2, max_order + 1, step_order)[::-1] # start with the highest order
nb_simu = len(order_array)
runtime = np.zeros([nb_simu])
r_diff_mag_max = np.zeros([nb_simu])
#file_out_error_runtime = open('./paper3_stk_results/runtime_error_order_gravity.txt', "w")
#print >> file_out_error_runtime, "Order Error(km) Runtime(s)"
for iorder in range(nb_simu):
            ## Write input file
    os.chdir("../run_paper3/")
    input_filename_propagator = 'alt_300_gravity_order_' + str(order_array[iorder]) 
    # if iorder > 0:
    #     os.system("cp input/main_input/alt_300_gravity_order_" + str(order_array[0]) + ".txt " + "input/main_input/" + input_filename_propagator + ".txt")
    # line_number = 0
    # for line in fileinput.input("input/main_input/" + input_filename_propagator + ".txt", inplace = True):
    #     line_number = line_number + 1
    #     if ( ( line_number != 21 ) & ( line_number != 28 ) ):
    #         print "%s" % (line),
    #     elif ( line_number == 21 ):
    #         string_here = str(order_array[iorder])
    #         print "%s\n" % string_here,
    #     elif ( line_number == 28 ):
    #         string_here = input_filename_propagator 
    #         print "%s\n" % string_here,
    # #Run propagator
    # input_filename_propagator_with_txt = input_filename_propagator + '.txt'
    # proc = subprocess.Popen(['/usr/bin/time','-p','/usr/local/bin/mpirun', '-np' ,'1' ,'spock',  input_filename_propagator_with_txt], stderr=subprocess.PIPE)
    # output=proc.stderr.read()
    # runtime[iorder] = output.split('\n')[0].split()[1]
    # Read position and velocity from our propagator
    input_filename_propagator_complete = get_prop_dir(1) + "run_paper3/input/main_input/"  + input_filename_propagator + '.txt'
    input_variables, order_input_variables = read_input_file(input_filename_propagator_complete)
    output_propagator_path = input_variables[6][0]; output_propgator_file = input_variables[7][0]
    to_output = ["position"]
    out, out_var = read_output_file(output_propagator_path + output_propgator_file, to_output)
    if iorder == 0: # highest order is reference
        r_eci_j2 = out[1]
    else:
        r_eci = out[1]
        r_diff = r_eci - r_eci_j2
        r_diff_mag = np.zeros([len(out[0])])
        for i in range(len(out[0])):
            r_diff_mag[i] = np.linalg.norm(r_diff[i]) # in km
        r_diff_mag_max[iorder] = np.max(r_diff_mag)
    # print r_diff_mag_max
    # print runtime
    # print >> file_out_error_runtime, order_array[iorder], r_diff_mag_max[iorder], runtime[iorder]
    #os.system("/usr/local/bin/mpirun -np 1 spock " + input_filename_propagator + ".txt")
#file_out_error_runtime.close()
os.chdir("../srcPython")

width_fig = 15.
height_fig = width_fig * 3 /4
fontsize_plot = 20 # 9
color_arr = ['k','b','r','g','m', 'y']

# ERROR
fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.97, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
#ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

ax1.semilogy(order_array[::-1], r_diff_mag_max[::-1], linewidth = 2, color = 'k')
ax1.set_xlim([0, 40])
ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

ax1.set_title('Error as a function of the gravity model order', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
ax1.legend( fontsize = fontsize_plot, loc = 4)
ax1.set_xlabel('Gravity model order', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Error (km)', fontsize = fontsize_plot, weight = 'bold')

fig_save_name = './paper3_stk_results/error_order_gravity_log.eps'
fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")


# # RUNTIME
# fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
# gs = gridspec.GridSpec(1, 1)
# gs.update(left=0.1, right=0.97, top = 0.90,bottom = 0.06)
# fig.suptitle('', fontsize = 22)
# plt.rc('font', weight='bold') ## make the labels of the ticks in bold
# ax1 = fig.add_subplot(gs[0, 0])
# #ax1 = fig.add_subplot(111)
# ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
# ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
# #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
# [i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure

# ax1.plot(order_array[::-1], runtime[::-1], linewidth = 2, color = 'k')
# ax1.set_xlim([0, 40])
# ax1.set_ylim([0, ax1.get_ylim()[1]])
# ax1.margins(0,0) # autoscale both axes(fist value is for the x axis, second value for the y axis)

# ax1.set_title('Computational time of the gravity model order', weight = 'bold', fontsize = fontsize_plot , y = 1.005) 
# ax1.legend( fontsize = fontsize_plot, loc = 4)
# ax1.set_xlabel('Gravity model order', fontsize = fontsize_plot, weight = 'bold')
# ax1.set_ylabel('Computational time (s)', fontsize = fontsize_plot, weight = 'bold')

# fig_save_name = './paper3_stk_results/runtime_order_gravity.eps'
# fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
# os.system("rsync -av " + fig_save_name + " srbwks2014-0008.engin.umich.edu:./")
