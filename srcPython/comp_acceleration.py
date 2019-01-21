# This script was written to compare the total accelerations between several runs


import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import numpy as np
from read_input_file import *
from read_output_file import *
from find_in_read_input_order_variables import *

filename_ref = "order4.txt"
filename = ["mapAccu.txt", "order10.txt", "order20.txt", "order50.txt"]
nfile = len(filename)

# Read acceleration for all files (ref and to compare)
var_to_read = ["acceleration"]
for ifile in range(nfile+1):
    # Get variable from input files 
    if ifile == 0:
        filename_now = filename_ref
    else:
        filename_now = filename[ifile - 1]
    var_in, var_in_order = read_input_file(filename_now)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
    nb_steps_now = var_in[find_in_read_input_order_variables(var_in_order, 'nb_steps')]; 
    # Read the acceleration from output files
    if ifile == 0:
        nb_steps = nb_steps_now
        acc = np.zeros([nfile, nb_steps, 3])
        delta_acc = np.zeros([nfile, nb_steps]) # magnitude of difference in acceleration between ref and run
        cumul_delta_acc = np.zeros([nfile]) # cumulative of magnitude of difference in acceleration
        var_out, var_out_order = read_output_file( output_file_path_list[0] + output_file_name_list[0], var_to_read )
        date = var_out[find_in_read_input_order_variables(var_out_order, 'date')]
        nb_seconds_since_start = var_out[find_in_read_input_order_variables(var_out_order, 'nb_seconds_since_start')]
        acc_ref = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration')]
    else:
        if nb_steps_now != nb_steps:
            sys.exit("!***!\nThe number of steps must be the same for all runs.\n!***!")
            var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
        var_out, var_out_order = read_output_file( output_file_path_list[0] + output_file_name_list[0], var_to_read )
        acc[ifile-1, :, :] = var_out[find_in_read_input_order_variables(var_out_order, 'acceleration')]
        delta_acc[ifile-1, :] = np.linalg.norm(acc[ifile-1, :, :] - acc_ref, axis = 1)
        cumul_delta_acc[ifile-1] = np.sum(delta_acc[ifile-1, :])

    
    
    

        
