# this uses read_collision_file_new.py for a buch of runs

import sys
sys.path.append("/home1/cbussy/Code/spock/srcPython")

from read_input_file import *
from read_collision_file_new_with_dca import *
from matplotlib import pyplot as plt

# read file including all run main input filename
filename_with_runs = "runlist_FM07_nov01_30000ens.txt"
file_with_runs = open(filename_with_runs)
read_file_with_runs = file_with_runs.readlines()
run_list = []
for irun in range(len(read_file_with_runs)):
    run_list.append(read_file_with_runs[irun].split()[0])

nb_run = len(run_list)
pc = []
collision_file_missing_names = []
filename_collision_missing = 'missing_again_' + filename_with_runs
file_collision_missing = open(filename_collision_missing, "w") # this file will be sent ot pleiades to run all these failed simu again (with encounter_geo_collision_run.py on pleiades)
filename_results = 'out_' + filename_with_runs
file_results = open(filename_results, "w")
nb_run_more_than_one_ca = 0
for irun in range(nb_run):
    print irun, nb_run-1
    input_filename =  run_list[irun]
    var_in, var_in_order = read_input_file(input_filename)
    output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
    collision_filename = '/'.join(output_file_path_list[0].split('/')[:-2]) + "/" + '/'.join(output_file_path_list[0].split('/')[:-2]) + "_collision.txt"
    if os.path.isfile(collision_filename): 
        date, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll  = read_collision_file_new_with_dca(collision_filename)
        if len(tca) > 0: # otherwise no close approach
            pc.append(cpc_final[0])
            if len(tca)>1:
                nb_run_more_than_one_ca = nb_run_more_than_one_ca + 1
            print >> file_results, input_filename, cpc_final[0], tca[0]
        else: # no close approach
            print >> file_results, input_filename, 0, -1
    else:
        collision_file_missing_names.append(input_filename)
        print >> file_collision_missing, input_filename



pc = np.array(pc)
    
file_collision_missing.close()
file_results.close()

# fig,ax = plt.subplots()
# ax.plot(pc)
# fig.savefig("pc.pdf")

