import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import os

from read_input_file import *
from read_output_file import *
from spock_main_input import *
from find_in_read_input_order_variables import *
from eci_to_lvlh import *
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
import pickle


msis_filename = '/Users/cbv/work/spockOut/density/FM03_20190415_msis.txt'

# Read observation ECI r/v
dir_simu = '/Users/cbv/work/spockOut/density' # directory where SpOCK simu are run (input and output files)
if dir_simu[-1] != '/':
    dir_simu = dir_simu + '/'
obs_rv_filename = dir_simu + 'HD_data/nadir/cyg03.ddmi.s20190415-000000-e20190415-235959.l1.power-brcs.a21.d21.txt'

obs_rv_filename_eci = obs_rv_filename.replace('.txt','_eci.txt')
obs_rv_file = open(obs_rv_filename_eci)
read_obs_rv_file = obs_rv_file.readlines()
nb_header = 0
while (read_obs_rv_file[nb_header].split()[0] != '#START'):
    nb_header = nb_header + 1
nb_header = nb_header + 1
nb_obs = len(read_obs_rv_file) - nb_header
date_obs = []
date_obs_str = []
r_obs = np.zeros([nb_obs, 3])
v_obs = np.zeros([nb_obs, 3])
ecc_obs = np.zeros([nb_obs])
for iobs in range(nb_obs):
    date_obs_str.append( read_obs_rv_file[iobs + nb_header].split()[0] )
    date_obs.append( datetime.strptime(date_obs_str[-1], "%Y-%m-%dT%H:%M:%S" ) )
    r_obs[iobs, 0] = np.float( read_obs_rv_file[iobs + nb_header].split()[1] ) 
    r_obs[iobs, 1] = np.float( read_obs_rv_file[iobs + nb_header].split()[2] ) 
    r_obs[iobs, 2] = np.float( read_obs_rv_file[iobs + nb_header].split()[3] ) 
    v_obs[iobs, 0] = np.float( read_obs_rv_file[iobs + nb_header].split()[4] ) 
    v_obs[iobs, 1] = np.float( read_obs_rv_file[iobs + nb_header].split()[5] ) 
    v_obs[iobs, 2] = np.float( read_obs_rv_file[iobs + nb_header].split()[6] ) 

# Run SpOCK: initial r/v is given by observations + ensemble with std given by x_sigma, y_sigma, etc
date_obs_start_str = date_obs_str[0]
date_obs_start= datetime.strptime(date_obs_start_str, "%Y-%m-%dT%H:%M:%S")
date_obs_end_str = date_obs_str[-1]
date_obs_end= datetime.strptime(date_obs_end_str, "%Y-%m-%dT%H:%M:%S")

# Read SpOCK
isc = 0
var_in, var_in_order = read_input_file(msis_filename)
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')]; 
var_to_read = ["position", "velocity"]
var_out, var_out_order = read_output_file( output_file_path_list[isc] + output_file_name_list[isc], var_to_read )
date_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date')])
date_datetime_round_sec_spock = np.array(var_out[find_in_read_input_order_variables(var_out_order, 'date_datetime_round_sec')])
r_spock = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
v_spock = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]

# Compare SpOCK and observation sample times 
nb_seconds_since_start = []
distance = []

date_start = date_obs_start
date_start_str = datetime.strftime(date_start, "%Y-%m-%dT%H:%M:%S")

index_obs_interval_start = 0 # !!!!!!!!!!! to change
index_obs_kept = []
index_spock_same_date_as_obs_pid = []
iobs = np.where(np.array(date_obs) == date_start)[0][0]
#ipdb.set_trace()
print 'iobs', iobs
while iobs < nb_obs:
    if date_obs[iobs] > date_datetime_round_sec_spock[-1]:
        break
    else:
        if len(index_spock_same_date_as_obs_pid) == 0:
            first_obs = iobs
        if len(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0]) != 0:#can be = 0 if an observation is missing at that time
            index_spock_same_date_as_obs_pid.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])
            index_obs_kept.append(iobs)
            iobs = iobs + 60
        else: # find next obs
            while len(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0]) == 0:
                iobs = iobs + 1
            index_spock_same_date_as_obs_pid.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])

            index_obs_kept.append(iobs)
            iobs = iobs + 60




n = len(index_spock_same_date_as_obs_pid) #!!!!!!!!!! j-index_interval[iinter]

# Compare SpOCK and data
date_datetime_round_sec_spock_ok_pid = date_datetime_round_sec_spock[index_spock_same_date_as_obs_pid]
r_spock_ok_pid = np.zeros([n, 3])
r_spock_ok_pid[:, 0] = r_spock[index_spock_same_date_as_obs_pid, 0]
r_spock_ok_pid[:, 1] = r_spock[index_spock_same_date_as_obs_pid, 1]
r_spock_ok_pid[:, 2] = r_spock[index_spock_same_date_as_obs_pid, 2]
v_spock_ok_pid = np.zeros([n, 3])
v_spock_ok_pid[:, 0] = v_spock[index_spock_same_date_as_obs_pid, 0]
v_spock_ok_pid[:, 1] = v_spock[index_spock_same_date_as_obs_pid, 1]
v_spock_ok_pid[:, 2] = v_spock[index_spock_same_date_as_obs_pid, 2]

distance_lvlh_pid_sub = []
iperiod_find = 0
for i in range(n):
    distance_here = r_obs[index_obs_kept][i, :] - r_spock_ok_pid[i, :]
    distance_lvlh_pid_sub.append( eci_to_lvlh(r_spock_ok_pid[i, :], v_spock_ok_pid[i, :], distance_here)[0] ) #[0]: along-track
    nb_seconds_since_start.append((date_datetime_round_sec_spock_ok_pid[i] - date_start).total_seconds())

distance_lvlh_pid_sub = np.array(distance_lvlh_pid_sub)
nb_seconds_since_start = np.array(nb_seconds_since_start)
pickle.dump([distance_lvlh_pid_sub, nb_seconds_since_start, date_start], open("pickle/distance_lvlh_nb_seconds_since_start_date_start_msis.pickle", "w"))    
# fig, ax = plt.subplots()
# ax.plot(distance_lvlh_pid_sub)
