# This script generates a SpOCK simulation and call cygnss_lakes for that simulation

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT                                                                                                                                
date_start_simu = '2019-08-28'
date_stop_simu = '2020-02-28'
# end of  PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT                                                                                     

# Algorithm
from datetime import datetime, timedelta
import sys
import os
import ipdb
sys.path.append("/Users/cbv/work/spock/srcPython")
from read_input_file import *
from cygnss_read_spock_spec_bin import *
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib import pyplot as plt
from cygnss_lakes import *

date_start_simu_date = datetime.strptime(date_start_simu, "%Y-%m-%d")
date_stop_simu_date = datetime.strptime(date_stop_simu, "%Y-%m-%d")
nday = (date_stop_simu_date - date_start_simu_date).days
date_arr_simu = np.array([date_start_simu_date + timedelta(days=i) for i in np.arange(0, nday, 1)])
visit_time = []; visit_which_sp = []; visit_which_sc = []
for iday in range(nday):
    date_start_date = date_arr_simu[iday]
    date_stop_date = date_start_date + timedelta(days = 1)
    date_start = datetime.strftime(date_start_date, "%Y-%m-%d") + "T00:00:00"
    spock_input_filename = date_start[:10].replace('-','') + '.txt'
    print '\niday', iday, nday-1,  spock_input_filename, str(datetime.now())[0:19]
    # Create the SpOCK main input file
    os.system("cp 082819.txt " + spock_input_filename)
    date_stop = datetime.strftime(date_stop_date, "%Y-%m-%d") + "T00:00:00" 
    with open(spock_input_filename) as f:
        lines = f.readlines()
        lines[1] = date_start + '\n'
        lines[2] = date_stop + '\n' 
    with open(spock_input_filename, "w") as f:
        f.writelines(lines)
    
    # Run SpOCK to predict the specular point positions
    os.system("mpirun -np 4 spock_dev " + spock_input_filename)
    os.system("mpirun -np 4 spec " + spock_input_filename + " -lon=0 -rot=0 -min")

    # Call cygnss_lakes to compute the revisit times
    visit_time_day, visit_which_sp_day, visit_which_sc_day = cygnss_lakes(spock_input_filename)
    print visit_time_day
    visit_time.append(visit_time_day); visit_which_sp.append(visit_which_sp_day); visit_which_sc.append(visit_which_sc_day)
