B1;95;0cimport subprocess
import fileinput
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
import numpy as np
import sys
from read_input_file import *
from read_output_file import *
from get_prop_dir import *

plt.ion()
raise Exception
# list_simu = ['storm_static_5000ens_min_dist_coll_5m.txt',
#              'storm_dynamic_5000ens_min_dist_coll_5m.txt',
#              'storm_static_10000ens_min_dist_coll_0_5m.txt',
#              'storm_dynamic_10000ens_min_dist_coll_0_5m.txt']
list_simu = ['storm_static_7872ens_min_dist_coll_5m.txt',
             'storm_dynamic_7872ens_min_dist_coll_5m.txt']


os.chdir("../run_collision/")
nb_simu = len(list_simu)
for irun in range(nb_simu):
    print '/usr/bin/time -p /usr/local/bin/mpirun -np 12 run_moat ' + list_simu[irun]
    os.system('/usr/bin/time -p /usr/local/bin/mpirun -np 12 run_moat ' + list_simu[irun])
