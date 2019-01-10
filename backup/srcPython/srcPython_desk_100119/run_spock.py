from cadre_read_last_tle import *
from get_prop_dir import *
import matplotlib.gridspec as gridspec
from read_input_file import *
from matplotlib.colors import LogNorm
import pickle
from eci_to_lvlh import *
import sys
import fileinput
import time
import datetime
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess


# os.chdir("../run_cygnss")
# os.system("/usr/local/bin/mpirun -np 10 run_moat cygnss_27_30_100.txt")
# os.chdir("../srcPython")
# os.system("python concatenate_proc.py run_cygnss cygnss_27_30_100.txt")
# os.system("python distance_ensemble_to_main_sc.py run_cygnss cygnss_27_30_100.txt 1")
# os.system("python plot_ensembles.py run_cygnss cygnss_27_30_100.txt rho along cross radial angle f107 f107a ap")


os.chdir("../run_cygnss")
os.system("/usr/local/bin/mpirun -np 10 run_moat cygnss_27_30_1000.txt")
os.chdir("../srcPython")
os.system("python concatenate_proc.py run_cygnss cygnss_27_30_1000.txt")
os.system("python distance_ensemble_to_main_sc.py run_cygnss cygnss_27_30_1000.txt 1")
os.system("python plot_ensembles.py run_cygnss cygnss_27_30_1000.txt rho along cross radial angle f107 f107a ap")
