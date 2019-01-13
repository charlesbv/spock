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

plt.ion()
plt.isinteractive()


os.chdir("../../")
os.system("/usr/local/mpi/bin/mpirun -np 10 run_moat SCION.txt")
os.chdir("./code/python_propagator/")
os.system("python concatenate_proc.py SCION.txt")
os.system("python distance_ensemble_to_main_sc.py SCION.txt 1")

