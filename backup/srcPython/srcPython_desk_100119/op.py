import datetime
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os # to exectute terminal commands from python
from  get_prop_dir import *

input_filename = "f107_for_81daverage_20130406_to_20130628.txt"
input_filename = get_prop_dir(2) + "input/density/density_NRLMSIS00e/" +  input_filename
file_to_read = open(input_filename, "r")
read_file = file_to_read.readlines()
n_h = 8
n = len(read_file) - n_h - 15
f107 = np.zeros(n)
for i in range(n):
    f107[i] = np.float(read_file[i + n_h].split()[3])

nb_day = 2
print np.mean(f107[nb_day : nb_day+81])


file_to_read.close()
