# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
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
