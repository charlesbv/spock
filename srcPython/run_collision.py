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
