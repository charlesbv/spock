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
from get_prop_dir import *
import sys
import fileinput
import time
import datetime
import numpy as np
import pylab as P
import os
import subprocess
from read_input_file import *
from mpi4py import MPI
comm = MPI.COMM_WORLD
iProc_python = comm.Get_rank()
nProcs_python = comm.Get_size()


# PARAMTERS TO SET
nb_ensemble = 
nProcs = 240
nb_ensemble_per_iProc_python = 


for i in range(iProc_python):
    iStart = iStart + nb_ensemble_per_iProc_python
    if ((i < nb_ensemble_to_plot_left) & (iProc_python > 0)):
        iStart = iStart + 1
iEnd = iStart+nb_ensemble_per_iProc_python
if (iProc_python < nb_ensemble_to_plot_left):
    iEnd = iEnd + 1
nb_ensemble_per_proc = int(np.floor(nb_ensembles / nProcs))
nb_ensembles_ok = nb_ensemble_per_proc * nProcs

