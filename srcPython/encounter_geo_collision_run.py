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

# This script runs the main input file included in the file filename_run_list. The main input file, as ell as all the necessary other SpOCK input file, needs to be alreeady created. This script only creates the job files and submit them. This script was inspired from the script  main_ensemble_f107_ap.py

# cbv modifed this script in july 2018 for paper CA3 (CYGNSS collision avoidance)
# ASSUMPTIONS
# - see section #PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT


import os.path
import sys
sys.path.append("/home1/cbussy/Code/spock/srcPython")

import time
from spock_job_file_new_july18 import *
import numpy as np
import os



#PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
folder_tar = 'tar_fm07_nov01_30000ens' # where to put all tar ru (the output fodler of SpOCK are tarred then poved here). Put the '/' or don't put, doesnt matter
filename_run_list = "missing_runlist_FM07_nov01_30000ens.txt"#runlist_FM07_nov01_30000ens.txt"
walltime_in = "02:00:00"
select_in = 23
ncpus_in = 24
only_devel_runs = 0 # if this is set to 1 then all simy will be run in devel mode. this script runs a simu and check every minute if the simu is done. if the simu is done, then it runs the new simu. and so on. this is because only one devel simu can be run at the same time. like this all simu are run in devel mode consecutively in the minimum amount of time
# end of parameters to set
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT  

# Algorithm
## Read the filename_run_list to get all the simu to run
file_run_list = open(filename_run_list)
run_list = []
read_file_run_list = file_run_list.readlines()
nb_run_temp = len(read_file_run_list)
for irun in range(nb_run_temp):
    run_list.append(read_file_run_list[irun].split()[0])
nb_run = len(run_list)


# if folder_tar doesnt exist then create it
if os.path.isdir(folder_tar) == False:
    os.system("mkdir " + folder_tar)

for irun in range(nb_run): # !!!!!!! should be nb_run
    main_input_filename = run_list[irun]
    print irun
    if only_devel_runs != 1:
        if irun == 0:
            filename_job = spock_job_file_new_july18(main_input_filename , "low",walltime_in,select_in, ncpus_in, folder_tar) # creates the job file for this run # s!!! should be devel queue here
        else:
            filename_job = spock_job_file_new_july18(main_input_filename, "low", walltime_in,select_in, ncpus_in, folder_tar) # creates the job file for this run   
        os.system("qsub "+ filename_job)
    else:

        filename_job = spock_job_file_new_july18(main_input_filename , "devel",walltime_in,select_in, ncpus_in, folder_tar) # creates the job file for this run 
        os.system("qsub "+ filename_job)    
        output_file_result_spock_run =  main_input_filename.replace(".txt", "_out.tgz")
        print output_file_result_spock_run
        while not os.path.exists(output_file_result_spock_run): # wait until previous simu is done (if simu is done then the file output_file_result_spock_run is created)
            time.sleep(60*3)

    # if irun == 1:
    #     raise Exception

    
# # below: every 3 minutes, check which simu is done (a simu is done is  main_input_filename.replace(".txt", "_out.tgz") exists). When all simu are done, then it performs tasks
# ## first wait the amount of time you think at le
# all_simu_done = 0
# while all_simu_done == 0:

#     for irun in range(nb_run): # !!!!!!! should be nb_run
        
