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
import os 

# On Saturday July 2nd, these were the simulations that were aldready done:
#run_names = ['SCION_4_0_deg_s_f107_120_ap_15.txt', 'SCION_5_0_deg_s_f107_120_ap_15.txt', 'SCION_0_8_deg_s_f107_120_ap_15.txt', 'SCION_0_8_deg_s_f107_200_ap_80.txt', 'SCION_0_8_deg_s_f107_90_ap_7.txt', 'SCION_2_5_deg_s_f107_120_ap_15.txt']

run_names = [ 'SCION_0_8_deg_s_f107_200_ap_80.txt']
n = len(run_names)
for i in range(n):
    print run_names[0:i+1]
    # os.chdir("../../")
    # os.system("time /usr/local/bin/mpirun -np 10 run_moat " + run_names[i])
    # os.chdir("code/python_propagator")
    os.system("python concatenate_proc.py " + run_names[i])
    os.system("python distance_ensemble_to_main_sc.py " + run_names[i] + " 1")








