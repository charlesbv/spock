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

#panel_conf = ['0panel', '1panel', '2panels', '3panels', '4panels']
panel_conf = ['1panel']


for ipanel in range(len(panel_conf)):
    os.chdir("../../")
    os.system("/usr/local/bin/mpirun -np 10 run_moat " + 'cadre_tumbling_' + panel_conf[ipanel] + '.txt')
    os.chdir("./code/python_propagator/")
    os.system("python concatenate_proc.py " + "cadre_tumbling_" + panel_conf[ipanel] + ".txt")
    os.system("python distance_ensemble_to_main_sc.py " + "cadre_tumbling_" + panel_conf[ipanel] + ".txt 1")
