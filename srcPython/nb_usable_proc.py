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

# This script finds the number of processors you can use when runninng mpirun. Other solutions (such as multiprocessing.cpu_count()) return more processors than actually available from mpirun.

import os
def nb_usable_proc(path_mpirun):
    os.system(path_mpirun + " hostname > nbcore_temp")
    file_nb_core = open("nbcore_temp")
    read_file_nb_core = file_nb_core.readlines()
    nb_core = len(read_file_nb_core)
    os.system("rm -f nbcore_temp")
    return nb_core

