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

# This function is pretty useless but just made if other users use this propagator and the python scripts I wrote

import os 

def get_prop_dir(dir_offset): # dir_offset is how deep the current directory is from the propagator directory
    propagator_directory_list = os.getcwd().split('/')[0:-dir_offset]
    propagator_directory = ""
    for i in range(len(propagator_directory_list)):
        if (i > 0):
            propagator_directory = propagator_directory + "/" + propagator_directory_list[i]
            
    propagator_directory = propagator_directory + '/'
    return propagator_directory