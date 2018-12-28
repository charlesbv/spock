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
from os import listdir
from os.path import isfile, join

dir_to_open = get_prop_dir(1) + "run_aerie/input/density/density_NRLMSIS00e/"
file_to_modify_arr = [f for f in listdir(dir_to_open) if isfile(join(dir_to_open, f))]

for ifile in range(len(file_to_modify_arr)):
    file_to_modify = open(dir_to_open + file_to_modify_arr[ifile])
    file_out = open(dir_to_open + 'modif_' + file_to_modify_arr[ifile], "w+")
    read_file_to_modify = file_to_modify.readlines()
    n_header = 3
    n = len(read_file_to_modify) - n_header - 1 #  -1 to remove last line '#ENDOFFILE'
    print >> file_out, "#BEGINNINGOFHEADER" 
    print >> file_out, "#ENDOFHEADER" 
    print >> file_out, "YEAR DOY HR    1" 
    for istep in range(n):
        print >> file_out, read_file_to_modify[istep + n_header].split()[0] + " " + read_file_to_modify[istep + n_header].split()[1] + " 0 " + read_file_to_modify[istep + n_header].split()[2]
    print >> file_out, "#ENDOFFILE" 
    file_out.close()
    file_to_modify.close()
