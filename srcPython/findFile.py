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

# This script finds the netcd file corresponidng to a FM for a given date
# Inputs:
# - date_netcdf: date of the netcdf file (dateimte object)
# - cygfm: which FM (1 to 8)
# - netcdf_dir: where to look for the file. ex:/Users/cbv/netcdf/

from datetime import datetime, timedelta
import os
from os import listdir
from os.path import isfile, join


def findFile(date_netcdf, cygfm, netcdf_dir):
    yy = datetime.strftime(date_netcdf, "%Y") 
    doy = datetime.strftime(date_netcdf, "%j").zfill(3)
    if netcdf_dir[-1] != '/':
        netcdf_dir = netcdf_dir + '/'
    path_netcdf_file = netcdf_dir + yy + '/' + doy + '/'
    list_file = [filename for filename in os.listdir(path_netcdf_file) if \
     filename.startswith('cyg0' + str(cygfm))]
    if len(list_file) == 0:
        print '***! No netcdf file for FM0' + str(cygfm) + ' on ' +\
                 str(date_netcdf) + " !***"
        return ''
    return path_netcdf_file + list_file[0]
    
