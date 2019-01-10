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

# This scriptdownloads L1 netcdf files from
# ftp://podaac.jpl.nasa.gov/allData/cygnss/L1/v2.1
# inputs: date_start, date_stop, cygfm, save_dir


import os
from datetime import datetime, timedelta

date_start = '2018-09-18' #'2018-10-26'
date_stop = '2018-09-25'#'2018-10-29'
podaac_path = "ftp://podaac.jpl.nasa.gov/allData/cygnss/L1/v2.1/"
save_dir = '/Volumes/Seagate_Expansion_Drive/netcdf/' 

#cygfm = 1 # 1 to 8
# '/Users/cbv/cygnss/netcdf/' # with slash at the end. the year
# will be added as a sub-directory


date_start_date = datetime.strptime(date_start, "%Y-%m-%d")
date_stop_date = datetime.strptime(date_stop, "%Y-%m-%d")


date_here = date_start_date
while date_here <= date_stop_date:
    yy = date_here.strftime('%Y')
    doy = date_here.strftime('%j')

    if (os.path.isdir(save_dir + yy + '/' + str(doy).zfill(3)) == False):
        os.system("mkdir " + save_dir + yy + '/' + str(doy).zfill(3))
    for cygfm in range(1, 9):
        if cygfm == 1:
            name_netcdf_podaac = []
            name_netcdf_local = []
            list_netcdf = []
            list_cyg = []
            filename_list =  save_dir + yy + '/'  + str(doy).zfill(3) +\
            '/index.html'
            os.system("wget " + podaac_path +
                      yy + '/'  + str(doy).zfill(3) + "/ -O " + filename_list)
            file_list = open(filename_list)
            r_file_list = file_list.readlines()
            n = len(r_file_list)
            for i in range(n):
                l = r_file_list[i]
                if '.nc</a>' in l:
                    name_netcdf_podaac.append(l.split('href=')[1].split('>')[0])
                    name_netcdf_local.append( save_dir + yy + '/'  + \
                        str(doy).zfill(3) + '/' + \
                        name_netcdf_podaac[-1].split('/')[-1][:-1] )
                    list_netcdf.append(name_netcdf_local)
                    list_cyg.append(\
                    name_netcdf_local[-1].split('cyg0')[1].split('.')[0])
        if (str(cygfm) in list_cyg) == True: # if file exists
            index_list = list_cyg.index(str(cygfm))
            # print "wget " + name_netcdf_podaac[index_list] +\
            #     " -O " + name_netcdf_local[index_list] 
            os.system( "wget " + name_netcdf_podaac[index_list] + \
                " -O " + name_netcdf_local[index_list] )
    date_here = date_here + timedelta(days = 1)
