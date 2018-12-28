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
# This script downloads the inputs for SIFT from the CYGNSS server to the SOC computers (kiosk account)

import sys
import os

start_date = sys.argv[1] # "2017-03-08T12:00:00"
end_date = sys.argv[2] #"2017-03-09T12:00:00"


path_spock_on_cygnss_server = '/data/cygnss/tools/spock/'
run_dir = 'run.cygnss.all' 

if ( 'T' in start_date ) == False:
    start_date = start_date + "T00:00:00"
if ( 'T' in end_date ) == False:
    end_date = end_date + "T00:00:00"

start_date_no_time = start_date[0:10]

main_input_file_name = path_spock_on_cygnss_server + run_dir + '/input/main_input/' + 'spock_spec_' + start_date + '_' + end_date + '.txt'
filename_gps_tle = path_spock_on_cygnss_server + run_dir + '/input/tle/constellation_gps_tle/gps_' + start_date_no_time + '.txt'


nb_sc = 8
for isc in range(nb_sc):
    output_file_path = "/data/cygnss/tools/spock/run.cygnss.all/output/spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_") + "/spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_") + str(isc+1) + "/"
    output_file_name = "spock_spec_start_" + start_date.replace(":","_") + "_end_" + end_date.replace(":","_")  + str(isc+1) + ".txt"
    sift_files = [output_file_path + "interpolated_position_LLA_ECEF_" + output_file_name,
                  output_file_path + "specular_" + output_file_name,
                  main_input_file_name,
                  filename_gps_tle] 
    nb_file = len(sift_files)
    input_sift_folder = ["sat_positions/", "spec_positions/", "main_input/", "tle/constellation_gps_tle/"]
    for ifile in range(nb_file):
        input_sift_folder[ifile] = '../input_sift/' + start_date.replace(":","_") + "_end_" + end_date.replace(":","_") + "/input/" + input_sift_folder[ifile]
    sift_files_str = " ".join(sift_files)
    ## Download
#    print "scp ./ cygnss-sftp-1.engin.umich.edu:" + sift_files_str
    print "\nXXXXXXXXXX"
    for ifile in range(nb_file):
#        print "scp " + input_sift_folder[ifile] + " cygnss-sftp-1.engin.umich.edu:" + sift_files[ifile]
#        print sift_files[ifile] +  ", " + input_sift_folder[ifile]
        sftp.get(sift_files[ifile],input_sift_folder[ifile])
#        os.system("scp " + input_sift_folder[ifile] + " cygnss-sftp-1.engin.umich.edu:" + sift_files_str[ifile])
