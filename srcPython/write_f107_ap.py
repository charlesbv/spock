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
# This script creates Ap and F107 file
from datetime import datetime, timedelta

date_start_str = "2017-04-04T00:00:00"
date_end_str =  "2017-04-05T00:00:00"

date_start = datetime.strptime( date_start_str, "%Y-%m-%dT%H:%M:%S")
date_start_minus_41_days = date_start - timedelta(days = 41) 
date_end = datetime.strptime( date_end_str, "%Y-%m-%dT%H:%M:%S")
date_end_plus_41_days = date_end + timedelta(days = 41) 


# F107 file
nb_hour_f107 = (int) (( date_end_plus_41_days - date_start_minus_41_days ).total_seconds()/3600./24) + 2
filename_f107 = "f107_" + date_start_str.replace(":", "_") + "_to_" + date_end_str.replace(":", "_") + ".txt"
file_f107 = open(filename_f107 , "w+")
print >> file_f107, "YEAR DOY HR 1"
for ihour in range(nb_hour_f107):
    date_now = date_start_minus_41_days + timedelta(days = ihour)
    date_now_str = datetime.strftime( date_now, "%Y %j %H")
    print >> file_f107, date_now_str, 100.0
print >> file_f107, "</pre><hr><HR>"
file_f107.close()
    
# AP file
nb_hour_ap = (int) (( date_end - date_start ).total_seconds()/3600./24) + 2
filename_ap = "ap_" + date_start_str.replace(":", "_") + "_to_" + date_end_str.replace(":", "_") + ".txt"
file_ap = open(filename_ap , "w+")
print >> file_ap, "YEAR DOY HR 1"
for ihour in range(nb_hour_ap):
    date_now = date_start + timedelta(days = ihour)
    date_now_str = datetime.strftime( date_now, "%Y %j %H")
    print >> file_ap, date_now_str, 15
print >> file_ap, "</pre><hr><HR>"
file_ap.close()
    

