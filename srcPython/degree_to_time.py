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
# THIS SCRIPT 
# Set twelve_hour_offset to 0 if you want to start at 12 am: 0 degree is midnight, 90 degrees is 6 am, 180 degrees is 12 pm, 270 degrees is 6 pm  
# Set twelve_hour_offset to 1 if you want to start at 12 pm and not 12 am: 0 degree is 12 pm, 90 degrees is 6 pm, 180 degrees is midnight, 270 degrees is 6 am

from datetime import datetime, timedelta

def degree_to_time(var_in_deg, twelve_hour_offset):
    # var_in_deg = [0, 90, 180, 270, 225]
    # twelve_hour_offset = 1
    n = len(var_in_deg)
    var_in_time = []
    if twelve_hour_offset == 0: # 0 degree is midnight, 90 degrees is 6 am, 180 degrees is 12 pm, 270 degrees is 6 pm
        time_start = datetime.strptime("00:00:00", "%H:%M:%S")
    if twelve_hour_offset == 1: # 0 degree is 12 pm, 90 degrees is 6 pm, 180 degrees is midnight, 270 degrees is 6 am
        time_start = datetime.strptime("12:00:00", "%H:%M:%S")
    for istep in range(n):
        var_in_hour = var_in_deg[istep] / 360. * 24
        var_in_hour_date_format = time_start + timedelta(hours=var_in_hour) 
        var_in_time.append( datetime.strftime(var_in_hour_date_format, "%H:%M:%S") )

    return var_in_time
