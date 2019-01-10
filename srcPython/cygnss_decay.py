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

# This script predicts the semi-major axis of the CYGNSS sc 6 months ahead (starting today at midnight). It initializes the orbits of the 8 sc with the latest TLE available at space-track.org and run SpOCK to make predict the sma. 
# To run it:
# python cygnss_decay.py
# Assumptions:
# - see assumptions of script "spock_cygnss_spec_set_spock_simu_parameters.py"

import os
from datetime import datetime, timedelta

# PARAMETERS TO SET BEFORE RUNNING THE SCRIPT
# -> none

# ALGORITHM
# Run SpOCK from the current date to current date + 6 months (the script spock_cygnss_spec.py also downloads the latest TLEs). Note that you might need to set up the path to mpirun in spock_cygnss_spec.py (for more details, see header in spock_cygnss_spec.py)
date_today = datetime.today()
date_today_str = datetime.strftime(date_today, "%Y-%m-%d")
six_month_in_days = 190
date_today_plus_6_month = datetime.today() + timedelta(days = six_month_in_days)
date_today_plus_6_month_str = datetime.strftime(date_today_plus_6_month, "%Y-%m-%d")
# Nadir
in_name_nadir = 'spock_spec_start_' + date_today_str + "T00_00_00" + '_end_' + date_today_plus_6_month_str + "T00_00_00" + '_nadir.txt'
#os.system("python spock_cygnss_spec_set_spock_simu_parameters.py " + date_today_str + " " + date_today_plus_6_month_str + " order=10 dt=5 geo=cygnss_geometry_2016.txt att=nadir in_name=" + in_name_nadir )
# High drag
in_name_high_drag = 'spock_spec_start_' + date_today_str + "T00_00_00" + '_end_' + date_today_plus_6_month_str + "T00_00_00" + '_high_drag.txt'
os.system("python spock_cygnss_spec_set_spock_simu_parameters.py " + date_today_str + " " + date_today_plus_6_month_str + " order=10 dt=5 geo=cygnss_geometry_2016.txt att='(82; 0; 0)(0; 0; 0)' in_name=" + in_name_high_drag )


# Plot the semi-major axis of the 8 CYGNSS as a function of time
## name of the simulation we just run with SpOCK
start_date = datetime.strftime( datetime.strptime(date_today_str, "%Y-%m-%d"), "%Y-%m-%dT%H:%M:%S" ) 
end_date = datetime.strftime( datetime.strptime(date_today_plus_6_month_str, "%Y-%m-%d"), "%Y-%m-%dT%H:%M:%S" ) 
main_input_filename = 'spock_spec_start_' + start_date.replace(":", "_") + '_end_' + end_date.replace(":", "_") + '.txt'
#os.system("python state.py run.cygnss " + main_input_filename + " plot sma")

