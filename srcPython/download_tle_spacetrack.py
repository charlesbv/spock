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
from get_prop_dir import *

# # IF USED CHOOSES DIRECTLY A LINK
# # PARAMETERS TO SET
# run_dir = "run.cygnss"
# link = "https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/NORAD_CAT_ID/41884,%2041885,%2041886,%2041887,%2041888,%2041889,%2041890,%2041891/orderby/TLE_LINE1%20ASC/format/tle"

# # ALGORITHM
# os.chdir("../" + run_dir + "/input/tle/")
# os.system("wget  --post-data='identity=cbv@umich.edu&password=cygnssisawesome' --cookies=on --keep-session-cookies --save-cookies=cookies.txt 'https://www.space-track.org/ajaxauth/login' -olog")
# os.system("wget --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt " + link)
# os.chdir("../../../srcPython")

# raise Exception
# IF USED CHOOSES A START AND AN END DATES
# PARAMETERS TO SET
run_dir = "run.cygnss"
date_start = "2016-12-20"
date_end = "2017-03-20"
link = "https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/" + date_start + "--" + date_end + "/NORAD_CAT_ID/41884,%2041885,%2041886,%2041887,%2041888,%2041889,%2041890,%2041891/orderby/TLE_LINE1%20ASC/format/tle"
#link = "https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/" + date_start + "--" + date_end + "/NORAD_CAT_ID/41891/orderby/TLE_LINE1%20ASC/format/tle"

# ALGORITHM
os.chdir("../" + run_dir + "/input/tle/")
os.system("wget  --post-data='identity=cbv@umich.edu&password=cygnssisawesome' --cookies=on --keep-session-cookies --save-cookies=cookies.txt 'https://www.space-track.org/ajaxauth/login' -olog")
os.system("wget --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt " + link + " -O " + date_start + "--" + date_end + ".txt")
os.system("cp " + date_start + "--" + date_end + ".txt ../../../srcPython/cygnss")
os.chdir("../../../srcPython")

