# This script downloads the TLE of a sc given by its NORAD ID norad_id between tle_epoch_start and tle_epoch_stop
# To run this script:
# python download_tle [norad_id] [tle_epoch_start] [tle_epoch_stop]
# note: tle_epoch must have the format YYYY-MM-DD 
# example: python download_tle 41885 2017-08-23 2017-08-29

import os
from datetime import datetime, timedelta
import sys
norad_id = sys.argv[1]
tle_epoch_start = sys.argv[2]
tle_epoch_stop = sys.argv[3] 

## put NORAD ID in link
link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle/NORAD_CAT_ID/" + norad_id
link_spacetrack = link_spacetrack + "," + norad_id
link_spacetrack = link_spacetrack + "/EPOCH/" + tle_epoch_start  + "--" + tle_epoch_stop + "/format/tle/"
link_spacetrack = link_spacetrack + "orderby/EPOCH/"
## Download with this link
os.system('wget --no-check-certificate  --post-data="identity=cbv@umich.edu&password=cygnssisawesome" --cookies=on --keep-session-cookies --save-cookies=cookies.txt "https://www.space-track.org/ajaxauth/login" -olog > /dev/null') #!!! can't redirect in windows version
name_tle = norad_id + "_" + tle_epoch_start + '_' + tle_epoch_stop  + ".txt"
os.system('wget --no-check-certificate --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt ' + link_spacetrack + ' -O ' + name_tle  + ' > /dev/null  2>&1')
os.system("rm -f login cookies.txt log  > /dev/null")
