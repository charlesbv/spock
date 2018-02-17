# This script downloads the TLEs at the epoch specified by tle_epoch for NORAD ID norad_id given in input (the TLE right "before" (older than) tle_epoch). The input are the NORAD ID, and tle epoch. These are given as argumeents of this script. If you add a fourth argument called 'latest_tle', it does not consider tle_epoch to download the TLE but download the latest TLE from space-track.org. The name of the TLE will still be norad_id + "_" + tle_epoch + ".txt".
# This script is similar to cygnss_tle.py except that cygnss_tle.py downloads TLEs specifically for CYGNSS and downloads the TLEs for all CYGNSS. Here it's only for on sc.
# Assumptions:
## - tle_epoch must have the format "YYYY-MM-DD"

import os
from datetime import datetime, timedelta
import sys

norad_id = str(sys.argv[1])
tle_epoch = sys.argv[2] # date of the TLEs to get # "2017-01-01" 

# download TLEs of this NORAD ID from space-track for the epoch tle_epoch 
if ( len(sys.argv) > 3 ) : # if the argument 'latest_tle' has been called than do not consider tle_epoch to download the TLE but download the latest TLE from space-track.org. The name of the TLE will still be norad_id + "_" + tle_epoch + ".txt"
    if ( sys.argv[3] == 'latest_tle' ):
        ## put NORAD ID in link
        link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/NORAD_CAT_ID/" + norad_id
        itle = 0
        link_spacetrack = link_spacetrack + "," + norad_id

        ## Take all TLEs from tle_epoch minus 10 days to tle_epoch. We do that because we want the tle right older than tle_epoch. so we get all these tles and save only the most recent. 10 days is arbitrary but it's to make sure that at least one tle was published (usually they are published about every day for LEO so technically don't need to go 10 days back)
        link_spacetrack = link_spacetrack +  "/format/tle/"
        ## Download with this link
        os.system("wget  --post-data='identity=cbv@umich.edu&password=cygnssisawesome' --cookies=on --keep-session-cookies --save-cookies=cookies.txt 'https://www.space-track.org/ajaxauth/login' -olog > /dev/null 2>&1")
        name_tle = norad_id + "_" + tle_epoch  + ".txt"
        os.system("wget --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt " + link_spacetrack + " -O " + name_tle + " > /dev/null 2>&1")

else:
    ## put NORAD ID in link
    link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle/NORAD_CAT_ID/" + norad_id
    itle = 0
    link_spacetrack = link_spacetrack + "," + norad_id

    ## Take all TLEs from tle_epoch minus 10 days to tle_epoch. We do that because we want the tle right older than tle_epoch. so we get all these tles and save only the most recent one. 10 days is arbitrary but it's to make sure that at least one tle was published (usually they are published about every day for LEO so technically don't need to go 10 days back)
    tle_epoch_minus_ten_days = datetime.strftime(datetime.strptime(tle_epoch, "%Y-%m-%d") - timedelta(days = 10), "%Y-%m-%d")
    link_spacetrack = link_spacetrack + "/EPOCH/" + tle_epoch_minus_ten_days  + "--" + tle_epoch + "/format/tle/"
    ## Order by NORAD ID (each NORAD has many TLEs because it's all TLEs since tle_epoch)
    link_spacetrack = link_spacetrack + "orderby/NORAD_CAT_ID/"
    ## Download with this link
    os.system("wget  --post-data='identity=cbv@umich.edu&password=cygnssisawesome' --cookies=on --keep-session-cookies --save-cookies=cookies.txt 'https://www.space-track.org/ajaxauth/login' -olog > /dev/null 2>&1")
    name_tle = norad_id + "_" + tle_epoch  + "_temp.txt"
    os.system("wget --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt " + link_spacetrack + " -O " + name_tle + " > /dev/null 2>&1")
    ## This TLE file contains too many TLEs. Indeed, we only want the most recent TLE  of the list (so the TLE right "before" (older than) tle_epoch). This file is arranged by epoch. So take only the last TLE
    tle_spacetrack = open(name_tle)
    read_tle_spacetrack = tle_spacetrack.readlines()
    nb_tle = len(read_tle_spacetrack) / 2
    name_new_tle = norad_id + "_" + tle_epoch  + ".txt"
    new_tle_spacetrack_name =  name_new_tle
    new_tle_spacetrack = open(name_new_tle, "w+")
    read_new_tle_spacetrack = new_tle_spacetrack.readlines()
    itle = 0
    cou = 0

    while itle < nb_tle:
        current_norad = read_tle_spacetrack[itle*2+1].split()[1]
        new_norad = current_norad
        while new_norad == current_norad: # this is from another script where there could be more than one sc
            itle = itle + 1
            if (itle < nb_tle):
                new_norad = read_tle_spacetrack[itle*2+1].split()[1]
            else:
                break
        if itle < nb_tle+1:
            print >> new_tle_spacetrack, read_tle_spacetrack[(itle-1)*2].replace("\r", "").replace("\n", "")
            print >> new_tle_spacetrack, read_tle_spacetrack[(itle-1)*2+1].replace("\r", "").replace("\n", "")


    tle_spacetrack.close()
    new_tle_spacetrack.close()

    # rm the spacetrack tle file with all tles (the one that has all tles between tle_epoch minus a month and tle_epoch)
    os.system("rm -f " + name_tle)
    os.system("rm -f login cookies.txt log")
