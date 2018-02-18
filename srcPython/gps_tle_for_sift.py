# This script downloads the GPS TLEs at the epoch specified by tle_epoch (the TLE right "before" (older than) tle_epoch). The input are the tle epoch and the SpOCK run directory where to store the TLEs. These are given as argumeents of this script
# Methodology: the TLEs are taken from space-track.org. But to do so, you need to put the name of the GPS that are currently operational. Since I don't know how to find which GPS are currently operational, I use celestrack.com, which has a link that directly gives the TLEs of the currently operational GPS. But using this link only forces us to use the latest TLEs. So I use this link only to get the list of the NORAD ID of the current operational GPS. Then I input this list on space-track.org to get the TLEs for this list for the tle_epoch chosen by the user
# Summary of methodology:
# 1- get latest TLEs of operational GPS at celestrak
# 2- from the TLEs at celestrak, get the NORAD IDs of the operational GPS
# 3- download TLEs of these GPS from space-track for the epoch tle_epoch
# Assumptions:
## - tle_epoch must have the format "YYYY-MM-DD"

import os
from datetime import datetime, timedelta
import sys

def gps_tle_for_sift(tle_epoch, latest_tle_or_not):
    #tle_epoch = sys.argv[1] # date of the GPS TLEs to get # "2017-01-01"

    # get latest TLEs of operational GPS at celestrak     
    log_filename = "log_gps_tle_" + tle_epoch + ".txt"
    os.system('wget --no-check-certificate https://www.celestrak.com/NORAD/elements/gps-ops.txt'+ " >> " + log_filename+ " 2>&1")# > /dev/null 2>&1") !!! can't redirect in windows version

    # from the TLEs at celestrak, get the NORAD IDs of the operational GPS 
    tle_celestrak = open("gps-ops.txt")
    read_tle_celestrak = tle_celestrak.readlines()
    nb_tle = len(read_tle_celestrak) / 3
    norad_id_operational_gps = []
    for itle in range(nb_tle):
        norad_id_operational_gps.append(read_tle_celestrak[itle*3+2].split()[1])
    tle_celestrak.close()

    # rm file form celestrak
    os.system("rm gps-ops.txt"+ " >> " + log_filename)

    # download TLEs of these GPS from space-track for the epoch tle_epoch 
    ## put NORAD ID in link
    if ( latest_tle_or_not == 'latest_tle' ): # if the argument 'latest_tle' has been called than do not consider tle_epoch to download the TLE but download the latest TLE from space-track.org. The name of the TLE will still be "gps_" + tle_epoch + ".txt"
        link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/NORAD_CAT_ID/" + norad_id_operational_gps[0]
        for itle in range(nb_tle):
            link_spacetrack = link_spacetrack + "," + norad_id_operational_gps[itle]
        link_spacetrack = link_spacetrack + "/predicates/OBJECT_NAME,TLE_LINE1,TLE_LINE2/format/3le/"
        ## Download with this link
        os.system('wget --no-check-certificate  --post-data="identity=cbv@umich.edu&password=cygnssisawesome" --cookies=on --keep-session-cookies --save-cookies=cookies.txt "https://www.space-track.org/ajaxauth/login" -olog'+ " >> " + log_filename)# !!! can't redirect in windows version > /dev/null 2>&1")
        name_tle = "gps_" + tle_epoch  + ".txt"
        os.system('wget --no-check-certificate --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt ' + link_spacetrack + ' -O ' + name_tle+ " >> " + log_filename+ " 2>&1")#!!! can't redirect in windows version + " > /dev/null 2>&1")
    else:
        link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle/NORAD_CAT_ID/" + norad_id_operational_gps[0]
        for itle in range(nb_tle):
            link_spacetrack = link_spacetrack + "," + norad_id_operational_gps[itle]

        ## Take all TLEs from tle_epoch minus one month to tle_epoch. We do that because we want the tle right older than tle_epoch. so we get all these tles and save only the most recent one for each gps. one month is arbitrary but it's to make sure that at least one tle for each gps was published (usualy they are published about every day so technically don't need to go back a month back)
        tle_epoch_minus_one_month = datetime.strftime(datetime.strptime(tle_epoch, "%Y-%m-%d") - timedelta(days = 10), "%Y-%m-%d")
        link_spacetrack = link_spacetrack + "/EPOCH/" + tle_epoch_minus_one_month + "--" + tle_epoch + "/predicates/OBJECT_NAME,TLE_LINE1,TLE_LINE2/format/3le/"
        ## Order by NORAD ID (each NORAD ahs many TLEs because it's all TLEs since tle_epoch)
        link_spacetrack = link_spacetrack + "orderby/NORAD_CAT_ID/"
        ## Download with this link
        os.system('wget --no-check-certificate  --post-data="identity=cbv@umich.edu&password=cygnssisawesome" --cookies=on --keep-session-cookies --save-cookies=cookies.txt "https://www.space-track.org/ajaxauth/login" -olog'+ " >> " + log_filename)#!!! can't redirect in windows version > /dev/null 2>&1")
        name_tle = "gps_" + tle_epoch  + "_temp.txt"
        os.system('wget --no-check-certificate --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt ' + link_spacetrack + ' -O ' + name_tle+ " >> " + log_filename+ " 2>&1")#!!! can't redirect in windows version + " > /dev/null 2>&1")
        ## This TLE file contains too many TLEs for each GPS. Indeed, for each GPS, we only want the most recent TLE  of the list (so the TLE right "before" (older than) tle_epoch). But this file is arranged by epoch (in addition to be arranged by NORAD ID). So for each NORAD ID, take only the last TLE
        tle_spacetrack = open(name_tle)
        read_tle_spacetrack = tle_spacetrack.readlines()
        nb_tle = len(read_tle_spacetrack) / 3
        name_new_tle = "gps_" + tle_epoch  + ".txt"
        new_tle_spacetrack_name = name_new_tle
        new_tle_spacetrack = open(name_new_tle, "w+")
        read_new_tle_spacetrack = new_tle_spacetrack.readlines()
        itle = 0
        cou = 0
        while itle < nb_tle:
            current_norad = read_tle_spacetrack[itle*3+2].split()[1]
            new_norad = current_norad
            while new_norad == current_norad:
                itle = itle + 1
                if (itle < nb_tle):
                    new_norad = read_tle_spacetrack[itle*3+2].split()[1]
                    norad_id_operational_gps.append(read_tle_spacetrack[itle*3+2].split()[1])
                else:
                    break
            if itle < nb_tle+1:
                print >> new_tle_spacetrack, read_tle_spacetrack[(itle-1)*3].replace("\r", "").replace("\n", "")
                print >> new_tle_spacetrack, read_tle_spacetrack[(itle-1)*3+1].replace("\r", "").replace("\n", "")
                print >> new_tle_spacetrack, read_tle_spacetrack[(itle-1)*3+2].replace("\r", "").replace("\n", "")


        tle_spacetrack.close()
        new_tle_spacetrack.close()

        # rm the spacetrack tle file with all tles (the one that for each gps has all tles between tle_epoch minus a month and tle_epoch)
        os.system("rm -f " + name_tle+ " >> " + log_filename)
        os.system("rm -f login cookies.txt log"+ " >> " + log_filename)
    return
