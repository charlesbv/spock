# This script downloads the CYGNSS TLEs at the epoch specified by tle_epoch (the TLE right "before" (older than) tle_epoch). The input is the tle epoch. These are given as argumeents of this script. If you add a third argument called 'latest_tle', it does not consider tle_epoch to download the TLE but download the latest TLE from space-track.org. The name of the TLE will still be "cygnss_" + tle_epoch + ".txt"
# Assumptions:
## - tle_epoch must have the format "YYYY-MM-DD"

import os
from datetime import datetime, timedelta
import sys

def cygnss_tle_for_sift(tle_epoch, latest_tle_or_not):
    #tle_epoch = sys.argv[1] # date of the CYGNSS TLEs to get # "2017-01-01" 
    log_filename = "log_cygnss_tle_" + tle_epoch + ".txt"
    norad_id_cygnss = ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
    # download TLEs of these CYGNSS from space-track for the epoch tle_epoch 
    if ( latest_tle_or_not == 'latest_tle' ):
        ## put NORAD ID in link
        link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle_latest/ORDINAL/1/NORAD_CAT_ID/" + norad_id_cygnss[0]
        nb_sc = len(norad_id_cygnss)
        for itle in range(nb_sc):
            link_spacetrack = link_spacetrack + "," + norad_id_cygnss[itle]

        ## Take all TLEs from tle_epoch minus 10 days to tle_epoch. We do that because we want the tle right older than tle_epoch. so we get all these tles and save only the most recent one for each cygnss. 10 days is arbitrary but it's to make sure that at least one tle for each cygnss was published (usually they are published about every day so technically don't need to go 10 days back)
        link_spacetrack = link_spacetrack +  "/format/tle/"
        ## Order by NORAD ID (each NORAD has many TLEs because it's all TLEs since tle_epoch)
        ## Download with this link
        os.system('wget --no-check-certificate  --post-data="identity=cbv@umich.edu&password=cygnssisawesome" --cookies=on --keep-session-cookies --save-cookies=cookies.txt "https://www.space-track.org/ajaxauth/login" -olog'+ " >> " + log_filename)# > /dev/null 2>&1") !!! can't redirect in windows version
        name_tle = "cygnss_" + tle_epoch  + ".txt"
        os.system('wget --no-check-certificate --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt ' + link_spacetrack + ' -O ' + name_tle + " >> " + log_filename + " 2>&1")# !!! can't redirect in windows version + " > /dev/null 2>&1")

    else:
        ## put NORAD ID in link
        link_spacetrack = "https://www.space-track.org/basicspacedata/query/class/tle/NORAD_CAT_ID/" + norad_id_cygnss[0]    
        nb_sc = len(norad_id_cygnss)
        for itle in range(nb_sc):
            link_spacetrack = link_spacetrack + "," + norad_id_cygnss[itle]
        ## Take all TLEs from tle_epoch minus 10 days to tle_epoch. We do that because we want the tle right older than tle_epoch. so we get all these tles and save only the most recent one for each cygnss. 10 days is arbitrary but it's to make sure that at least one tle for each cygnss was published (usually they are published about every day so technically don't need to go 10 days back)
        tle_epoch_minus_ten_days = datetime.strftime(datetime.strptime(tle_epoch, "%Y-%m-%d") - timedelta(days = 10), "%Y-%m-%d")
        link_spacetrack = link_spacetrack + "/EPOCH/" + tle_epoch_minus_ten_days  + "--" + tle_epoch + "/format/tle/"
        ## Order by NORAD ID (each NORAD has many TLEs because it's all TLEs since tle_epoch)
        link_spacetrack = link_spacetrack + "orderby/NORAD_CAT_ID/"
        ## Download with this link
        os.system('wget --no-check-certificate  --post-data="identity=cbv@umich.edu&password=cygnssisawesome" --cookies=on --keep-session-cookies --save-cookies=cookies.txt "https://www.space-track.org/ajaxauth/login" -olog' + " >> " + log_filename)# > /dev/null 2>&1") !!! can't redirect in windows version
        name_tle = "cygnss_" + tle_epoch  + "_temp.txt"
        os.system('wget --no-check-certificate --limit-rate=100K --keep-session-cookies --load-cookies=cookies.txt ' + link_spacetrack + ' -O ' + name_tle+ " >> " + log_filename + " 2>&1") # + " >> " + log_filename)#+ " > /dev/null 2>&1") !!! can't redirect in windows version
        ## This TLE file contains too many TLEs for each CYGNSS. Indeed, for each CYGNSS, we only want the most recent TLE  of the list (so the TLE right "before" (older than) tle_epoch). This file is arranged by epoch (in addition to be arranged by NORAD ID). So for each NORAD ID, take only the last TLE
        tle_spacetrack = open(name_tle)
        read_tle_spacetrack = tle_spacetrack.readlines()
        nb_tle = len(read_tle_spacetrack) / 2
        name_new_tle = "cygnss_" + tle_epoch  + ".txt"
        new_tle_spacetrack_name =  name_new_tle
        new_tle_spacetrack = open(name_new_tle, "w+")
        read_new_tle_spacetrack = new_tle_spacetrack.readlines()
        itle = 0
        cou = 0

        while itle < nb_tle:
            current_norad = read_tle_spacetrack[itle*2+1].split()[1]
            new_norad = current_norad
            while new_norad == current_norad:
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

        # rm the spacetrack tle file with all tles (the one that for each cygnss has all tles between tle_epoch minus a month and tle_epoch)
        os.system("rm -f " + name_tle + " >> " + log_filename)
        os.system("rm -f login cookies.txt log" + " >> " + log_filename)
    return
