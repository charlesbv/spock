# This script gathers all steps necessary to compute the density correction factor to apply to NRLMSIS
# Inputs
# - eng_pvt (ECEF r/v) and eng_adcs (attitude quaternions) of CYGNSS provided by SwRI

# 


import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
import os
from convert_cygnss_obs_ecef_to_eci import *


# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
eng_pvt = 'FM4_20171216_eng_pvt_query-13525.txt'
eng_adcs = 'FM4_20171216_eng_adcs_query-13526.txt'
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT


# cygnss_convert_swri_att_to_spock.py changes the format of eng_pvt and eng_adcs (readible by SpOCK)
if ('/' in eng_pvt) == False: # for some reason need ./ if in current directory
    eng_pvt = './' + eng_pvt
    eng_adcs = './' + eng_adcs
os.system("python cygnss_convert_swri_att_to_spock.py " + eng_pvt + " " + eng_adcs)

# convert_cygnss_obs_ecef_to_eci converts ECEF r/v (eng_pvt by SwRI) into ECI r/v
eng_pvt_no_root = eng_pvt.split('/')[-1]
eng_pvt_only_root = '/'.join(eng_pvt.split('/')[:-1]) + '/'
spock_pvt_ecef = eng_pvt_only_root + 'spock_' + eng_pvt_no_root
convert_cygnss_obs_ecef_to_eci(spock_pvt_ecef)
spock_pvt_eci = spock_pvt_ecef.replace('.txt','_eci.txt')

