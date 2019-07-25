# This script:
# 1- runs SpOCK with all 8 CYGNSS satellites (no GPS) between date_start and date_stop
# 2- runs cygnss_camp2ex.py to output the lon/lat/alt and solar zenith angle
# To run it:
# python cygnss_main_camp2ex.py date_start date_stop
# where date_start and date_stop have the format YYYY-mm-dd
# ASSUMPTIONS:
# - see sections PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
# - none
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT



import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
from spock_main_input import *

date_start = sys.argv[1] #'2019-07-21'#sys.argv[1]
date_stop = sys.argv[2] #'2019-07-22'#sys.argv[2]

print date_start, date_stop
# Download the CYGNSS TLEs for date_start
os.system("cygnss_tle_web.py " +   date_start)
cygnss_tle_filename = 'cygnss_' + date_start + '.txt'

# Run SpOCK between date_start and date_stop
date_start_with_time = date_start + 'T00:00:00'
date_stop_with_time = date_stop + 'T00:00:00'
dt_simu = 1.
order_gravity = 20 # !!!!!!! put 20
forces = 'drag moon_gravity sun_gravity'
main_input_filename = 'camp2ex.txt'
spock_main_input(
    main_input_filename,
    # for TIME section
    date_start_with_time,
    date_stop_with_time,
    dt_simu,
    # for SPACECRAFT section
    8,
    '0',
    29.,
    'cygnss_geometry_2016_acco09.txt',
    # for ORBIT section
    cygnss_tle_filename,
    # for FORCES section
    order_gravity,
    forces,
    'static', # !!!!!!!!make sure the static values are consistent with the current trend of f107 and ap
    # for OUTPUT section
    '~/cygnss/camp2ex/out/out',
    dt_simu,
    # for ATTITUDE section
    "nadir",
    # for GROUNDS_STATIONS section
    "0",#"my_ground_stations.txt"
    # for SPICE section
    '/Users/cbv/cspice/data',
    # for DENSITY_MOD section
    1
)
os.system("mpirun -np 4 spock_grav_read_bin_earth_map " + main_input_filename) # !!!!!!! make sure this is the latestest executable (spock on my laptop...)

# Run cygnss_camp2ex.py to output the lon/lat/alt and solar zenith angle
os.system('python cygnss_camp2ex.py ' + main_input_filename)

# Send the files to Roman Kowch at JPL
output_dir = 'Users/cbv/cygnss/camp2ex/out/' + date_start.replace('-','') + '_to_' + date_stop.replace('-','')
output_dir_no_path = date_start.replace('-','') + '_to_' + date_stop.replace('-','')
os.chdir('out')
filename_out = '/Users/cbv/cygnss/camp2ex/toshare/' + output_dir_no_path + '.tgz'
os.system('tar -zcvf ' + filename_out + ' ' +  output_dir_no_path)

os.system('/usr/bin/uuencode ' + filename_out + ' ' + filename_out + ' | /usr/bin/mail -s "me5 CYGNSS positions by SpOCK - Camp2Ex - ' + date_start + ' to ' + date_stop + '" cbv@umich.edu')

#os.system('/usr/bin/uuencode ' + filename_out + ' ' + filename_out + ' | /usr/bin/mail -s "CYGNSS positions by SpOCK - Camp2Ex - ' + date_start + ' to ' + date_stop + '" roman.s.kowch@nasa.gov')

#os.system('scp -p /Users/cbv/cygnss/camp2ex/' + filename_out.replace(".txt", ".tgz") + ' cygnss-sftp-1.engin.umich.edu:/data/temp/piston')

os.chdir('../')
