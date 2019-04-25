#!/Users/cbv/Library/Enthought/Canopy_64bit/User/bin/python
# This script:
# - calls cygnss_tle.py
# - create 8 TLE file - one TLE per file
# - add the name of the cygnss as the first line of the TLE -> the file is 3 lines
from datetime import datetime, timedelta
import sys
sys.path.append('/Users/cbv/work/spock/srcPyton')
import os
from norad_id_to_cygnss_id import *

tle_epoch = sys.argv[1] # date of the CYGNSS TLEs to get # "2017-01-01" 
os.system("cygnss_tle.py " + tle_epoch)
filename_all = 'cygnss_' + tle_epoch + '.txt'
file_all = open(filename_all)
read_file_all = file_all.readlines()
for itle in range(8):
    norad = read_file_all[itle*2+1].split()[1]
    fm = 'CYG' + norad_id_to_cygnss_id(norad)
    filename_fm = fm + '_' + tle_epoch + '.txt'
    file_fm = open(filename_fm, 'w')
    file_fm.writelines(fm + '\n')
    file_fm.writelines(read_file_all[itle*2])
    file_fm.writelines(read_file_all[itle*2+1])
    file_fm.close()

os.system('rm ' + filename_all)
os.system('rm ' + 'log_cygnss_tle_' + tle_epoch + '.txt')
