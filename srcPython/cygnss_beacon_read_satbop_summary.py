# This scripts reads the summary.txt output file created by sat-bop
# Inputs:
# - satbop_out_dir
# Outputs:
# - overp_time = list of overpass times (0 start, 1 end)

import numpy as np
import ipdb
# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
satbop_out_dir = '/Users/cbv/Google Drive/Work/PhD/Research/Code/cygnss/beacon/bruce/tds-bop V1.2.3/output/'
# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT

if satbop_out_dir[-1] != '/':
    satbop_out_dir = satbop_out_dir + '/'

summary_filename = satbop_out_dir + 'summary.txt'
summary_file = open(summary_filename)
read_summ = summary_file.readlines()
nline = len(read_summ)
iline = 0

# Find first overpass
## skip 'Current GPS TLE ages:' section. The firs overpasss is right after
str_interest = 'Current GPS TLE ages:'
while (str_interest in read_summ[iline]) == False:
    iline = iline + 1
while (len(read_summ[iline].replace('\r','').replace('\n','')) > 0 ):
    iline = iline + 1
iline = iline + 1

fm = []
prn = []
ipass = []
max_elev = []
dur = []
rise_time = []
set_time = []
while iline < nline:
    prn_sub = []
    fm_sub = read_summ[iline].split()[0]
    ipass_sub = (int)(read_summ[iline].split(':')[0].split()[-1])
    max_elev_sub = np.float(read_summ[iline+1].split(':')[1].split('deg')[0])
    dur_sub = np.float(read_summ[iline+4].split(':')[1].split('s')[0]) 
    rise_time_sub = read_summ[iline+5].split(':                ')[1].replace('\r','').replace('\n','')
    set_time_sub = read_summ[iline+9].split(':                ')[1].replace('\r','').replace('\n','')
    iline = iline + 12
    while (len(read_summ[iline].replace('\r','').replace('\n','')) > 0 ):
        prn_sub.append((int)(read_summ[iline].split('PRN')[1].split()[0]))
        iline = iline + 1
    iline = iline + 1

    fm.append(fm_sub)
    prn.append(prn_sub)
    ipass.append(ipass_sub)
    max_elev.append(max_elev_sub)
    dur.append(dur_sub)
    rise_time.append(rise_time_sub)
    set_time.append(set_time_sub)    
