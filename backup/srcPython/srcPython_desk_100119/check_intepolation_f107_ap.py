from get_prop_dir import *
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta

plt.ion()
# Note: if ap historical or if it comes from omniweb then this script won't compute Ap
source_file = "omniweb" # external 

# FILES REFERENCE
if source_file == "external":
    ## F10.7
    filename_f107_ref = get_prop_dir(1) + "run_aerie/input/density/density_NRLMSIS00e/modif_f107_high_not_ok.txt"
    file_f107_ref = open(filename_f107_ref)
    read_file_f107_ref = file_f107_ref.readlines()
    n_header = 3
    n = len(read_file_f107_ref) - n_header - 1
    f107_ref = np.zeros([n])
    for istep in range(n):
        f107_ref[istep] = np.float( read_file_f107_ref[istep + n_header].split()[3] )

    ## Ap
    filename_ap_ref = get_prop_dir(1) + "run_aerie/input/density/density_NRLMSIS00e/modif_ap_high_not_ok.txt"
    file_ap_ref = open(filename_ap_ref)
    read_file_ap_ref = file_ap_ref.readlines()
    n_header = 3
    n = len(read_file_ap_ref) - n_header - 1
    ap_ref = np.zeros([n])
    for istep in range(n):
        ap_ref[istep] = np.float( read_file_ap_ref[istep + n_header].split()[3] )

if source_file == "omniweb":
    ## F10.7
    filename_f107_ref = get_prop_dir(1) + "run_aerie/input/density/density_NRLMSIS00e/f107_for_81daverage_20130621_to_20140209_not_ok.txt"
    file_f107_ref = open(filename_f107_ref)
    read_file_f107_ref = file_f107_ref.readlines()
    n_header = 8
    n_end = 15
    n = len(read_file_f107_ref) - n_header - n_end
    f107_ref = np.zeros([n])
    for istep in range(n):
        f107_ref[istep] = np.float( read_file_f107_ref[istep + n_header].split()[3] )


# FILE OUTPUT
filename_output = get_prop_dir(1) + "run_aerie/output/aerie_omniweb_test1/aerie_omniweb_test11/given_output_aerie_omniweb_test11.txt"
file_output = open(filename_output)
read_file_output = file_output.readlines()
n_out = len(read_file_output)
f107_out = np.zeros([n_out])
f107A_out = np.zeros([n_out])
ap_out = np.zeros([n_out])
date_elapsed = np.zeros([n_out])
for istep in range(n_out):
    var_for_time = read_file_output[istep].split()[0] + 'T' + read_file_output[istep].split()[2][0:-5] 
    if istep == 0:
        date_out_start = datetime.strptime( var_for_time, "%Y-%jT%H:%M:%S")
    f107_out[istep] = np.float( read_file_output[istep].split()[3] )    
    f107A_out[istep] = np.float( read_file_output[istep].split()[4] )    
    ap_out[istep] = np.float( read_file_output[istep].split()[5] )   # see note at top of script
    time_elapsed = (datetime.strptime( var_for_time, "%Y-%jT%H:%M:%S") - date_out_start)
    date_elapsed[istep] = time_elapsed.days * 24 * 3600 + time_elapsed.seconds

# COMPARE WITH PLOT
## F10.7 and F10.7A
fig, ax = plt.subplots()
if source_file == "external": # daily
    x_ref = np.arange(0, n)*24*3600
if source_file == "omniweb": # hourly
    x_ref = np.arange(0, n)*3600

ax.plot(date_elapsed, f107_out, linewidth = 2, color = 'b', linestyle = '', marker = '.')
ax.plot(date_elapsed, f107A_out, linewidth = 2, color = 'k', linestyle = '', marker = '.')
ax.plot(x_ref, f107_ref, linewidth = 2, color = 'r',marker = 'o', linestyle = '', markersize = 10)
ax.margins(0,0)
plt.show()

if source_file == "omniweb":
    raise Exception

## Ap daily
fig, ax = plt.subplots()
x_ref = np.arange(0, n)*24*3600
ax.plot(date_elapsed, ap_out, linewidth = 2, color = 'b', linestyle = '', marker = '.')
ax.plot(x_ref, ap_ref, linewidth = 2, color = 'r',marker = 'o', linestyle = '', markersize = 10)
ax.margins(0,0)
plt.show()
