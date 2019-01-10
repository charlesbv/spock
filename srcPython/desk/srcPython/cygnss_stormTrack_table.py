# This script creates the table of stormss for the webpage Hurricane Tracking -> wind measurements

import os
from os import listdir
from os.path import isfile, join
import glob, os
import ipdb
import gzip
import numpy as np
list_storm = [f for f in listdir('data') if isfile(join('data', f))]
n = len(list_storm)

filename_out = 'stormTrack.txt'
file_out = open(filename_out, "w+")
name_all = []
id_all = []
start_all = []
stop_all = []
ani_all = []
for i in range(n):
    filename = list_storm[i]
    forecastLines = []
    with gzip.open('data/' + filename,'r') as fin:        
        for line in fin:        
            forecastLines.append(line)
    nline = len(forecastLines)
    iline = 0
    cast = [x.strip() for x in forecastLines[iline].split(',')]
    idstorm = cast[0] + cast[1]
    start_date_temp = cast[2]
    start_date = start_date_temp[:4] + '/' + start_date_temp[4:6] + '/' +  start_date_temp[6:8] + 'T' + start_date_temp[8:10] + ':00'
    year = start_date_temp[:4]
    iline = nline - 1
    cast = [x.strip() for x in forecastLines[iline].split(',')]
    stop_date_temp = cast[2]
    stop_date = stop_date_temp[:4] + '/' + stop_date_temp[4:6] + '/' +  stop_date_temp[6:8] + 'T' + stop_date_temp[8:10] + ':00'
    name_storm = cast[27].title()
 
    ani = idstorm + year + '.mp4'# AL112018
    name_all.append(name_storm)
    id_all.append(idstorm)
    start_all.append(start_date)
    stop_all.append(stop_date)
    ani_all.append(ani)

index_sort = sorted(range(len(name_all)), key=lambda k: name_all[k])
name_all_sorted = np.array(name_all)[index_sort]
id_all_sorted = np.array(id_all)[index_sort]
start_all_sorted = np.array(start_all)[index_sort]
stop_all_sorted = np.array(stop_all)[index_sort]
ani_all_sorted = np.array(ani_all)[index_sort]
for i in range(n):
    print >> file_out, name_all_sorted[i] + ' ' + id_all_sorted[i] + ' ' + start_all_sorted[i] + ' ' + stop_all_sorted[i] + ' ' + ani_all_sorted[i]

file_out.close()
