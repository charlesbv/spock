import os
spock_filename_array = ['spec_spock_today_with_IIF.txt', 'spec_spock_2020_with_IIF.txt', 'spec_spock_2021_with_IIF.txt','spec_spock_2022_with_IIF.txt', 'spec_spock_2023_with_IIF.txt', 'spec_spock_2024_with_IIF.txt']
nrun = len(spock_filename_array)
for irun in range(nrun):
    #os.system("mpirun -np 4 spock_dev " + spock_filename_array[irun])
    #os.system("mpirun -np 4 spec " + spock_filename_array[irun] + " -lon=0 -rot=0 -min")
    os.system("python cygnss_coverage_5years.py " + spock_filename_array[irun])
                        


raise Exception
cell_width = 0.25#!!!!!!! should be 0.25 # in degrees
exclude_revisit_dt = cell_width*np.sqrt(2)*110/7.5 # if two visits are within exclude_revisit_dt seconds, then count the two visits as one visit
lat_max = 35
lat_min = -35
print_hour = 1. # print time every print_hour hour(s)
ncell_lon = (int)(360./cell_width) + 1
ncell_lat = (int)((lat_max - lat_min)/cell_width) + 1
ncell = ncell_lon*ncell_lat
cov_total = np.zeros([ncell_lon, ncell_lat])
time_visit = np.zeros([ncell_lon, ncell_lat, nsc*nspec]) -1 
date_ref = date_spock[0]
cov_time = np.zeros([ntime])
cell_covered = np.zeros([ncell_lon, ncell_lat])
ncell_covered_time = np.zeros([ntime])
itime = 0
while ((itime < ntime)):
    for isc in range(nsc):
        for ispec in range(nspec):
            icell_lon = (int)(lon_spec[itime, isc, ispec] / cell_width)
            if ((lat_spec[itime, isc, ispec] >= lat_min) & (lat_spec[itime, isc, ispec] <= lat_max)): # theya re SPs at latitudes > 35 deg
                if gain_spec[itime, isc, ispec] >= 4: #!!!!!!! sure you want to keep the condition "if gain_spec[itime, isc, ispec]	> 0"
                    icell_lat = (int)((lat_spec[itime, isc, ispec] - lat_min) / cell_width)
                    time_visit[icell_lon, icell_lat, isc*nspec + ispec] = (date_spock[itime] - date_ref).total_seconds()
                    cell_covered[icell_lon, icell_lat] = cell_covered[icell_lon, icell_lat] + 1
    ncell_covered_time[itime] = len(np.where(cell_covered != 0)[0])
    cov_time[itime] = ncell_covered_time[itime] * 100./ncell
    if cov_time[itime]  >= 70.:
        break
    if np.mod(itime, print_hour*3600) == 0:
        print itime, ntime, str(itime/3600) + 'h', format(cov_time[itime], ".0f") + '%'    
    itime = itime + 1
revisit_time = []

for icell_lon in range(ncell_lon):
    #print icell_lon, ncell_lon
    for icell_lat in range(ncell_lat):
        visiter = np.where(time_visit[icell_lon, icell_lat, :] != -1)[0]
        if len(visiter) > 0:
            time_visit_temp_not_sorted = time_visit[icell_lon, icell_lat, visiter]
            sort_time_visit = np.argsort(time_visit_temp_not_sorted)
            time_visit_temp = time_visit_temp_not_sorted[sort_time_visit]
            nvisit = len(time_visit_temp)
            time_visit_valid = []
            time_visit_valid.append(time_visit_temp[0])
            for ivisit in range(nvisit-1):
                if time_visit_temp[ivisit+1] - time_visit_temp[ivisit] > exclude_revisit_dt: #  if two visits are within exclude_revisit_dt seconds, then count the two visits as one visit 
                    time_visit_valid.append(time_visit_temp[ivisit+1])
                    revisit_time.append( time_visit_valid[-1] - time_visit_valid[-2] )
            if len(time_visit_valid) > 0: # this cell has been visited at least once
                cov_total[icell_lon, icell_lat] = cov_total[icell_lon, icell_lat] + 1
revisit_time = np.array(revisit_time)
ilon_cell_not_covered = np.where(cov_total == 0)[0]
ilat_cell_not_covered = np.where(cov_total == 0)[0]
ncell_not_covered = len(ilon_cell_not_covered)
ncell_covered = ncell - ncell_not_covered
perc_cov = ncell_covered * 100./ncell
mean_revisit_time = np.mean(revisit_time) / 3600. # in hours

print 'Percentage coverage: ' + format(perc_cov, ".1f")
print 'Mean revisit time: ' + format(mean_revisit_time, ".1f")
