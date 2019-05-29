# Compares netcdf prn/antenna selection as read in newer_validate_sift.py and the netcdf prn/antenna reported in Scott's CYGNSS_yaw_test_DMR_satselection.dat
## First you need to run (in the same terminal where you run newer_validate_sift.py: run -i cygnss_read_scott_prn_ant.py -> fill out variables nb_seconds_since_start_scott, gps_scott and which_ant_scott
# ASSUMPTIONS:
# - only work for idate = nb_date - 1

idate = nb_date - 1
# Firgure out when to start the comparison in Scott's file (ie when SpoCK and Scott start at the same date)
date_start_scott_str = '2018-09-26T00:00:01'
date_start_scott = datetime.strptime(date_start_scott_str, "%Y-%m-%dT%H:%M:%S")
date_when_spock_is_same_as_scott_start_date_temp = np.where(date_spock_same_time_as_netcdf == date_start_scott_str)[0]
iscott_start = 0
if len(date_when_spock_is_same_as_scott_start_date_temp) == 0: # SpOCK date starts later than Scott
    while ( ( date_start_scott + timedelta(seconds = iscott_start) ) != date_spock_same_time_as_netcdf[0]):
        iscott_start = iscott_start + 1

#iscott_start = 0 # !!!!!!!! remove this line
# Now compare gps_netcdf_all_date[idate]/gps_scott, which_ant_netcdf_all_date[idate]/which_ant_scott
## the netcdf file data is truncated from newer_validate_sift.py because tehre are times in the netdf file when there are masks. So filter Scott's data so it's taken at the same dates
inetcdf = 0
indices_scott_same_as_netcdf = []
scott_diff_netcdf = []
for iscott in range(iscott_start, nscott):
    print iscott, nscott
    date_scott_now = date_start_scott  + timedelta(seconds = nb_seconds_since_start_scott[iscott])
    if inetcdf < len(date_spock_same_time_as_netcdf):
        if (date_scott_now == date_spock_same_time_as_netcdf[inetcdf]):
            indices_scott_same_as_netcdf.append(iscott)
            if ((False in (gps_scott[iscott, :] == gps_netcdf_all_date[idate][inetcdf])) | (False in (which_ant_scott[iscott, :] == which_ant_netcdf_all_date[idate][inetcdf]))):
                scott_diff_netcdf.append(iscott)
                ipdb.set_trace()
            inetcdf = inetcdf + 1



    
