# This script reanmes every file and subfolder of a SpOCK simulaltion (it makes  a copy of the previous simulation before renaming, it doesn't just overwrite it)


import os


# PARAMETERS TO SET UP BEFORE RUNNNING THIS SCRIPT
root_old = '/Users/cbv/work/spockOut/beacon/spockMerged/out/newlong'#newfm02SepYaw_minus90_1s_out_which_ant_debug' # without last '/'
root_new = '/Users/cbv/work/spockOut/beacon/spockMerged/out/newlong_62degfilter' # without last '/'
# end of PARAMETERS TO SET UP BEFORE RUNNNING THIS SCRIPT

os.system('cp -Rp ' + root_old + ' ' + root_new)
dir_old = root_old.split('/')[-1]
dir_new = root_new.split('/')[-1]

sat_dir_list = next(os.walk(root_new))[1]
nsc = len(sat_dir_list)
if 'constellation_GPS' in sat_dir_list:
    new_sat_dir_list = []
    for isc in range(nsc):
        if sat_dir_list[isc] != 'constellation_GPS':
            new_sat_dir_list.append(sat_dir_list[isc])
    sat_dir_list = new_sat_dir_list
    nsc = len(sat_dir_list)

os.system('mv ' + root_new + '/' + 'CONSTELLATION_GPS_for_run_' + dir_old + '1.txt ' + root_new + '/' + 'CONSTELLATION_GPS_for_run_' + dir_new + '1.txt ')
os.system('mv ' + root_new + '/' + 'CONSTELLATION_CYGNSS_for_run_' + dir_old + '1.txt ' + root_new + '/' + 'CONSTELLATION_CYGNSS_for_run_' + dir_new + '1.txt ')
for isc in range(nsc):
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/LLA_' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/LLA_' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/ECEF_' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/ECEF_' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/TLE_' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/TLE_' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/attitude_' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/attitude_' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/density_' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/density_' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/given_output_' + sat_dir_list[isc] + '.txt' + ' ' + root_new + '/' + sat_dir_list[isc] + '/given_output_' + dir_new + str(isc+1) + '.txt')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + '/specular_' + sat_dir_list[isc] + '.bin' + ' ' + root_new + '/' + sat_dir_list[isc] + '/specular_' + dir_new + str(isc+1) + '.bin')
    os.system('mv ' + root_new + '/' + sat_dir_list[isc] + ' ' + root_new + '/' + dir_new + str(isc+1))
    