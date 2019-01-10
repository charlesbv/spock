# This scirpt is a copy of cygnss_coverage_temporal_spatial_relation_equal_cell_area2 on Sep 16 2018. It was made to answer Chris Ruf's analysis questions by email on Sep 13 2018.                                                          

pleiades = 0 # set this to 1 if running this script on pleiades
import os
import sys 
if pleiades == 1:
    sys.path.append("/home1/cbussy/Code/spock/srcPython")
else:
    sys.path.append("/Users/cbv/Google Drive/Work/PhD/Research/Code/spock/srcPython") 
    #sys.path.append("/home/cbv/spock_development_new_structure_kalman_dev/srcPython")
    import matplotlib
    matplotlib.use("Agg") # without this line, when running this script from a cronjob we get an error "Unable to access the X Display, is $DISPLAY set properly?"
    from matplotlib import pyplot as plt
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.basemap import Basemap, shiftgrid
    from collections import *
    from matplotlib import colors
    import matplotlib.ticker as ticker
    import ipdb
from read_input_file import *
from find_in_read_input_order_variables import *
from read_spec_spock import *
from cygnss_read_spock_spec_bin import *
import pickle
re = 6371.0 # mean Earth radius 
input_filename = sys.argv[1]
# PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT
#load = 0 # set this to 1 to load a previsouly pickle of running all this script 
max_lat_grid_array = [35]#latitude of the top of the grid. There's a good reason to take the top and not the bottom (when we calculate the number of grids in the longitudinal direction, it depends of the latitude. To avoid overlap after one round of the Earth, we need to take the top of the grid, not the bottom). See more explnaation further in the code
# for cyggnss extended: lat_grid_center_array =  [5, 10, 15, 20, 25, 30]         
# for inclinatino at 90 deg: lat_grid_center_array = [60,65,70,75,80,85]           
pickle_root = input_filename.replace(".txt","")
max_lat_grid_array = np.array(max_lat_grid_array) * np.pi / 180.
ilat_center = 0
rlat_width = re * max_lat_grid_array[ilat_center]#500 # latitude width of the grid, in km
rlon_width = 500.#2 * np.pi * re * np.cos(max_lat_grid_array[ilat_center])#500 # longitude width of the grid, in km
pickle_folder = "/Users/cbv/coverage_temporal_spatial_relation/oct18" # with or without the last '/' (doesnt matter). also the directory doesnt need to exist 
# pickle_folder = "/nobackup/cbussy/coverage_temporal_spatial_relation/sep18/pickle" # with or without the last '/' (doesnt matter). also the directory doesnt need to exist 
pickle_visu_subfolder = "2d_visu" # with or without the last '/' (doesnt matter). also the directory doesnt need to exist 
# end of PARAMETERS TO SET UP BEFORE RUNNING THIS SCRIPT

if os.path.isdir(pickle_folder) == False:
    os.system("mkdir " + pickle_folder)
if pickle_folder[-1] != '/':
    pickle_folder = pickle_folder + '/'

if os.path.isdir(pickle_folder + pickle_visu_subfolder) == False:
    os.system("mkdir " + pickle_folder + pickle_visu_subfolder)
if pickle_visu_subfolder[-1] != '/':
    pickle_visu_subfolder = pickle_visu_subfolder + '/'




width_cell_array = [20.] # in km
coverage_array = [70.] # in percentage in ASCENDING ORDER


color_array = ['k','cornflowerblue','r','g', 'm', 'gold', 'cyan', 'fuchsia', 'lawngreen', 'darkgray', 'green', 'chocolate']


nb_type_cell = len(width_cell_array)
nb_coverage = len(coverage_array)

width_cell_array = np.array(width_cell_array)
coverage_array = np.array(coverage_array)
    

#if load != 1:
    # Read the ECEF position of the specular points. Note: if less tha 4 spec for a time then ecef components are equal to 99999999


var_in, var_in_order = read_input_file(input_filename)
dt_spock_output = var_in[find_in_read_input_order_variables(var_in_order, 'dt_output')]; 
output_file_path_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_path_list')]; 
output_file_name_list = var_in[find_in_read_input_order_variables(var_in_order, 'output_file_name_list')];
gps_name_list_spock = var_in[find_in_read_input_order_variables(var_in_order, 'gps_name')];
nb_sc = var_in[find_in_read_input_order_variables(var_in_order, 'nb_sc')];
cygfm_to_spock_nb = [4,3,8,2,1,6,7,5] # ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']


lon_spec_list = []; lat_spec_list = [];

nb_spec = 4 # if less than 4 sp (eg: removing block IIF GPS, happens 0.05% of the time to have less than 4 SPs) then the lon, lat and gain will be set to 99999999

print 'Reading SP positions...'
for isc in range(nb_sc):
    print 'isc', isc
    spec_spock_filename = output_file_path_list[isc] + "specular_" + output_file_name_list[isc]
    data_spec = cygnss_read_spock_spec_bin(spec_spock_filename.replace('.txt','.bin'), gps_name_list_spock, dt_spock_output, 0) 
    date_spock = data_spec[0]; lon_spec_sc = data_spec[1]; lat_spec_sc = data_spec[2];   # !!!!!!! to take into account all gains: lon_spec_sc = data_spec[1]; lat_spec_sc = data_spec[2]; lon_spec_sc = data_spec[5]; lat_spec_sc = data_spec[6];  to take into account only SPs with non-0 gain
    gain_spec_sc = data_spec[3]; gps_name_sc = data_spec[4]; 
    if isc == 0:
        nb_time = len(lon_spec_sc)
        lon_spec = np.zeros([nb_spec, nb_sc, nb_time]) + 99999999*180/np.pi
        lat_spec = np.zeros([nb_spec, nb_sc, nb_time]) + 99999999*180/np.pi
    for itime in range(nb_time):
        for ispec in range(len(lon_spec_sc[itime])):
            lon_spec[ispec, isc, itime] = lon_spec_sc[itime][ispec]
            lat_spec[ispec, isc, itime] = lat_spec_sc[itime][ispec]

print 'Done reading the SP positions.'
# lon_spec_temp = np.array(lon_spec_list)
# lat_spec_temp = np.array(lat_spec_list)
# lon_spec = np.moveaxis(lon_spec_temp,[0,1,2],[1,2,0]) # reoganization of shape to be [nb_spec, nb_sc, nb_time]
# lat_spec = np.moveaxis(lat_spec_temp,[0,1,2],[1,2,0])
lon_spec = lon_spec * np.pi / 180.
lat_spec = lat_spec * np.pi / 180.
print ""
#print lon_spec

##################### =
# Sel<ecting only the specular ponts in the window (also called grid) defined by lon/min_lat_grid and lon/rlat_width, count the number of specular points in each cell (defined by width_cell) as a function of time
nb_lat_center = len(max_lat_grid_array) # contratily to nb_lon_center, here we move the grid with respect to a speciafic array chosen by the user.
#if load != 1:




#nb_cell_covered = np.zeros([nb_lat_center, nb_lon_center, nb_type_cell,  nb_time])
cell_array = []
nb_cell_covered_average_over_lon = np.zeros([nb_lat_center, nb_type_cell, nb_time])
percentage_cell_at_least_one_revisit_average_over_lon = np.zeros([nb_lat_center, nb_type_cell, nb_time])
#time_to_reach_coverage_average_over_lon = np.zeros([nb_lat_center, nb_type_cell, nb_coverage])
dt_visu_hour = 1 # time step for the visu in hour
dt_visu_spec_hour = 1 # how long the track of spec stay on the visu
dt_visu = dt_visu_hour * 3600. # !!!need to be a double here
dt_visu_spec = dt_visu_spec_hour * 3600
nb_dt_visu = (int) (nb_time / dt_visu)
nb_day_visu_with_spec_plotted = 2# beyond this limit, the path of spec wont be pltted on the animation, just the grid cells                                 
nb_second_visu = nb_day_visu_with_spec_plotted * 24*3600
save_lon_spec = lon_spec[:,:,:nb_second_visu] # save only the first nb_day_visu day of spec
save_lat_spec = lat_spec[:,:,:nb_second_visu]
mean_revisit_time_average_over_lon = np.zeros([nb_lat_center, nb_type_cell])
for itype in range( nb_type_cell): # !!!!!!!!!!
    #print itype, nb_type_cell-1
    width_cell = width_cell_array[itype]
    no_revisit_dt = 1.25#95*width_cell*np.sqrt(2)/(2*np.pi*re)*60 # this is the max time in seconds that the subsatellite point stays in a cell!!!!!!!! assumes the sc orbital period is 95 min
    nb_cell_lat = (int) (rlat_width / width_cell) # nb cell bins in one grid in the latitudinal direction. always the same, whatever latitude 
    nb_cell_lon = (int) (rlon_width / width_cell) # nb cell bins in one grid in the longitudinal direction. always the same, whatever latitude
    nb_cell_in_grid = nb_cell_lon * nb_cell_lat
    for ilat_center in range(nb_lat_center):# !!!!!!!!!!!!!
        icell_max_lat_grid = (int) ( re * max_lat_grid_array[ilat_center]  / width_cell ) # cell# that corresponds to the top of the grid
        icell_min_lat_grid = icell_max_lat_grid - nb_cell_lat  # cell# that corresponds to the bottom of the grid
        nb_lon_center = (int) (2 * np.pi * re * np.cos(max_lat_grid_array[ilat_center]) / rlon_width) # since the cells are equal area, the number of grids in the longitudinal direction depends on the latitude. That's why here we take the top and the bottom because if we took the bottom, then the top could overlap after one round around the Earth
        #time_to_reach_coverage = np.zeros([ nb_lon_center, nb_coverage]) - 1
        nb_cell_covered = np.zeros([ nb_lon_center, nb_time]) - 1
        mean_revisit_time = np.zeros([nb_lon_center])
        percentage_cell_at_least_one_revisit = np.zeros([nb_lon_center, nb_time]) - 1
        time_since_last_visit_of_this_cell_list = []
        for ilon_center in range(nb_lon_center): #!!!!!!!!!!!
            time_since_last_visit_of_this_cell_list_lon = []
            icell_min_lon_grid = ilon_center*nb_cell_lon
            icell_max_lon_grid = ilon_center*nb_cell_lon + nb_cell_lon
            rlon_spec = re * np.cos(lat_spec) * lon_spec # distance from Greenwhich meridian to spec in the longitudinal direction
            rlat_spec = re * lat_spec # distance from Greenwhich meridian to spec in the latitudinal direction
            icell_lon = (rlon_spec / width_cell).astype(int) # cell# corresponding toe the spec in the longitudinal direction
            icell_lat =  (rlat_spec / width_cell).astype(int) # cell# corresponding toe the spec in the latitudinal direction
            where_in_grid = np.where( ( icell_lon >= icell_min_lon_grid ) & ( icell_lon < icell_max_lon_grid ) & ( icell_lat >= icell_min_lat_grid ) & ( icell_lat < icell_max_lat_grid  ) & ( lat_spec < 9999999 ) )
            indices_time_sorted = [i[0] for i in sorted(enumerate(where_in_grid[2]), key=lambda x:x[1])]
            which_spec = where_in_grid[0][indices_time_sorted]#[x for _,x in sorted(zip(where_in_grid[2],where_in_grid[0]))]#where_in_grid[0]
            which_sc = where_in_grid[1][indices_time_sorted]#[x for _,x in sorted(zip(where_in_grid[2],where_in_grid[1]))]#where_in_grid[1]
            which_time = where_in_grid[2][indices_time_sorted] #[x for _,x in sorted(zip(where_in_grid[2],where_in_grid[2]))]#where_in_grid[2]
            nb_time_spec_in_grid = len(which_time)

            icoverage = 0
            print itype, nb_type_cell-1, ilat_center, nb_lat_center-1, ilon_center,nb_lon_center-1,nb_time_spec_in_grid
            icell_save = []
            cell_ilon_ilat_itype = np.zeros([nb_cell_lat, nb_cell_lon])
            nb_cell_filled = np.zeros([nb_time])
            nb_cell_at_least_one_revisit = np.zeros([nb_time])
            if ilon_center == 0: 
                if ((itype == 0) | (itype == nb_type_cell-1)):
                    cell_array_type_lon_lat_cumul = np.zeros([nb_dt_visu, nb_cell_lat, nb_cell_lon])# # 2d visu only for one position of grid in longitudinal direction, for the smallest type cell and largest type cell, for all lat center
            #itime_visu_previous = 0
            cell_ilon_ilat_itype_save_time = np.zeros([nb_cell_lat, nb_cell_lon]) # for each cell, record the time when a spec was in.  This is to calculate the revisit time. 
            revisit_time = 0 # the revisit time is the time difference between all measurements in the bin.
            nb_revisit_time = 0 # total number of time all cells have been revisited
            first_time_above_70_revisit = 0
            for itime in range( nb_time_spec_in_grid ):
                #print itime, nb_time_spec_in_grid
                if itime > 0:
                    nb_cell_filled[which_time[itime-1]+1:which_time[itime]+1] = nb_cell_filled[which_time[itime-1]]
                    nb_cell_at_least_one_revisit[which_time[itime-1]+1:which_time[itime]+1] = nb_cell_at_least_one_revisit[which_time[itime-1]]
                    #                         if ilon_center == 0: # 2d visu only for one position of grid in longitudinal direction
                    #                             if ((itype == 0) | (itype == nb_type_cell-1)):
                    #                                 if which_time[itime]+1 < nb_second_visu-1:
                    #                                     cell_array_type_lon_lat_cumul[which_time[itime-1]+1:which_time[itime]+1, :, :] = cell_array_type_lon_lat_cumul[which_time[itime-1], :, :]


                icell_lat_here = icell_lat[which_spec[itime], which_sc[itime], which_time[itime]]
                icell_lon_here = icell_lon[which_spec[itime], which_sc[itime], which_time[itime]]
                if cell_ilon_ilat_itype[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] == 0: # if this cell has not been filled yet. We remove icell_min_lat_grid and icell_min_lon_grid because cell_ilon_ilat_itype is defiend from 0 to nb_cell_lat and 0 to nb_cell_lon
                    nb_cell_filled[which_time[itime]] = nb_cell_filled[which_time[itime-1]] + 1
                    cell_ilon_ilat_itype[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] = cell_ilon_ilat_itype[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] + 1
                else: # this is not the first time is cell has been visited 
                    time_since_last_visit_of_this_cell = which_time[itime] - cell_ilon_ilat_itype_save_time[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid]
                    if (((nb_cell_at_least_one_revisit[which_time[itime-1]] * 100. / nb_cell_in_grid) >= 70.) & (first_time_above_70_revisit == 0)) :
                        print  which_time[itime] / 3600.
                        first_time_above_70_revisit = 1
                    if ((time_since_last_visit_of_this_cell > no_revisit_dt) & ((nb_cell_at_least_one_revisit[which_time[itime-1]] * 100. / nb_cell_in_grid) < 70. )): # only record the revisit time if it's larger than the time it takes for a sc to fly over the cell and don't go revisits once the coverage of the grid reached 70% (that's how Dorina's algo is defined
                        revisit_time = revisit_time + time_since_last_visit_of_this_cell
                        nb_revisit_time = nb_revisit_time + 1
                        time_since_last_visit_of_this_cell_list_lon.append(time_since_last_visit_of_this_cell)
                        if cell_ilon_ilat_itype[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] == 1: # if this cell has only been visited once so far
                            nb_cell_at_least_one_revisit[which_time[itime]] = nb_cell_at_least_one_revisit[which_time[itime-1]] + 1
                            cell_ilon_ilat_itype[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] = cell_ilon_ilat_itype[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] + 1
                        

                cell_ilon_ilat_itype_save_time[icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] = which_time[itime]
                #                     raise Exception
                #                     print which_time[itime]/dt_visu, (int)(which_time[itime]/dt_visu)


                if ilon_center == 0: 
                    if ((itype == 0) | (itype == nb_type_cell-1)):
                        itime_visu = (int)(which_time[itime]/dt_visu) 
                        if which_time[itime]/dt_visu == itime_visu: # this time falls exactly at a time of the animation so take a "screenshot" of the grid
                            cell_array_type_lon_lat_cumul[itime_visu, :, :] = cell_ilon_ilat_itype 
                            #                         if itime+1 < len(which_time):
                            #                             if  which_time[itime+1]/dt_visu != itime_visu: #move to next visu dt only if no more spec in grid at this time ( which_time[itime+1] =  which_time[itime] if more than one spec at this time that falls in the grid)
                            #                                 itime_visu = itime_visu + 1

                        else: # this time does not fall at a time of the animation...
                            if itime+1 < len(which_time):
                                if which_time[itime+1]/dt_visu > itime_visu+1:  #... but the following time jumps further in the future than the time step of the animation so all time steps in the animation until this time should show the grid as it is at this time
                                    cell_array_type_lon_lat_cumul[itime_visu+1:(int)(which_time[itime+1]/dt_visu)+1, :, :] = cell_ilon_ilat_itype
    #                                itime_visu = (int)(which_time[itime+1]/dt_visu)+1

                                    #                     if np.mod(which_time[itime], dt_visu) == 0: # for the animation, record which cell is filled only every dt_visu 
                                    #                         if ilon_center == 0:  # 2d visu only for one position of grid in longitudinal direction
                                    #                             if ((itype == 0) | (itype == nb_type_cell-1)):
                                    #                                 if ((int) (which_time[itime] / dt_visu)) < nb_dt_visu:
                                    #                                     itime_visu = (int) (which_time[itime] / dt_visu)
                                    #                                     print itime_visu_previous, itime_visu,which_time[itime]/3600.
                                    #                                     cell_array_type_lon_lat_cumul[itime_visu_previous:itime_visu, :, :] = cell_array_type_lon_lat_cumul[itime_visu_previous, :, :]
                                    #                                     cell_array_type_lon_lat_cumul[itime_visu, :, :] = cell_ilon_ilat_itype
                                    #                                     itime_visu_previous = itime_visu
                                    #                                        cell_array_type_lon_lat_cumul[which_time[itime], icell_lat_here-icell_min_lat_grid, icell_lon_here-icell_min_lon_grid] =  1 #We remove icell_min_lat_grid and icell_min_lon_grid becausecell_array_type_lon_lat_cumul is defiend from 0 to nb_cell_lat and 0 to nb_cell_lon
                                #itime_visu_save = itime
                                #                         if icoverage < nb_coverage:
                                #                             if nb_cell_filled[which_time[itime]] >= coverage_array[icoverage]/100. * nb_cell_in_grid:
                                #                                 time_to_reach_coverage[ilon_center, icoverage] = which_time[itime] # since the specular file is every second, itime is moving one second by one second
                                #                                 #time_to_reach_coveragze[ilat_center, ilon_center, itype, icoverage] = which_time[itime] # since the specular file is every second, itime is moving one second by one second
                                #                                 icoverage = icoverage + 1
            time_since_last_visit_of_this_cell_list.append(time_since_last_visit_of_this_cell_list_lon)
            pickle.dump(time_since_last_visit_of_this_cell_list, open( pickle_folder +  pickle_root +  "_time_since_last_visit_of_this_cell_list_ilat_" + str((int)(max_lat_grid_array[ilat_center]*180./np.pi)) + '_cell_' + str(width_cell_array[itype])   + "_gain_not0.pickle", "w"))
            if nb_time_spec_in_grid > 0 :
                if nb_revisit_time > 0:
                    mean_revisit_time[ilon_center] = revisit_time / nb_revisit_time
                nb_cell_filled[which_time[itime]+1:] = nb_cell_filled[which_time[itime]]
                nb_cell_at_least_one_revisit[which_time[itime]+1:] = nb_cell_at_least_one_revisit[which_time[itime]]
                nb_cell_covered[ ilon_center, :] = nb_cell_filled * 100. / nb_cell_in_grid                    
                percentage_cell_at_least_one_revisit[ ilon_center, :] = nb_cell_at_least_one_revisit * 100. / nb_cell_in_grid                    
                if ilon_center == 0:  # 2d visu only for one position of grid in longitudinal direction  
                    if ((itype == 0) | (itype == nb_type_cell-1)):
                        #                            cell_array_type_lon_lat_cumul[which_time[itime_visu_save]+1: ,:, :] = cell_array_type_lon_lat_cumul[which_time[itime_visu_save], :, :]
                        pickle.dump(cell_array_type_lon_lat_cumul, open( pickle_folder + pickle_visu_subfolder +  pickle_root +  "_cell_array_type_lon_lat_cumul_ilat_" + str((int)(max_lat_grid_array[ilat_center]*180./np.pi)) + '_cell_' + str(width_cell_array[itype])   + "_gain_not0.pickle", "w"))
                        #             for icoverage in range(nb_coverage):
                        #                 where_lon_reached_coverage = np.where( time_to_reach_coverage[:,icoverage] != 1 )[0] #time_to_reach_coverage[ilon_center, icoverage] = -1 if coverage_array[icoverage] has not been reached at this lon ilon_center
                        #                 time_to_reach_coverage_average_over_lon[ilat_center, itype, icoverage] = np.mean( time_to_reach_coverage[where_lon_reached_coverage, icoverage])
        nb_cell_covered_average_over_lon[ilat_center, itype, :] = np.mean( nb_cell_covered[ :, :], axis = 0 )
        percentage_cell_at_least_one_revisit_average_over_lon[ilat_center, itype, :] = np.mean( percentage_cell_at_least_one_revisit[ :, :], axis = 0 )
        mean_revisit_time_average_over_lon[ilat_center, itype] = np.mean(mean_revisit_time)

pickle.dump(nb_cell_covered_average_over_lon, open( pickle_folder + pickle_root + "_nb_cell_covered_average_over_lon_gain_not0.pickle", "w"))
pickle.dump(percentage_cell_at_least_one_revisit_average_over_lon, open( pickle_folder + pickle_root + "_percentage_cell_at_least_one_revisit_average_over_lon_gain_not0.pickle", "w"))
pickle.dump(mean_revisit_time_average_over_lon, open( pickle_folder + pickle_root + "_mean_revisit_time_average_over_lon_gain_not0.pickle", "w"))
#    pickle.dump(time_to_reach_coverage_average_over_lon, open(  pickle_folder + pickle_root + "_time_to_reach_coverage_average_over_lon_gain_not0.pickle", "w"))
pickle.dump(save_lon_spec, open( pickle_folder + pickle_visu_subfolder + pickle_root + "_save_lon_spec_gain_not0.pickle", "w"))
pickle.dump(save_lat_spec, open( pickle_folder + pickle_visu_subfolder + pickle_root + "_save_lat_spec_gain_not0.pickle", "w"))

raise Exception

#else:
    #time_to_reach_coverage_average_over_lon = pickle.load( open( pickle_folder + pickle_root + "_time_to_reach_coverage_average_over_lon.pickle", "r"))
#    nb_cell_covered_average_over_lon = pickle.load( open( pickle_folder + pickle_root + "_nb_cell_covered_average_over_lon.pickle", "r"))




if pleiades != 1:
    ## Parameters for the figure
    height_fig = 11.  # the width is calculated as height_fig * 4/3.
    fontsize_plot = 20 
    ratio_fig_size = 4./3

    raise Exception

    # 2D visualization
    itype = 2#3 # grid bin type
    width_cell = width_cell_array[itype]
    max_lat_grid_array = np.array(max_lat_grid_array)
    ilat_center = np.where(max_lat_grid_array == 30)[0][0]#2
    ilon_center = 9
    min_lat_grid = max_lat_grid_array[ilat_center]
    lon_grid_center = lon_grid_center_array[ilon_center]
    lat_min_grid = min_lat_grid - lat_width / 2.
    lat_max_grid = min_lat_grid + lat_width / 2.
    nb_cell_lat = (int) ( ( lat_max_grid - lat_min_grid ) / width_cell )
    lat_max_grid_cell = lat_min_grid + nb_cell_lat * width_cell # latitude of the cell of the grid that correspond to the highest latitude. 
    lon_min_grid = (lon_grid_center - lon_width / 2.)%360
    lon_max_grid = (lon_grid_center + lon_width / 2. )%360
    nb_cell_lon = (int) ( ( lon_max_grid - lon_min_grid ) / width_cell ) 
    lon_max_grid_cell = lon_min_grid + nb_cell_lon * width_cell # longitude of the cell of the grid that correspond to the highest longitude.
    nb_cell_in_grid = nb_cell_lon * nb_cell_lat
    cell_array_type_lon_lat_cumul = np.zeros([nb_time, nb_cell_lat, nb_cell_lon])# +1 for safety
    cell_array_type_lon_lat = np.zeros([nb_time, nb_cell_lat, nb_cell_lon])# +1 for safety


    where_in_grid = np.where(( lon_spec[:, :, :] >= lon_min_grid ) & ( lon_spec[:, :, :] < lon_max_grid_cell ) & ( lat_spec[:, :, :] >= lat_min_grid ) & ( lat_spec[:, :, :] < lat_max_grid_cell  ) & ( lat_spec[:, :, :] < 9999999 ) )
    indices_time_sorted = [i[0] for i in sorted(enumerate(where_in_grid[2]), key=lambda x:x[1])]
    which_spec = where_in_grid[0][indices_time_sorted]#[x for _,x in sorted(zip(where_in_grid[2],where_in_grid[0]))]#where_in_grid[0]
    which_sc = where_in_grid[1][indices_time_sorted]#[x for _,x in sorted(zip(where_in_grid[2],where_in_grid[1]))]#where_in_grid[1]
    which_time = where_in_grid[2][indices_time_sorted] #[x for _,x in sorted(zip(where_in_grid[2],where_in_grid[2]))]#where_in_grid[2]
    nb_time_spec_in_grid = len(which_time)
    cell_ilon_ilat_itype = np.zeros([nb_cell_lat, nb_cell_lon])
    nb_cell_filled = np.zeros([nb_time])

    for itime in range( nb_time_spec_in_grid ):
        print itime, nb_time_spec_in_grid
        if itime > 0:
            cell_array_type_lon_lat_cumul[which_time[itime-1]+1:which_time[itime]+1, :, :] = cell_array_type_lon_lat_cumul[which_time[itime-1], :, :]

        lon_spec_in_grid = lon_spec[which_spec[itime], which_sc[itime], which_time[itime]]  
        lat_spec_in_grid = lat_spec[which_spec[itime], which_sc[itime], which_time[itime]]
        icell_lat = (int) ( ( lat_spec_in_grid - lat_min_grid ) / width_cell )
        icell_lon = (int) ( ( lon_spec_in_grid - lon_min_grid ) / width_cell )


        if cell_ilon_ilat_itype[icell_lat, icell_lon] == 0: # if this cell has not been filled yet
            cell_ilon_ilat_itype[icell_lat, icell_lon] = 1

            cell_array_type_lon_lat_cumul[which_time[itime], icell_lat, icell_lon] =  1


    cell_array_type_lon_lat_cumul[which_time[itime]+1:] = cell_array_type_lon_lat_cumul[which_time[itime]]

    y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
    x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'
    cmap = colors.ListedColormap(['azure', 'lightskyblue'])

    #itime = 9 * 3600-1
    itime_count = -1
    for itime in np.arange( 0,  3 * 24 *3600 , 3600):
    #for itime in np.arange( 0, 2, 1):
        itime_count = itime_count + 1
        print itime,  12 * 3600
        fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
        fig.suptitle('', y = 0.96,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold
        gs = gridspec.GridSpec(1, 1)
        gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
        ax_2d = fig.add_subplot(gs[0, 0])


        nb_days_from_start = itime / 3600 / 24
        time_since_start_day = (int) ( np.floor(itime / 3600. / 24) )
        time_since_start_hour =(int) (  np.floor(  (itime / 3600. / 24 - time_since_start_day )*24) )
        time_since_start_min = (int) ( np.floor(( (itime / 3600. / 24 - time_since_start_day )*24 - time_since_start_hour ) * 60) )
        time_since_start_sec = (int) ( np.round( ( ( (itime / 3600. / 24 - time_since_start_day )*24 - time_since_start_hour ) * 60 - time_since_start_min ) * 60 ) )
        time_since_start = str(time_since_start_day) + 'd ' + str(time_since_start_hour) + 'h' + str(time_since_start_min) + "'"  + str(time_since_start_sec) + "''"
    #    percentage_covered = nb_cell_covered[ilat_center, ilon_center, itype, itime]

        ax_2d_title = '2D visualization of coverage (lon = ' + str(lon_grid_center) + u'\N{DEGREE SIGN}'+ ', lat = ' + str(min_lat_grid) + u'\N{DEGREE SIGN})\n' + time_since_start + "- Coverage: " + format(percentage_covered, ".1f") + "%"

        ax_2d.set_title(ax_2d_title, weight = 'bold', fontsize  = (int)(fontsize_plot*1.1), y = 1.0)
        ax_2d.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
        ax_2d.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

        [i.set_linewidth(2) for i in ax_2d.spines.itervalues()] # change the width of the frame of the figure
        ax_2d.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
        plt.rc('font', weight='bold') ## make the labels of the ticks in bold

        # if lon_min_grid > 180:
        #     lon_min_grid_converted = lon_min_grid - 360
        # else:
        #     lon_min_grid_converted = lon_min_grid 
        # if lon_max_grid > 180:
        #     lon_max_grid_converted = lon_max_grid - 360
        # else:
        #     lon_max_grid_converted = lon_max_grid
        lon_min_grid_converted = lon_min_grid 
        lon_max_grid_converted = lon_max_grid 
        m = Basemap( projection       = 'cyl',
                     llcrnrlon        =  lon_min_grid_converted, #Lower Left  CoRNeR Longitude
                     urcrnrlon        =  lon_max_grid_converted , #Upper Right CoRNeR Longitude
                     llcrnrlat        = lat_min_grid  , #Lower Left  CoRNeR Latitude
                     urcrnrlat        = lat_max_grid,   #Upper Right CoRNeR Latitude
                     resolution       = 'i'  , # 'i', 'h'
                     suppress_ticks   = False,
                     ax = ax_2d,
                     )

        lon_arr = np.arange(lon_min_grid_converted,lon_max_grid_converted,width_cell)
        lat_arr = np.arange(lat_min_grid,lat_max_grid,width_cell)

        # draw grid
        data = cell_array_type_lon_lat_cumul[itime, : ,:]

        m.pcolormesh(lon_arr, lat_arr, data, cmap=cmap)

        specular_list = []
        point = namedtuple('point', ['x', 'y'])
        color = namedtuple('color', 'red green blue')
        specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))

        for isc in range(nb_sc):
            # Add on plot the specular points over one orbit
            for k in range(nb_spec):
                specular_list.append(specular)
                if itime > 3600:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,itime-3600:itime], lat_spec[k,isc,itime-3600:itime])
                else:
                    specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y =  m(lon_spec[k,isc,:itime], lat_spec[k,isc,:itime])
                specular_list[k+isc*nb_spec].point_plot = m.scatter(specular_list[k+isc*nb_spec].x, specular_list[k+isc*nb_spec].y, marker='o', s = 50, color = 'b')


    #    date_map = ax_2d.text((lon_min_grid_converted + lon_max_grid_converted)/2.,lat_min_grid + (lat_max_grid - lat_min_grid)/30., time_since_start + "\nCoverage: " + format(percentage_covered, ".1f") + "%", fontsize = fontsize_plot, weight = 'bold', horizontalalignment = 'center')




        m.drawcoastlines(linewidth=0.7, color='blue')
        m.drawmeridians(np.arange(lon_min_grid_converted,lon_max_grid_converted,width_cell) )
        m.drawparallels(np.arange(lat_min_grid,lat_max_grid,width_cell) )

        fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/2d_visu/' + input_filename.replace(".txt", "") + '_2d_coverage_itime_' + str(itime_count) + '.png'
        fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')



    raise Exception
    # AVERAGE OVER ALL LONGITUDES

    # Coverage vs time for different bin sizes
    #for ilat_center in range(nb_lat_center):#for ilat_center in range(nb_lat_center): # !!!!!!!!!!

    ilat_center = np.where(max_lat_grid_array == 30)[0][0]
    min_lat_grid = max_lat_grid_array[ilat_center]
    fig_title = 'Long. avg. temporal coverage VS time for different spatial resolutions - Grid at latitude ' + str(min_lat_grid) + u'\N{DEGREE SIGN}'
    y_label = 'Coverage (%)'
    x_label = 'Time (days)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    # for itype in range(nb_type_cell):#for itype in range(nb_type_cell): # !!!!!!!
    #     ax.plot(np.arange(0, nb_time, 1) / 3600. / 24., nb_cell_covered_average_over_lon[ilat_center, itype, :],  linewidth = 4, color = color_array[itype], label = str(width_cell_array[itype]) + u'\N{DEGREE SIGN}')


    ax.minorticks_off()
    xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
    xticks_label = []
    for i in range(len(xticks)):
        xticks_label.append( format( xticks[i], ".0f" ) )
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

    ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Width of bins", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))
    fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/new' + input_filename.replace(".txt", "") + '_coverage_vs_time_latitude_' + str(min_lat_grid).replace(".","_") + '_average_over_lon.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Width of bin vs Time to reach coverage goal for different latitude center of the window, for one coverage
    icoverage = 0
    fig_title = 'Long. avg. spatial VS temporal coverage for different latitudes of the grid - Coverage goal ' + str(coverage_array[icoverage]) + '%'
    y_label = 'Width of bin (lat/lon deg)'
    x_label = 'Time to reach coverage goal (days)'

    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')
    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    max_x = 0
    for ilat_center in range(nb_lat_center): # !!!!!!!!!!!# a value is negative if the coverage never reached the coverage goal. SO we don't want it on the plot
        ax.scatter(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage]/3600. / 24. > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage]/3600. / 24. > 0)[0]], linewidth = 2, color = color_array[ilat_center], marker = 'o', s = 100)
        ax.loglog(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage]/3600. / 24. > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage]/3600. / 24. > 0)[0]], linewidth = 3, color = color_array[ilat_center], label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN}')
        if np.max(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage]/3600. / 24. > 0)[0], icoverage]/3600. / 24.) > max_x:
            max_x = np.max(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage]/3600. / 24. > 0)[0], icoverage]/3600. / 24.)

    ax.minorticks_off()
    xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
    xticks_label = []
    for i in range(len(xticks)):
            xticks_label.append( format( xticks[i], ".0f" ) )
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

    ylim_save = [ax.get_ylim()[0], ax.get_ylim()[1]]
    yticks = width_cell_array
    yticks_label = []
    for i in range(len(yticks)):
        yticks_label.append( str( width_cell_array[i] ) )
    ax.yaxis.set_ticks(yticks)
    ax.yaxis.set_ticklabels(yticks_label, fontsize = fontsize_plot)#, rotation='vertical')

    ax.set_xlim([0, max_x*1.1])
    ax.set_ylim(ylim_save)

    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Grid lat.", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))

    fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/new/' + input_filename.replace(".txt", "") + '_diff_latitude_coverage_' + str(coverage_array[icoverage]).replace(".","_") + '_average_over_lon.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Width of bin vs Time to reach coverage goal for different coverage goals, at one latittude
    ilat_center = np.where(max_lat_grid_array == 30)[0][0]
    min_lat_grid = max_lat_grid_array[ilat_center]
    fig_title = 'Long. avg. spatial VS temporal coverage for different coverage goals - Grid at latitude ' + str(min_lat_grid) + u'\N{DEGREE SIGN}'
    y_label = 'Width of bin (lat/lon deg)'
    x_label = 'Time to reach coverage goal (days)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    max_x = 0
    for icoverage in range(nb_coverage):  # a value is negative if the coverage never reached the coverage goal. SO we don't want it on the plot        
        ax.scatter(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage] > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage] > 0)[0]], linewidth = 2, color = color_array[icoverage], marker = 'o', s = 100)
        ax.loglog(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage] > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage] > 0)[0]], linewidth = 3, color = color_array[icoverage], label = str(coverage_array[icoverage]) + "%")
        if np.max(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage] > 0)[0], icoverage]/3600. / 24.) > max_x:
            max_x = np.max(time_to_reach_coverage_average_over_lon[ilat_center, np.where(time_to_reach_coverage_average_over_lon[ilat_center, :, icoverage] > 0)[0], icoverage]/3600. / 24.)



    ax.minorticks_off()
    xticks = np.arange(0,np.max(nb_time) / 3600. / 24., 1)
    xticks_label = []
    for i in range(len(xticks)):
            xticks_label.append( format( xticks[i], ".0f" ) )
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticks_label, fontsize = fontsize_plot)

    ylim_save = [ax.get_ylim()[0], ax.get_ylim()[1]]
    yticks = width_cell_array
    yticks_label = []
    for i in range(len(yticks)):
        yticks_label.append( str( width_cell_array[i] ) )
    ax.yaxis.set_ticks(yticks)
    ax.yaxis.set_ticklabels(yticks_label, fontsize = fontsize_plot)#, rotation='vertical')

    ax.set_xlim([0, max_x*1.1])
    ax.set_ylim(ylim_save)



    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Cov. goal", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))


    fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/new/' + input_filename.replace(".txt", "") + '_diff_coverage_latitude_' + str(min_lat_grid).replace(".","_") + '_average_over_lon.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    raise Exception
    # ONE LONGITUDE IN PARTICULAR
    # Coverage vs time for different latitude center of the grid
    ilat_center = 3
    ilon_center = 0
    min_lat_grid = max_lat_grid_array[ilat_center]
    fig_title = 'Temporal coverage VS time for different spatial resolutions - Grid at latitude ' + str(min_lat_grid) + u'\N{DEGREE SIGN}'
    y_label = 'Coverage (%)'
    x_label = 'Time (days)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    # for itype in range(nb_type_cell):
    #     ax.plot(np.arange(0, nb_time, 1) / 3600. / 24., nb_cell_covered[ilat_center, itype, :],  linewidth = 4, color = color_array[itype], label = str(width_cell_array[itype]) + u'\N{DEGREE SIGN}')

    ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Width of bins", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))
    fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/coverage_vs_time_latitude_' + str(min_lat_grid).replace(".","_") + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Width of bin vs Time to reach coverage goal for different latitude center of the window, for one coverage
    icoverage = 3
    fig_title = 'Spatial VS temporal coverage for different latitudes of the grid - Coverage goal ' + str(coverage_array[icoverage]) + '%'
    y_label = 'Width of bin (lat/lon deg)'
    x_label = 'Time to reach coverage goal (days)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    for ilat_center in range(nb_lat_center): # a value is negative if the coverage never reached the coverage goal. SO we don't want it on the plot
        ax.scatter(time_to_reach_coverage[ilat_center, np.where(time_to_reach_coverage[ilat_center, :, icoverage]/3600. / 24. > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage[ilat_center, :, icoverage]/3600. / 24. > 0)[0]], linewidth = 2, color = color_array[ilat_center], marker = 'o', s = 100)
        ax.plot(time_to_reach_coverage[ilat_center, np.where(time_to_reach_coverage[ilat_center, :, icoverage]/3600. / 24. > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage[ilat_center, :, icoverage]/3600. / 24. > 0)[0]], linewidth = 3, color = color_array[ilat_center], label = str(max_lat_grid_array[ilat_center]) + u'\N{DEGREE SIGN}')

    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Cov. goal", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))
    ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])

    fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/diff_latitude_coverage_' + str(coverage_array[icoverage]).replace(".","_") + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  


    # Width of bin vs Time to reach coverage goal for different coverage goals, at one latittude
    min_lat_grid = max_lat_grid_array[ilat_center]
    fig_title = 'Spatial VS temporal coverage for different coverage goals - Grid at latitude ' + str(min_lat_grid) + u'\N{DEGREE SIGN}'
    y_label = 'Width of bin (lat/lon deg)'
    x_label = 'Time to reach coverage goal (days)'
    fig = plt.figure(num=None, figsize=(height_fig * ratio_fig_size, height_fig), dpi=80, facecolor='w', edgecolor='k')

    fig.suptitle(fig_title, y = 0.965,fontsize = (int)(fontsize_plot*1.1), weight = 'bold',)
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold
    gs = gridspec.GridSpec(1, 1)
    gs.update(left = 0.11, right=0.87, top = 0.93,bottom = 0.12, hspace = 0.01)
    ax = fig.add_subplot(gs[0, 0])

    ax.set_ylabel(y_label, weight = 'bold', fontsize  = fontsize_plot)
    ax.set_xlabel(x_label, weight = 'bold', fontsize  = fontsize_plot)

    [i.set_linewidth(2) for i in ax.spines.itervalues()] # change the width of the frame of the figure
    ax.tick_params(axis='both', which='major', labelsize=fontsize_plot, size = 10, width = 2, pad = 7) 
    plt.rc('font', weight='bold') ## make the labels of the ticks in bold

    for icoverage in range(nb_coverage):  # a value is negative if the coverage never reached the coverage goal. SO we don't want it on the plot        
        ax.scatter(time_to_reach_coverage[ilat_center, np.where(time_to_reach_coverage[ilat_center, :, icoverage] > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage[ilat_center, :, icoverage] > 0)[0]], linewidth = 2, color = color_array[icoverage], marker = 'o', s = 100)
        ax.plot(time_to_reach_coverage[ilat_center, np.where(time_to_reach_coverage[ilat_center, :, icoverage] > 0)[0], icoverage]/3600. / 24., width_cell_array[np.where(time_to_reach_coverage[ilat_center, :, icoverage] > 0)[0]], linewidth = 3, color = color_array[icoverage], label = str(coverage_array[icoverage]) + "%")

    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Cov. goal", fontsize = fontsize_plot)
    legend.get_title().set_fontsize(str(fontsize_plot))
    ax.set_xlim([0, np.max(nb_time) / 3600. / 24.])

    fig_save_name = '/Users/cbv/cygnss/coverage_temporal_spatial_relation/plot/diff_coverage_latitude_' + str(min_lat_grid).replace(".","_") + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  









    # for ilat_center in range(nb_lat_center):
    #     cell_array_lat = []
    #     for ilon_center in range(nb_lon_center):
    #         lon_grid_center = lon_grid_center_array[ilon_center]
    #         if lon_grid_center < 0:
    #             lon_grid_center  = lon_grid_center + 360 # convert to from 0 to + 360

    #         min_lat_grid = max_lat_grid_array              [ilat_center] # latitude of the center of the grid
    #         lat_min_grid = min_lat_grid - lat_width / 2.
    #         lat_max_grid = min_lat_grid               + lat_width / 2.
    #         lon_min_grid = (lon_grid_center - lon_width / 2.)%360
    #         lon_max_grid = (lon_grid_center + lon_width / 2. )%360
    #         cell_array_lon_lat = []
    #         for itype in range(nb_type_cell):
    #             print 'type cell ' + str(itype) + ' out of ' + str( nb_type_cell - 1 ) + ' / ilon ' + str(ilon_center) + ' out of ' + str(nb_lon_center-1) + ' / ilat ' + str(ilat_center)  + ' out of ' + str(nb_lat_center-1) 
    #             #itype = 1
    #             #icoverage = 0
    #             width_cell = width_cell_array[itype]

    #             nb_cell_lon = (int) ( ( lon_max_grid - lon_min_grid ) / width_cell ) 
    #             nb_cell_lat = (int) ( ( lat_max_grid - lat_min_grid ) / width_cell ) 
    #             nb_cell_in_grid = nb_cell_lon * nb_cell_lat
    #             cell_array_type_lon_lat = np.zeros([nb_cell_lat+1, nb_cell_lon+1])# +1 for safety
    #             icell_save = []
    #             icoverage = 0
    #             for itime in range(nb_time):
    #                 #print itime, nb_time
    #                 for isc in range(nb_sc):
    #                     for ispec in range(nb_spec):
    #                         if ( ( lon_spec[ispec, isc, itime] >= lon_min_grid ) & ( lon_spec[ispec, isc, itime] < lon_max_grid ) & ( lat_spec[ispec, isc, itime] >= lat_min_grid ) & ( lat_spec[ispec, isc, itime] < lat_max_grid  ) & ( lat_spec[ispec, isc, itime] < 9999999 ) ): # if the specular point in the window we're looking at
    #                         #if ( ( lon_spec[ispec, isc, itime] >= lon_min_grid ) & ( lon_spec[ispec, isc, itime] < lon_max_grid - width_cell ) & ( lat_spec[ispec, isc, itime] >= lat_min_grid ) & ( lat_spec[ispec, isc, itime] < lat_max_grid - width_cell ) & ( lat_spec[ispec, isc, itime] < 9999999 ) ): # if the specular point in the window we're looking at
    #                             icell_lat = (int) ( ( lat_spec[ispec, isc, itime] - lat_min_grid ) / width_cell )
    #                             icell_lon = (int) ( ( lon_spec[ispec, isc, itime] - lon_min_grid ) / width_cell )
    #                             if ([icell_lat, icell_lon] in icell_save) == False:
    #                                 icell_save.append([icell_lat, icell_lon])
    #                             cell_array_type_lon_lat[icell_lat, icell_lon] = cell_array_type_lon_lat[icell_lat, icell_lon] + 1

    #                 # nb_cell_covered_per_lat = np.zeros([nb_cell_lat])
    #                 # for icell_lat in range(nb_cell_lat):
    #                 #     nb_cell_covered_per_lat[icell_lat] = len(np.where( cell_array_type_lon_lat[icell_lat, :] > 0)[0])
    #                 nb_cell_covered[ilat_center, ilon_center, itype, itime] = len(icell_save) *100. / nb_cell_in_grid   #np.sum( nb_cell_covered_per_lat) / nb_cell_in_grid * 100
    #                 if icoverage < nb_coverage:
    #                     if  nb_cell_covered[ilat_center, ilon_center, itype, itime] >= coverage_array[icoverage]:
    #                         time_to_reach_coverage[ilat_center, ilon_center, itype, icoverage] = itime # since the specular file is every second, itime is moving one second by one second
    #                         icoverage = icoverage + 1

    #             cell_array_lon_lat.append( cell_array_type_lon_lat )
    #         cell_array_lat.append( cell_array_lon_lat )
    #     cell_array.append( cell_array_lat ) 
