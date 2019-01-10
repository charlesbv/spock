from matplotlib import pyplot as plt
import numpy as np
import sys
from read_input_file import *
from read_output_file import *
from get_prop_dir import *
# This scrypt compares the ECI position in J2000 computed with the propagator to the one computed with HPOP in STK. The propagation is made for one day, for different satellites. We took the same satellites as the one in Gaylo et al., 2006 ('Testing of the Java Astrodynamics Toolkit Propagator').

# File with the results of the difference in positions and velocities between both propagators
file_results = open("./paper3_stk_results/results_comparison.txt", "w+")
print >> file_results, "name_run max_diff_pos max_diff_vel (distances are in meters)\n"

# Runs
input_filename_prop_array = ['iss','iss_msis','iss_gravity','iss_third_body', 'iss_solar_pressure', 'geo','geo_msis','geo_gravity','geo_third_body', 'geo_solar_pressure', 'molniya','molniya_msis','molniya_gravity','molniya_third_body', 'molniya_solar_pressure'] 
filename_hpop_array = ['ISS','ISS_msis','ISS_gravity','ISS_third_body', 'ISS_solar_pressure', 'geo','geo_msis','geo_gravity','geo_third_body', 'geo_solar_pressure', 'molniya','molniya_msis','molniya_gravity','molniya_third_body', 'molniya_solar_pressure'] 

nb_run = len(input_filename_prop_array)
##### BEGINNING OF LOOP
for irun in range(nb_run):
    # Read position and velocity from our propagator
    input_filename_prop = input_filename_prop_array[irun] + '.txt'
    input_filename_prop_complete = get_prop_dir(1) + "run_paper3/input/main_input/"  + input_filename_prop
    input_variables, order_input_variables = read_input_file(input_filename_prop_complete)
    output_propagator_path = input_variables[6][0]; output_propgator_file = input_variables[7][0]
    to_output = ["position", "velocity"]
    out, out_var = read_output_file(output_propagator_path + output_propgator_file, to_output)
    r_eci_prop = out[1]; v_eci_prop = out[2]

    # Read position and velocity from HPOP
    filename_hpop = filename_hpop_array[irun] + '.txt'
    file_hpop = open("./paper3_stk_results/" + filename_hpop)
    n_header = 7
    read_file_hpop = file_hpop.readlines()
    n_hpop = len(read_file_hpop) - n_header
    r_eci_hpop = np.zeros([n_hpop, 3]); v_eci_hpop = np.zeros([n_hpop, 3]); 
    for i in range(n_hpop):
        r_eci_hpop[i, 0] = read_file_hpop[i + n_header].split()[4]
        r_eci_hpop[i, 1] = read_file_hpop[i + n_header].split()[5]
        r_eci_hpop[i, 2] = read_file_hpop[i + n_header].split()[6]
        v_eci_hpop[i, 0] = read_file_hpop[i + n_header].split()[7]
        v_eci_hpop[i, 1] = read_file_hpop[i + n_header].split()[8]
        v_eci_hpop[i, 2] = read_file_hpop[i + n_header].split()[9]

    # Difference between positions and velocities
    r_diff = r_eci_hpop - r_eci_prop
    v_diff = v_eci_hpop - v_eci_prop

    r_diff_mag = np.zeros([n_hpop])
    v_diff_mag = np.zeros([n_hpop])
    for i in range(n_hpop):
        r_diff_mag[i] = np.linalg.norm(r_diff[i]) * 1000. # in m
        v_diff_mag[i] = np.linalg.norm(v_diff[i]) * 1000. # in m/s

    # Write results in file
    print >> file_results, input_filename_prop.split('.')[0] + ' ' + '{0:.3f}'.format(np.max(r_diff_mag)) + ' ' + '{0:.3f}'.format(np.max(v_diff_mag))

##### END OF LOOP

# Close file with results
file_results.close()

################################################################################# 
###################################### ISS ###################################### 
# INITIAL POSITION AND VELOCITY IN ECI J2000
  # CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = -4467.3083710453;  CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] =  -5053.5032161387; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] =  -427.6792780851;
  # CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = 3.8279470371;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = -2.8757493806;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] =  -6.0045254246;


# #2 body
# us = [770.586232, 5360.928715, 4042.807030] #770.585406 ,5360.928442 ,4042.807536]
# hpop = [770.586239 ,    5360.928711,    4042.807034] #770.585412,     5360.928443,     4042.807534] #[770.585413   ,  5360.928438  ,  4042.807540 ] #

# 2 body with msis but starting at 00:00  on june 1st and ending at 00:00 on june 2nd
us = [ 759.678696 ,5357.288996, 4049.450931 ]
hpop = [ 759.670701   ,  5357.286312  ,  4049.455800 ]

# # Spherical harmonics order and degree equal to 20
# us = [39.464148, 4885.624882 ,4653.642994]
# hpop = [39.463656    , 4885.624838  ,  4653.643086]

# # Third-body the Sun and the Moon
# us = [770.530914, 5360.902563, 4042.847865 ]
# hpop = [770.530920   ,  5360.902559  ,  4042.847869 ]

# # Solar pressure BUT eb careful, this time for our propagator we are Sun-pointing and the one plate is (0; 0; 1) (geometry file). DON'T KNOW YET WHICH COEFF TO PUT IN THE GEOMETRY FIOLE FOR THE REFLECTIVITIES
# us = [770.583878 ,5360.926564, 4042.807318]
# hpop = [770.585744  ,   5360.927296,    4042.806237]

# us = np.array(us)
# hpop = np.array(hpop)
# print np.linalg.norm(hpop - us) * 1000

################################################################################# 
###################################### GEO ###################################### 

# INITIAL POSITION AND VELOCITY IN ECI J2000  
  # CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[0] = 36607.3548590000;  CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[1] =  -20921.7217480000; CONSTELLATION->spacecraft[ii][eee].r_i2cg_INRTL[2] = -0.0000000000;
  # CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[0] = 1.5256360000;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[1] = 2.6694510000;  CONSTELLATION->spacecraft[ii][eee].v_i2cg_INRTL[2] =  0.0000000000;

# # 2 body
# us = [36961.885619, -20288.811600, -0.000000]
# hpop = [36961.885619 ,   -20288.811600 ,   -0.000076]

# # Third body perturbation
# us = [36955.051967, -20302.211961, -3.236235]
# hpop = [36955.051967 ,   -20302.211963 ,   -3.236310]

# # Spherical harmonics
# us = [36971.203521, -20271.891260, 0.002186]
# hpop = [36971.203534   , -20271.891207 ,   0.002111]

# # Solar radiation pressure
# us = [36962.030751, -20288.804515, -0.000027]
# hpop = [36961.970200 ,   -20288.807414  ,  -0.000085]

us = np.array(us)
hpop = np.array(hpop)
print np.linalg.norm(hpop - us) * 1000
