# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
from get_prop_dir import *
import os
from read_input_file import *
from read_output_file import *
from os import listdir
from os.path import isfile, join


# ASSUMPTIONS:
# - mpirun is in /data/cygnss/tools/spock/mpi_installation/bin/mpirun (if not then change the line that runs the propagator (line 'os.system("/data/cygnss/tools/spock/mpi_installation/bin/mpirun -np 1 spock " + orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + ".txt")')
# - the following STK results have to be in ./Examples/Validation:
    # - stk_leo_2body.txt,
    # - stk_leo_2body_drag.txt,
    # - stk_leo_aspherical_earth.txt,
    # - stk_leo_2body_third_body.txt,
    # - stk_leo_2body_solar_pressure.txt,
    # - stk_geo_2body.txt,
    # - stk_geo_2body_drag.txt,
    # - stk_geo_aspherical_earth.txt,
    # - stk_geo_2body_third_body.txt,
    # - stk_geo_2body_solar_pressure.txt,
    # - stk_heo_2body.txt,
    # - stk_heo_2body_drag.txt
    # - stk_heo_aspherical_earth.txt,
    # - stk_heo_2body_third_body.txt,
    # - stk_heo_2body_solar_pressure.txt
    # - baseline_results_validation_propagator.txt

# List simulations to validate the propagator
orbit_arr = ["leo","geo","heo"] # heo: high elliptical orbit
forces_for_orbit_arr = ["2body", "2body_drag", "aspherical_earth", "2body_third_body", "2body_solar_pressure"]

eci_orbit = ["(-4467.3083710453; -5053.5032161387; -427.679278085) (3.8279470371; -2.8757493806; -6.0045254246)",\
                 "(36607.3548590000; -20921.7217480000; -0.0000000000) (1.5256360000; 2.6694510000; 0.0000000000)",\
                 "(-1529.8942870000; -2672.8773570000; -6150.1153400000) (8.7175180000; -4.9897090000; 0.0000000000)"]
attitude_forces = ["nadir", "nadir", "nadir", "nadir", "sun_pointed"]
geometry_forces =  ["one_plate.txt", "one_plate.txt", "one_plate.txt", "one_plate.txt", "one_plate_for_solar_pressure.txt"]
gravity_forces = ["0", "0", "20", "0", "0"]
forces_forces = ["none", "drag", "none", "sun_gravity moon_gravity", "solar_pressure"]

# # Create run diectory where the validation is made
run_dir = get_prop_dir(1) + "run_validation"
if os.path.exists(run_dir):
    print "***! It looks like a run_validation directory already exists in " + get_prop_dir(1) + ". If you want to validate the propagator, first delete this directory and then run validation_propagator.py again. !***\n"; raise Exception

os.system("mkdir " + run_dir)
os.chdir(run_dir)

# Set up the run directory 
os.system("mkdir input")
os.system("mkdir input/main_input")
os.system("mkdir input/geometry")
os.system("cp ../Examples/geometry/One_plate/* input/geometry/")
os.system("mkdir output")
os.system("ln -s " + get_prop_dir(1) + "src/egm96_to360_not_norm.txt " + run_dir + "/input/egm96_to360_not_norm.txt")
os.system("ln -s ../spock .")

# Open file in which the results of the comparison with HPOP are written
file_results = open("./output/results_validation_propagator.txt", "w")
print >> file_results, "orbit_force max_diff_pos max_diff_vel (distances are in meters)\n"

# Run the different simulations
nb_orbit = len(orbit_arr)
nb_runs_for_orbit = len(forces_for_orbit_arr)
max_r_diff_mag = np.zeros([nb_orbit, nb_runs_for_orbit])
max_v_diff_mag = np.zeros([nb_orbit, nb_runs_for_orbit])
for iorbit in range(nb_orbit):
    for irun in range( nb_runs_for_orbit ):
        print orbit_arr[iorbit] + ": " + forces_for_orbit_arr[irun]
        ## Create the main input file
        main_input_filename = "input/main_input/" + orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + ".txt" 
        main_input_file = open( main_input_filename, "w")
        print >> main_input_file, \
        "#TIME\n\
2004-06-01T00:00:00\n\
2004-06-02T00:00:00\n\
5\n\
\n\
#SPACECRAFT\n\
1\n\
0\n\
100\n\
-1\n"\
+ geometry_forces[irun] + "\n\
\n\
#ORBIT\n\
state_eci\n"\
+ eci_orbit[iorbit] + "\n\
\n\
#ATTITUDE\n"\
+ attitude_forces[irun] + "\n\
\n\
#FORCES\n"\
+ gravity_forces[irun] + "\n" \
+ forces_forces[irun] + "\n\
static\n\
100\n\
100\n\
15\n\
\n\
#OUTPUT\n"\
+ orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + "\n\
60"
        main_input_file.close()
        ## Run the propagator
        os.system("/data/cygnss/tools/spock/mpi_installation/bin/mpirun -np 1 spock " + orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + ".txt")

        # Compare the position and velocity with HPOP (from STK)
        input_filename_prop_complete = run_dir + "/input/main_input/" +  orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + ".txt"
        input_variables, order_input_variables = read_input_file(input_filename_prop_complete)
        output_propagator_path = input_variables[6][0]; output_propgator_file = input_variables[7][0]
        to_output = ["position", "velocity"]
        out, out_var = read_output_file(output_propagator_path + output_propgator_file, to_output)
        r_eci_prop = out[1]; v_eci_prop = out[2]

        # Read position and velocity from HPOP
        filename_hpop = get_prop_dir(1) + "Examples/Validation/stk_" + orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + ".txt"
        file_hpop = open(filename_hpop)
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
        print >> file_results, orbit_arr[iorbit] + "_" + forces_for_orbit_arr[irun] + ' ' + '{0:.3f}'.format(np.max(r_diff_mag)) + ' ' + '{0:.3f}'.format(np.max(v_diff_mag))

        # save results to compare with baseline results
        max_r_diff_mag[iorbit, irun] = np.max(r_diff_mag)
        max_v_diff_mag[iorbit, irun] = np.max(v_diff_mag)

# Compare new results with baseline results
baseline_results_filename = get_prop_dir(1) + "Examples/Validation/baseline_results_validation_propagator.txt"
baseline_results_file = open(baseline_results_filename)
read_baseline_results_file = baseline_results_file.readlines()
max_r_diff_mag_baseline = np.zeros([nb_orbit, nb_runs_for_orbit])
max_v_diff_mag_baseline = np.zeros([nb_orbit, nb_runs_for_orbit])
print >> file_results, "\nComparison with baseline (" + baseline_results_filename + "):"
success = 1
for iorbit in range(nb_orbit):
    for irun in range( nb_runs_for_orbit ):
        max_r_diff_mag_baseline[iorbit, irun] = np.float( read_baseline_results_file[iorbit * nb_runs_for_orbit + irun + 2].split()[1] )
        max_v_diff_mag_baseline[iorbit, irun] = np.float( read_baseline_results_file[iorbit * nb_runs_for_orbit + irun + 2].split()[2] )
        if ( max_r_diff_mag_baseline[iorbit, irun] < np.float( '{0:.3f}'.format(max_r_diff_mag[iorbit, irun]) ) ):
            success = 0
            print >> file_results, "Run " + read_baseline_results_file[iorbit * nb_runs_for_orbit + irun+ 2].split()[0] + " has a position difference with HPOP of " + '{0:.3f}'.format(max_r_diff_mag[iorbit, irun]) + " m but the baseline comparison has a difference of only " + '{0:.3f}'.format(max_r_diff_mag_baseline[iorbit, irun]) + " m."
        if ( max_v_diff_mag_baseline[iorbit, irun] < np.float( '{0:.3f}'.format(max_v_diff_mag[iorbit, irun] ) ) ):
            success = 0
            print >> file_results, "Run " + read_baseline_results_file[iorbit * nb_runs_for_orbit + irun+ 2].split()[0] + " has a velocity difference with HPOP of " + '{0:.3f}'.format(max_v_diff_mag[iorbit, irun]) + " m/s but the baseline comparison has a difference of only " + '{0:.3f}'.format(max_v_diff_mag_baseline[iorbit, irun]) + " m/s."

if ( success == 1 ):
    print >> file_results, "GOOD NEWS! All the position and velocity differences with HPOP are at least smaller than the baseline results."
    print "\nPropagator validated with success!"
else:
    print "Oh oh! Bigger differences in position or velocity than the baseline have been found. Please have a look at " + get_prop_dir(1) + "run_validation/output/results_validation_propagator.txt for more information." 

file_results.close()
os.chdir("../srcPython")
