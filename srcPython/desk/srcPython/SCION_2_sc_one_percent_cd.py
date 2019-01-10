import sys
from get_prop_dir import *
from read_input_file import *
from read_output_file import *
from matplotlib import pyplot as plt

plt.ion()

if len(sys.argv) > 2:
    main_input_file_name = get_prop_dir(1) + sys.argv[1] + '/input/main_input/' + sys.argv[2]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'
else:
    main_input_file_name = get_prop_dir(1) + 'run/input/main_input/' + sys.argv[1]  # ex of sys.argv[1]: cadre_tumbling_0panel.txt'    

# read input file
input_variables, order_input_variables = read_input_file(main_input_file_name)
dt = input_variables[2]; nb_steps = input_variables[3]; satellite_to_plot_path = input_variables[6]; satellite_to_plot = input_variables[7];
nb_satellites = input_variables[4]

# read output files
var_out_choose = ["position", "velocity"]
r = np.zeros([nb_satellites, nb_steps, 3])
v = np.zeros([nb_satellites, nb_steps, 3])
for isat in range(nb_satellites):
    var_out, var_out_name = read_output_file(satellite_to_plot_path[isat] + satellite_to_plot[isat], var_out_choose)
    r[isat, :, :] = var_out[1]
    v[isat, :, :] = var_out[2]

# Distances between the 2 sc
dist = np.zeros([nb_steps])
for istep in range(nb_steps):
    dist[istep] = np.linalg.norm(r[1, istep] - r[0, istep])

# Plot
fig, ax = plt.subplots()
ax.plot(dist)
