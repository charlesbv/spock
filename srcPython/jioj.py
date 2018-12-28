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







width_fig = 13.
height_fig = 1 width_fig * 3 /4 

fig = plt.figure(num=None, figsize=(width_fig, height_fig), dpi=80, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(1, 1)
gs.update(left=0.1, right=0.82, top = 0.90,bottom = 0.06)
fig.suptitle('', fontsize = 22)
plt.rc('font', weight='bold') ## make the labels of the ticks in bold
ax1 = fig.add_subplot(gs[0, 0])
        #ax1 = fig.add_subplot(111)
ax1.set_title('', weight = 'bold', fontsize = 20,  y = 1.008) 
ax1.tick_params(axis='both', which='major',  width = 2, pad = 6, labelsize=fontsize_plot, size = 10)
        #ax1.tick_params(axis='both', which='major', labelsize=16, size = 10, width = 2, pad = 7)
[i.set_linewidth(2) for i in ax1.spines.itervalues()] # change the width of the frame of the figure


# ax1.xaxis.set_major_locator(FixedLocator(np.arange(min_lon_array[zoom], max_lon_array[zoom]+1, step_lon[zoom])))
# ax1.yaxis.set_major_locator(FixedLocator(np.arange(min_lat_array[zoom], max_lat_array[zoom]+1, step_lat[zoom])))
# ax1.set_xticklabels(array_lon, color = 'k', fontsize = fontsize_plot, weight = 'bold')
# ax1.set_yticklabels(array_lat, color = 'k', fontsize = fontsize_plot, weight = 'bold')
ax1.set_xlabel('Longitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')
ax1.set_ylabel('Latitude (' + u'\N{DEGREE SIGN}'+ ')', fontsize = fontsize_plot, weight = 'bold')

# add_lat = 2.175e1
# add_lon = 1.56659e2 
# min_lon_collision = 0.003 + add_lon
# max_lon_collision = 0.006 + add_lon
# min_lat_collision = 0.002 + add_lat
# max_lat_collision = 0.0035 + add_lat

add_lat = 2.1751e1
add_lon = 1.56662e2
min_lon_collision = 0.0007 + add_lon 
max_lon_collision = 0.0026 + add_lon
min_lat_collision = 0.0020 + add_lat
max_lat_collision = 0.0024 + add_lat

# add_lat = 2.1751e1
# add_lon = 1.56661e2 
# min_lon_collision = 0.0016 + add_lon
# max_lon_collision = 0.0035 + add_lon
# min_lat_collision = 0.0022 + add_lat
# max_lat_collision = 0.00225 + add_lat


m = Basemap( projection       = 'cyl',
	     llcrnrlon        = min_lon_collision , #Lower Left  CoRNeR Longitude
	     urcrnrlon        = max_lon_collision  , #Upper Right CoRNeR Longitude
	     llcrnrlat        = min_lat_collision  , #Lower Left  CoRNeR Latitude
	     urcrnrlat        = max_lat_collision,   #Upper Right CoRNeR Latitude
	     resolution       = 'l'  ,
	     suppress_ticks   = False,
	     ax = ax1 ,
	     )

# ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
# ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
marker_array = ['o', 's']
color_array = ['b', 'r']

for irun in range(nb_runs):
    input_filename = get_prop_dir(1) + run_dir + 'input/main_input/' + run_list[irun]
    param_in, param_in_name = read_input_file(input_filename)
    output_file_path_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_path_list')]
    output_file_name_list = param_in[find_in_read_input_order_variables(param_in_name, 'output_file_name_list')]
    dt = param_in[find_in_read_input_order_variables(param_in_name, 'dt')]
    nb_satellites = param_in[find_in_read_input_order_variables(param_in_name, 'nb_satellites')]
    nb_steps = param_in[find_in_read_input_order_variables(param_in_name, 'nb_steps')]
    date_start = param_in[find_in_read_input_order_variables(param_in_name, 'date_start')]
    output_run_dir_name = output_file_path_list[0].split('/')[-3]
    collision_filename = get_prop_dir(1) + run_dir + 'output/' + output_run_dir_name + '/' + output_run_dir_name + '_collision.txt'
    date_collision, nb_collisions_each_dt, cpc, cpc_final, tca = read_collision_file( collision_filename )
    dt_coll = ( date_collision[0][1] - date_collision[0][0] ).total_seconds()
    nb_ca = cpc.shape[0]
    if nb_ca > 1:
        print "For now script works only if there is one single close approach. Here there are " + str(nb_ca) + " close approaches. The program will stop."; raise Exception

    # Figure out which time step in propagation the TCA corresponds to
    tca = datetime.strptime(tca[0], "%Y-%m-%dT%H:%M:%S.%f")
    index_tca_in_propagation = (int) ( ( tca - date_start ).total_seconds() / dt )
    #!!!!!!!!!!!!!! comment block below if you don't understand it. this is here because for these 2 runs the time step is one below the time step for all other runs
    if ( ( irun == 7 ) | ( irun == 8 ) ):
        index_tca_in_propagation = index_tca_in_propagation + 1
    #!!!!!!!!!!!!!! end of comment block below if you don't understand it

    # Read position of 2 sc 
    list_var = ["position", "latitude", "longitude"]
    eci_r = np.zeros([nb_satellites, nb_steps, 3])
    lat = np.zeros([nb_satellites, nb_steps])
    lon = np.zeros([nb_satellites, nb_steps])
    for isat in range(nb_satellites):
        out_var, order_var = read_output_file( output_file_path_list[isat] + output_file_name_list[isat], list_var  )
        date = out_var[find_in_read_input_order_variables(order_var, 'date')]
        latitude = out_var[find_in_read_input_order_variables(order_var, 'latitude')]
        longitude = out_var[find_in_read_input_order_variables(order_var, 'longitude')] 
        position = out_var[find_in_read_input_order_variables(order_var, 'position')]
        eci_r[isat, :, :] = position
        lat[isat, :] = latitude
        lon[isat, :] = longitude

        colorVal = scalarMap.to_rgba(irun)
        # if irun == index_median_run:
        #     # x, y =  m(lon[isat, index_tca_in_propagation-1:index_tca_in_propagation+2] , lat[isat, index_tca_in_propagation-1:index_tca_in_propagation+2])
        #     # m.plot(x, y,  marker='s', markersize=10,color = colorVal,markeredgecolor = 'none',linestyle="None", label = run_list[irun].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%"  )[0]
        # else:
    #!!!!!!!!!!!!!! comment block below if you don't understand it. this is here because for these 2 runs the time step is one below the time step for all other runs
        if ( ( irun == 7 ) | ( irun == 8 ) ):
            x, y =  m(lon[isat, index_tca_in_propagation-1: index_tca_in_propagation+1] , lat[isat,index_tca_in_propagation-1: index_tca_in_propagation+1])
    #!!!!!!!!!!!!!! end of comment block below if you don't understand it
        else:
            x, y =  m(lon[isat, index_tca_in_propagation: index_tca_in_propagation+2] , lat[isat,index_tca_in_propagation: index_tca_in_propagation+2])
        print isat, lon[isat, index_tca_in_propagation], lat[isat, index_tca_in_propagation]
        if isat == 0:
            m.plot(x, y,  marker='s', markersize=10,color = colorVal,markeredgecolor = 'none', label = run_list[irun].split('f107_')[1].split('_ap')[0].replace("_quartile", "") + "0%", linewidth = 2  )[0]
        else:
            m.plot(x, y,  marker='o', markersize=10,color = colorVal,markeredgecolor = 'none' , linewidth = 2 )[0]
legend = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="Quartile of\nF10.7/Ap dist.", fontsize = fontsize_plot)
legend.get_title().set_fontsize(str(fontsize_plot))

        # x, y =  m(lon[isat, index_tca_in_propagation] , lat[isat, index_tca_in_propagation])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'r')[0]
        # x, y =  m(lon[isat, index_tca_in_propagation+1] , lat[isat, index_tca_in_propagation+1])
        # m.plot(x, y,  marker=marker_array[isat], markersize=15,color = 'b')[0]



if save_results == 1:
    fig_save_name = 'zoom_position_sc_before_collision'
    fig_save_name = root_save_fig_name + fig_save_name + '.pdf'
    fig.savefig(fig_save_name, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')  
    os.system("rsync -av " + fig_save_name +" srbwks2014-0008.engin.umich.edu:./" + name_mission + '/' )

