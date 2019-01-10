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

import numpy as np
from datetime import datetime
#import ipdb

def read_collision_file_new_with_dca( collision_filename ):
    #collision_filename = '/home/cbv/Propagator/run_collision/output/dca/dca_collision.txt'
    collision_file = open(collision_filename, "r")
    read_collision_file = collision_file.readlines()

    nb_ensembles = (int)( read_collision_file[0].split(":")[1].replace("\n", "")[:-1] )

    find_dt_collision = 0

    if read_collision_file[find_dt_collision][:17] != "No close approach":
        while ( ( read_collision_file[find_dt_collision].split(' ')[0].replace("\n","") != 'Time' ) | ( read_collision_file[find_dt_collision].split(' ')[1].replace("\n","") != 'step' ) ):
            find_dt_collision = find_dt_collision + 1
        dt_collision = np.float( read_collision_file[find_dt_collision].split(':')[1].split("s")[0][:-1] )

        find_nb_ca = 0
        while (read_collision_file[find_nb_ca].split(' ')[0].replace("\n","") != 'Number' ):
            find_nb_ca = find_nb_ca + 1
        nb_ca = (int)( read_collision_file[find_nb_ca].split(':')[1].split('(')[0].replace("\n","") )
        nb_ca_no_coll = (int)( read_collision_file[find_nb_ca].split('+')[1].split('w')[0].replace("\n","") )

        
        if nb_ca > 0:
            date = []
            nb_collisions_each_dt = []
            cpc_final = np.zeros([nb_ca])
            find_ca_section = 0
            tca = []
            dca = []
            tca_no_coll = []
            dca_no_coll = []
            for ica in range(nb_ca):
                date_sub = []
                nb_collisions_each_dt_sub = []
                while (read_collision_file[find_ca_section][0:9] != '#Detailed' ):
                    find_ca_section = find_ca_section + 1
                tca.append(read_collision_file[find_ca_section + 1].split(': ')[1].replace("\n","")[:-1] )
                dca.append(float(read_collision_file[find_ca_section + 2].split(': ')[1].split()[0] ))
                #cpc_final[ica] = np.float( read_collision_file[find_ca_section + 2].split(':')[1].replace("\n","")[:-1] ) # [:-1] to get rid of the final period
                if ica == 0:
                    start_span_tca = read_collision_file[find_ca_section + 5].split(': ')[1].split(' - ')[0].replace("\n","").replace(" ","")
                    start_span_tca = datetime.strptime(start_span_tca, "%Y-%m-%dT%H:%M:%S.%f")
                    end_span_tca = read_collision_file[find_ca_section + 5].split(': ')[1].split(' - ')[1].split('(')[0].replace("\n","").replace(" ","")[:-1]
                    end_span_tca = datetime.strptime(end_span_tca, "%Y-%m-%dT%H:%M:%S.%f")
                    nb_time_steps_collision = (int) ( (end_span_tca - start_span_tca).total_seconds() / dt_collision ) 
                    cpc = np.zeros([nb_ca, nb_time_steps_collision])
                for idt in range(nb_time_steps_collision):
                    date_temp = read_collision_file[idt + find_ca_section + 7].split()[0]
                    date_sub.append( datetime.strptime( date_temp, "%Y-%m-%dT%H:%M:%S.%f" ) )
                    nb_collisions_each_dt_sub.append( np.float( read_collision_file[idt + find_ca_section + 7].split()[3] ) )
                    if idt == 0:
                        cpc[ica, 0] = nb_collisions_each_dt_sub[-1]              
                    else:
                        cpc[ica, idt] = cpc[ica, idt-1] + nb_collisions_each_dt_sub[-1]             


                date.append(date_sub)
                nb_collisions_each_dt.append(nb_collisions_each_dt_sub)
                find_ca_section = find_ca_section + 6 + nb_time_steps_collision
                cpc_final[ica] = cpc[ica, -1]/ nb_ensembles**2
            collision_file.close()    
            cpc = cpc / nb_ensembles**2
            nb_collisions_each_dt = np.array(nb_collisions_each_dt)
            date = np.array(date)
        elif nb_ca_no_coll > 0:
            date = []
            nb_collisions_each_dt = []
            cpc = []
            cpc_final = []
            tca = []
            dca = []
            find_ca_section = 0
            tca_no_coll = []
            dca_no_coll = []
            for ica in range(nb_ca_no_coll):
                while (read_collision_file[find_ca_section][0:9] != '#Detailed' ):
                    find_ca_section = find_ca_section + 1
                tca_no_coll.append(read_collision_file[find_ca_section + 1].split(': ')[1].replace("\n","")[:-1] )
                dca_no_coll.append(float(read_collision_file[find_ca_section + 2].split(': ')[1].split()[0] ) )
            
        else:
            date = []
            nb_collisions_each_dt = []
            cpc = []
            cpc_final = []
            tca = []
            tca_no_coll = []
            dca = []
            dca_no_coll = []

    else: # no close approach
        date = []
        nb_collisions_each_dt = []
        cpc = []
        cpc_final = []
        tca = []
        tca_no_coll = []
        dca = []
        dca_no_coll = []

    return date, nb_collisions_each_dt, cpc, cpc_final, tca, tca_no_coll, dca, dca_no_coll
