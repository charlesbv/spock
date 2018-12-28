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
# this function creates the collision input file
# ASSUMPTIONS:
# - r and v in m and m/s
# r and v are both [2,3] vectors

def spock_collision(collision_filename, r, v, nb_ens, min_dist_ca, min_dist_coll): 
    collision_file = open(collision_filename, "w")

    print >> collision_file, "#STATE_ECI"
    print >> collision_file, '(' + str(r[0,0]) + '; ' + str(r[0,1]) + '; ' + str(r[0,2]) + ') (' +  str(v[0,0]) + '; ' + str(v[0,1]) + '; ' + str(v[0,2]) + ')'
    print >> collision_file, '(' + str(r[1,0]) + '; ' + str(r[1,1]) + '; ' + str(r[1,2]) + ') (' +  str(v[1,0]) + '; ' + str(v[1,1]) + '; ' + str(v[1,2]) + ')'
    print >> collision_file, ''
    print >> collision_file, '#COVARIANCE'
    print >> collision_file, '((4.7440894789163000000000; -1.2583279067770000000000; -1.2583279067770000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(-1.2583279067770000000000; 6.1279552605419000000000; 2.1279552605419000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(-1.2583279067770000000000; 2.1279552605419000000000; 6.1279552605419000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000010000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000001; 0.0000010000000000000000; -0.0000000000000000000001);\n(0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000001; -0.0000000000000000000001; 0.0000010000000000000000))\n((4.7439512624715000000000; -1.2582550000046000000000; -1.2582079265320000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(-1.2582550000046000000000; 6.1281039832865000000000; 2.1280243672750000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(-1.2582079265320000000000; 2.1280243672750000000000; 6.1279447542420000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n(0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000010000000000000000; 0.0000000000000000000000; 0.0000000000000000000000);\n (0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000010000000000000000; -0.0000000000000000000002);\n(0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; 0.0000000000000000000000; -0.0000000000000000000002; 0.0000010000000000000000))'

    print >> collision_file, ''
    print >> collision_file, '#NB_ENSEMBLES_COLLISION\n' + str(nb_ens) + '\n'
    print >> collision_file, '#MIN_DISTANCE_CLOSE_APPROACH\n' + str(min_dist_ca) + '\n'
    print >> collision_file, '#MIN_DISTANCE_COLLISION\n' + str(min_dist_coll)

    collision_file.close()

    return 0 
