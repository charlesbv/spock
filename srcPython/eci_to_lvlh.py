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
import numpy.linalg as li

def eci_to_lvlh(position_eci, velocity_eci, vector_eci):
    "This function transforms a vector from the intertial ECI coordinates to the LVLH coordinates "
    # Radial (1/2)
    radial_minus = -position_eci / li.norm(position_eci)

    # Cross-track
    cross_track_temp = np.cross(radial_minus, velocity_eci)
    cross_track = -cross_track_temp / li.norm(cross_track_temp)

    # Along-track
    along_track_temp = np.cross(cross_track, radial_minus)
    along_track = - along_track_temp / li.norm(along_track_temp)

    # Radial (2/2)
    radial = -radial_minus

    # Matrix conversion
    T_eci_2_lvlh = np.zeros([3, 3])
    T_eci_2_lvlh[0, :] = along_track
    T_eci_2_lvlh[1, :] = cross_track 
    T_eci_2_lvlh[2, :] = radial

    # Matrix * vector_eci
    vector_lvlh = np.dot(T_eci_2_lvlh , vector_eci)

    return vector_lvlh;
