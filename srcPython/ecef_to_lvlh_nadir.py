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

# same as ecef_to_lvlh but radial is pointing toward the Earth, ie nadir, not away from it (zenith)
def ecef_to_lvlh_nadir(position_ecef, velocity_ecef, vector_ecef):
    "This function transforms a vector from the interial ECEF coordinates to the LVLH coordinates "
    # Radial (1/2)
    radial_minus = -position_ecef / li.norm(position_ecef)

    # Cross-track
    cross_track_temp = np.cross(radial_minus, velocity_ecef)
    cross_track = cross_track_temp / li.norm(cross_track_temp)

    # Along-track
    along_track_temp = np.cross(cross_track, radial_minus)
    along_track = along_track_temp / li.norm(along_track_temp)

    # Radial (2/2)
    radial = radial_minus

    # Matrix conversion
    T_ecef_2_lvlh = np.zeros([3, 3])
    T_ecef_2_lvlh[0, :] = along_track
    T_ecef_2_lvlh[1, :] = cross_track 
    T_ecef_2_lvlh[2, :] = radial

    # Matrix * vector_ecef
    vector_lvlh = np.dot(T_ecef_2_lvlh , vector_ecef)

    return vector_lvlh;
