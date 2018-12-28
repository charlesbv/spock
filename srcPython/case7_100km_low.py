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

earth_flattening    = 1/298.257223560; # Earth flattening coefficient (no unit)
earth_radius        = 6378.137; # mean equatorial radius (km)
earth_mu    = 398600.4418; # gravitational parameter (km^3/s^2)
earth_j2    = 1.081874e-3; # J2 zonal harmonic coefficient (no unit)
rad2deg = 180./np.pi
deg2rad = 1 / rad2deg
second2year = 3600 * 24 * 365.25


r1_old = np.array( [6337966.5365122000, 1889309.8033002000, 1889309.8033002000] ) / 1000.
v1_old = np.array( [-2957.1994197397, 4960.1786006789, 4960.1786006789] ) / 1000.

r2_old = np.array( [6338670.5812339000, 1888217.7447937000, 1888150.3615517000] ) / 1000.
v2_old = np.array( [-2955.2934327179, 4960.8101452466, 4960.6245630827] ) / 1000.

r1_new = r1_old - 100
r1_new_mag = np.linalg.norm(r1_new)
r1_old_mag = np.linalg.norm(r1_old)
v1_new = v1_old * np.sqrt( r1_old_mag / r1_new_mag )

r2_new = r2_old - 100
r2_new_mag = np.linalg.norm(r2_new)
r2_old_mag = np.linalg.norm(r2_old)
v2_new = v2_old * np.sqrt( r2_old_mag / r2_new_mag )

print r1_new * 1000. ,v1_new * 1000.
print r2_new * 1000. ,v2_new * 1000.

#cae 5
v5 = np.array( [-6384206.8367291000, -1809788.8923854000, -1809788.8923854000] ) / 1000
np.linalg.norm(v5) - 6378

# case 6
v6 = np.array([-6345736.9319327000, -1876218.4567378000, -1876218.4567378000]) / 1000
np.linalg.norm(v6) - 6378
