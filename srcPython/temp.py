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
isc = 6
nb_tle_for_sc = len(ecc[isc])
nb_sec_since_date_ref = np.zeros([nb_tle_for_sc])
for itle in range(nb_tle_for_sc):
    nb_sec_since_date_ref[itle] = ( datetime.strptime( date[isc][itle].split('.')[0], "%y%j" ) - date_ref ).total_seconds() + np.float( '0.' + date[isc][itle].split('.')[1] ) * 24 * 3600
save_last_tle_date[isc] = nb_sec_since_date_ref[-1]
x = nb_sec_since_date_ref 
y = ( np.array(sma[isc]) ) / ( np.array(sma[isc][0]) )

nslope = len(x) - 1
yslope = np.zeros([nslope])
xslope = x[:-1]
for itle in range(nslope):
    yslope[itle] = ( y[itle+1] - y[itle] ) / ( x[itle+1] - x[itle] )
