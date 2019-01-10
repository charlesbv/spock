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

file_switch = open("../../run_paper3/input/switch_for_msis.txt", "w+")
n = 24
switch_for_msis = np.zeros([n])
switch_for_msis = switch_for_msis.astype(int)
switch_for_msis[ 1 ] = 1 # F10.7 effect on mean remove increase
switch_for_msis[ 2 ] = 1 # time independent remove decrease
switch_for_msis[ 3 ] = 1 # symmetrical annual remove increase
switch_for_msis[ 4 ] = 1 # symmetrical semiannual remove decrease
switch_for_msis[ 5 ] = 1 # asymmetrical annual remove decrease
switch_for_msis[ 6 ] = 1 # asymmetrical semiannual remove decrease
switch_for_msis[ 7 ] = 1 # diurnal remove decrease
switch_for_msis[ 8 ] = 1 # semidiurnal remove increase
switch_for_msis[ 9 ] = 1 # daily ap [when this is set to -1 (!) the pointer ap_a in struct nrlmsise_input must point to a struct ap_array]
switch_for_msis[ 10 ] = 1 # all UT/long effects remove increase
switch_for_msis[ 11 ] = 1  # longitudinal remove increase DO NOT REMOVE (DOES NOT MAKE SENSE TO IGNORE IT)
switch_for_msis[ 12 ] = 1 # UT and mixed UT/long remove increase
switch_for_msis[ 13 ] = 1 # mixed AP/UT/LONG remove increase
switch_for_msis[ 14 ] = 1 # terdiurnal remove increase
switch_for_msis[ 15 ] = 1 # departures from diffusive equilibrium remove decrease meters
switch_for_msis[ 16 ] = 1 # all TINF var remove increases a lot
switch_for_msis[ 17 ] = 1  # all TLB var remove decrease
switch_for_msis[ 18 ] = 1 # all TN1 var remove nothing
switch_for_msis[ 19 ] = 1 # all S var remove increase
switch_for_msis[ 20 ] = 1 # all TN2 var remove nothing
switch_for_msis[ 21 ] = 1 # all NLB var remove increase
switch_for_msis[ 22 ] = 1 # all TN3 var remove nothing
switch_for_msis[ 23 ] = 1  # turbo scale height var remove nothing

for iswitch in range(1,n):
    print >> file_switch, switch_for_msis[ iswitch ] 

file_switch.close()
