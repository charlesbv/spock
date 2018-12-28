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
# this script converts the FM01, FM02, ..., FM08 designation to the NORAD id of CYGNSS

def cygnss_name_to_norad_id(cygnss_name):
    cygnss_norad_array = ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
    cygnss_name_array = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
    cygnss_norad = cygnss_norad_array[cygnss_name_array.index(cygnss_name)]#cygnss_name_array[conversion_norad_to_cygnss_name[cygnss_norad.index(norad_id)]]
    return cygnss_norad

