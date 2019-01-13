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

# this cript converts the norad id of CYGNSS to its FM01, FM02, ..., FM08 designation

def norad_id_to_cygnss_name(norad_id):
    cygnss_norad = ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
    cygnss_name_array = ['FM01', 'FM02', 'FM03', 'FM04', 'FM05', 'FM06', 'FM07', 'FM08']
    conversion_norad_to_cygnss_name = [4, 3, 1, 0, 7, 5, 6, 2]
    cygnss_name = cygnss_name_array[conversion_norad_to_cygnss_name[cygnss_norad.index(norad_id)]]
    return cygnss_name  
