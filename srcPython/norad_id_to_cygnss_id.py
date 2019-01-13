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

def norad_id_to_cygnss_id(norad_id):
    if norad_id == '41884':
        cygnss_id = 'FM05'
    elif norad_id == '41885':
        cygnss_id = 'FM04'
    elif norad_id == '41886':
        cygnss_id = 'FM02'
    elif norad_id == '41887':
        cygnss_id = 'FM01'
    elif norad_id == '41888':
        cygnss_id = 'FM08'
    elif norad_id == '41889':
        cygnss_id = 'FM06'
    elif norad_id == '41890':
        cygnss_id = 'FM07'
    elif norad_id == '41891':
        cygnss_id = 'FM03'
    else:
        print "***! Not a CYGNSS spacecraft !***"
        
    return cygnss_id
