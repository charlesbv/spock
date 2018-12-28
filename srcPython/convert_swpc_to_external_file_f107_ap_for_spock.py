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
# This script converts a file from SWPC into the format for SpOCK for external files (Ap and F10.7). The script figures out itself if the input file is F10.7 or Ap (as long as it's a file downloaded from SWPC)
# ASSUMPTIONS

filename_to_convert = "2017Q3_DGD.txt"
file_to_convert = open(filename_to_convert)
read_file_to_convert = file_to_convert.readlines()
file_to_convert.close()

# Figure out if the file is Ap or F10.7: 1st line has 'Geomagnetic' if Ap, 'Solar' if F10.7
if 'Geomagnetic' in read_file_to_convert[0]:
    param = 'ap'
else:
    param = 'f107'

# skip header (first line not header starts with a '2' (for the year))
nb_header = 0
while (read_file_to_convert[nb_header][0] != '2'):
    nb_header = nb_header + 1
nb_header = nb_header + 1

