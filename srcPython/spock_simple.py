# This script creates a basic SpOCK main input file.
# The name of the file is given as arguement of the script
# If a file with the same name already exists,
# the prefix "spock_" is added to the name of the file

import sys
import os.path

filename = sys.argv[1]
if os.path.isfile(filename):
    filename = 'spock_' + filename
file = open(filename, "w+")

print >> file, " #TIME"
print >> file, "2019-01-01T00:00:00"
print >> file, "2019-01-02T00:00:00"
print >> file, "10"
print >> file, ""
print >> file, "#SPACECRAFT"
print >> file, "1"
print >> file, "0"
print >> file, "29"
print >> file, "-1"
print >> file, "cygnss_geometry_2016_acco09.txt"
print >> file, ""
print >> file, "#ORBIT"
print >> file, "oe"
print >> file, "500 35 0 0 0 0"
print >> file, ""
print >> file, "#ATTITUDE"
print >> file, "nadir"
print >> file, ""
print >> file, "#FORCES"
print >> file, "20"
print >> file, "drag solar_pressure moon_gravity sun_gravity"
print >> file, "static"
print >> file, "100"
print >> file, "100"
print >> file, "15"
print >> file, ""
print >> file, "#OUTPUT"
print >> file, "out/out"
print >> file, "60"
print >> file, ""
print >> file, "#SPICE"
print >> file, "/Users/cbv/cspice/data"
print >> file, ""
print >> file, "#DENSITY_MOD"
print >> file, "1"


file.close()
