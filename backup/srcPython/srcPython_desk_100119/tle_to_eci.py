import numpy as np
from sgp4.earth_gravity import wgs72old
from sgp4.io import twoline2rv

filename_tle = "fm5.txt"
file_tle = open(filename_tle)
read_file_tle = file_tle.readlines()
line1 = read_file_tle[0].replace("\r", "").replace("\n","") 
line2 = read_file_tle[1].replace("\r", "").replace("\n","")
sc_tle = twoline2rv(line1, line2, wgs72old)
r_eci_tle, v_eci_tle = sc_tle.propagate( sc_tle.epoch.year, sc_tle.epoch.month, sc_tle.epoch.day, sc_tle.epoch.hour, sc_tle.epoch.minute, sc_tle.epoch.second ) 

print r_eci_tle
