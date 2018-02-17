import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
import h5py
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.animation as animation
from netCDF4 import Dataset
import os
from collections import *
import time
import datetime
from itertools import groupby
from operator import itemgetter

import urllib
import gdal
from gdalconst import *


############# USER INPUTS ###############
show_clouds = 1 # 0 for no clouds, 1 for clouds
show_specular_points = 1 # 0 for no specular points, 1 for specular points
zoom = 1  # 0 for full view, 1 for zoom
show_storm = 1 # this is for the cone of uncertainty and the circles to represent to storm (0: don't show, 1: show)
which_sat = [1,2,3,4]#[6,7,8]
nb_storms = 1

step_visualization = 1 # in seconds (time step to read the satellite files, IF the satellite files are written with one second time step). For example, 60 means that it's like if the satellite files were written with a time step of 60s. BE CAREFUL: if for example you set it to 35 then make sure that when_start_visu and when_stop_visu (see below) are consistent with this value!
when_start_visu = '2015-10-03T08:41:20'# 57925 # number of seconds after the start of the satellites files when to start the simulation
#35430 - 300  #- 300 # 39941 35430 index when to start the simulation (40685 start of sat1)
when_stop_visu = '2015-10-03T08:42:30' # number of seconds after the start of the satellites files when to stop the simulation
#  37016 (end of sat 6)  (40725 end of sat1) 
############# END OF USER INPUTS ###############

nb_satellites = len(which_sat)


# ######################################################
# used for the circle for the hurricanes:
def radius_for_tissot(dist_km):
    return np.rad2deg(dist_km/6378.) # this is calculating using the Haversine formula

def get_ellipse_coords(a=0.0, b=0.0, x=0.0, y=0.0, angle=0.0, k=2):
    """ Draws an ellipse using (360*k + 1) discrete points; based on pseudo code
    given at http://en.wikipedia.org/wiki/Ellipse
    k = 1 means 361 points (degree by degree)
    a = major axis distance,
    b = minor axis distance,
    x = offset along the x-axis
    y = offset along the y-axis
    angle = clockwise rotation [in degrees] of the ellipse;
        * angle=0  : the ellipse is aligned with the positive x-axis
        * angle=30 : rotated 30 degrees clockwise from positive x-axis
    """
    pts       = np.zeros((360*k+1, 2))

    beta      = -angle * np.pi/180.0
    sin_beta  = np.sin(beta)
    cos_beta  = np.cos(beta)
    alpha     = np.radians(np.r_[0.:360.:1j*(360*k+1)])
 
    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    
    pts[:, 0] = x + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = y + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts



######################################################
# Set up the 2D map
fig = plt.figure()#plt.figure(figsize=(10, 8))
background_color = (0/555,76./255,153/255.)
fig.set_facecolor(background_color)

# ax0 = fig.add_subplot(122)
# ax2 = fig.add_subplot(121)
#[left, bottom, width, height]
width_half = 0.40
ax0 = plt.axes([0.5+(0.5-width_half)/2,0.05,width_half,0.9])
#ax2 = plt.axes([0.05,0.5,width_half,0.3])

# ax1 = plt.axes([0.05,0.05,width_half,0.1125])
# ax2 = plt.axes([0.05,0.05+0.1125,width_half,0.1125])#ax2 = plt.axes([0.05,0.5,width_half,0.3])#ax2 = plt.axes([0.05,0.05+0.1125,width_half,0.1125])
# ax3 = plt.axes([0.05,0.05+0.1125*2,width_half,0.1125])
# ax4 = plt.axes([0.05,0.05+0.1125*3,width_half,0.1125])
# ax5 = plt.axes([0.05,0.05+0.1125*4,width_half,0.1125])
# ax6 = plt.axes([0.05,0.05+0.1125*5,width_half,0.1125])
# ax7 = plt.axes([0.05,0.05+0.1125*6,width_half,0.1125])
# ax8 = plt.axes([0.05,0.05+0.1125*7,width_half,0.1125])

# plt.setp( ax1.get_xticklabels(), visible=False)
# plt.setp( ax2.get_xticklabels(), visible=False)
# plt.setp( ax3.get_xticklabels(), visible=False)
# plt.setp( ax4.get_xticklabels(), visible=False)
# plt.setp( ax5.get_xticklabels(), visible=False)
# plt.setp( ax6.get_xticklabels(), visible=False)
# plt.setp( ax7.get_xticklabels(), visible=False)
# plt.setp( ax8.get_xticklabels(), visible=False)

# plt.setp( ax1.get_yticklabels(), visible=False)
# plt.setp( ax2.get_yticklabels(), visible=False)
# plt.setp( ax3.get_yticklabels(), visible=False)
# plt.setp( ax4.get_yticklabels(), visible=False)
# plt.setp( ax5.get_yticklabels(), visible=False)
# plt.setp( ax6.get_yticklabels(), visible=False)
# plt.setp( ax7.get_yticklabels(), visible=False)
# plt.setp( ax8.get_yticklabels(), visible=False)

# This should be called after all axes have been added
#fig.tight_layout()

min_lon = [-180, -90] #0
min_lat = [-90, 10] #-10
max_lon = [180, -50] #60
max_lat = [90,50]#40
step_lon = [60,10]#10
step_lat = [30, 10]#10
array_lon = [['180W', '120W', '60W', '0', '60E', '120E', '180E'],['90W', '80W', '70W', '60W', '50W']] #['0', '10E', '20E', '30E', '40E', '50E', '60E']]
array_lat = [['90S', '60S', '30S', 'EQ', '30N', '60N', '90N'],['10N', '20N', '30N','40N','50N']]#['10S', 'EQ', '10N', '20N', '30N', '40N']]
ax0.xaxis.set_major_locator(FixedLocator(np.arange(min_lon[zoom], max_lon[zoom]+1, step_lon[zoom])))
ax0.yaxis.set_major_locator(FixedLocator(np.arange(min_lat[zoom], max_lat[zoom]+1, step_lat[zoom])))
ax0.set_xticklabels(array_lon[zoom])
ax0.set_yticklabels(array_lat[zoom])
 

m = Basemap( projection       = 'cyl',
             llcrnrlon        = min_lon[zoom] , #Lower Left  CoRNeR Longitude
             urcrnrlon        = max_lon[zoom]  , #Upper Right CoRNeR Longitude
             llcrnrlat        = min_lat[zoom]  , #Lower Left  CoRNeR Latitude
             urcrnrlat        = max_lat[zoom],   #Upper Right CoRNeR Latitude
             resolution       = 'l'  ,
             suppress_ticks   = False,
             ax = ax0,
             )

m.drawcoastlines(linewidth=0.7, color='blue')
#m.bluemarble()                                                                                                               

########################################################
########################################################
if (show_clouds == 1):
# # Everything between the double lines of '#' are the insert from other code
# # The newest GOES images are always called 'latest':

    # #filnam_e = 'GoesEast1V_latest.tif'
    # # Now get the most current files from the website:
    # #testfile_e = urllib.URLopener()
    # #testfile_e.retrieve("ftp://satepsanone.nesdis.noaa.gov/GIS/GOESeast/"+filnam_e, filnam_e)

    filnam_e = 'IRCOMP_20151002_2100.tif' #'GoesEast1V2781845.tif'#'COMPOSITE0001.tif'

# ################################################################################
    # # PLOTTING GOES EAST
    ds = gdal.Open(filnam_e)
    
    #Inquire about the projection (may not be required)
    proj = ds.GetProjection()

    # Use GetGeoTransform() to obtain the geoinformation
    gt = ds.GetGeoTransform()

# Examine the extents of the image based on pixel resolution
    width  = ds.RasterXSize
    height = ds.RasterYSize
    xres   = gt[1]
    yres   = gt[5]

# Get the edge coordinates (note, values are already mid-pixel, so there is no need to shift by xres, yres)
    minx = gt[0]
    miny = gt[3] + width*gt[4] + height*yres 
    maxx = gt[0] + width*xres + height*gt[2]
    maxy = gt[3] 

# Now read in the pixel values of the whole image
    band = ds.GetRasterBand(1)
    image = band.ReadAsArray(0,0,width,height)
# value = image[42,94] <-- For example, the pixel brightness value at 42,94

# Establish a lat lon array for this image
    lonE = np.arange(minx,maxx,xres)  
    latE = np.arange(miny,maxy,-1*yres) # <-- Has to be -1*yres to be positive

#Generate a lon x lat grid for the tif file
    xE, yE = np.meshgrid(lonE, latE)

#Now the image has to be rotated 180 degrees before plotting
    image_r = np.zeros([height,width])
    for i in range(height):
        image_r[i,:] = image[(height-1)-i,:]

    Tf  = m.pcolormesh(xE,yE,image_r,cmap='Greys_r')

# ########################################################
# ########################################################

# ######################################################
if (show_storm == 1):
# # Read the storm interpolated files.  ASSUMPTION: all stroms are run at the same epoch and for the same amount of time

# ##
    directory_file = './storm_files/'
    storm_filenames = ['joaquin_interpolated.dat']
#nb_storms = len(storm_filenames)
        
## Loop over the storms
### Read first filename to initialize the size of x, y, z
    file = open(directory_file+storm_filenames[0],'r')
    a = file.readlines()
    nb_lines_storm = len(a) 
    line = ["" for x in range( nb_lines_storm )]
    time_storm_temp = ["" for x in range( (int)( nb_lines_storm / step_visualization )+1 )]
    time_storm = ["" for x in range( (int)( nb_lines_storm / step_visualization )+1 )]
    lon_storm,lat_storm,radius_uncertainty_storm = np.zeros([(int)(nb_lines_storm/step_visualization)+1,nb_storms]), np.zeros([(int)(nb_lines_storm/step_visualization)+1,nb_storms]), np.zeros([(int)(nb_lines_storm/step_visualization)+1,nb_storms])
    file.close()

### Now read each file
    for i in range(nb_storms):
	file = open(directory_file+storm_filenames[i],'r')
	a = file.readlines()
	for line_count in range(0,nb_lines_storm,step_visualization ):
		lon_storm_index_line = 3
		lat_storm_index_line = 2
		radius_storm_index_line = 7
		line[line_count] = (a[line_count]).split(' ')
		time_storm_temp[(int)(line_count / step_visualization)] = line[line_count][0]
		while (line[line_count][lon_storm_index_line] == ''):
			lon_storm_index_line = lon_storm_index_line + 1
			lat_storm_index_line = lat_storm_index_line + 1
			radius_storm_index_line = radius_storm_index_line + 1
		while (line[line_count][lat_storm_index_line] == ''):
			lat_storm_index_line = lat_storm_index_line + 1
			radius_storm_index_line = lat_storm_index_line + 1
		lon_storm[(int)(line_count / step_visualization),i] = float( line[line_count][lon_storm_index_line]  )
		if (lon_storm[(int)(line_count / step_visualization),i] > 180.0):
			lon_storm[(int)(line_count / step_visualization),i] = lon_storm[(int)(line_count / step_visualization),i] - 360.0
		lat_storm[(int)(line_count / step_visualization),i] = float( line[line_count][lat_storm_index_line]  )
		while (line[line_count][radius_storm_index_line] == ''):
			radius_storm_index_line = lat_storm_index_line + 1
		radius_uncertainty_storm[(int)(line_count / step_visualization),i] = float( line[line_count][radius_storm_index_line]  )
	file.close()

# just make some timef format conversion
    for i in range(0,nb_lines_storm,step_visualization):
        time_storm[(int)(i/step_visualization)] = time_storm_temp[(int)(i/step_visualization)][0:4]+'/'+time_storm_temp[(int)(i/step_visualization)][5:7]+'/'+time_storm_temp[(int)(i/step_visualization)][8:10]+' '+time_storm_temp[(int)(i/step_visualization)][11:19]


# Build the tuples for the visualization of the storms
    point = namedtuple('point', ['x', 'y'])
    point_circle = namedtuple('point_circle', ['x_circle', 'y_circle'])
    color = namedtuple('color', 'red green blue')
    storm_list = []
    name_storm = ["" for x in range(nb_storms)]
    for i in range(nb_storms):
	storm = namedtuple('storm',('name',) +  point._fields + point_circle._fields + color._fields + ('point_plot',) + ('point_circle_plot',) + ('marker_storm',))
	storm_list.append(storm)
    # name
	name_temp = 'Igor' # !!!! valid only if there is the only storm there is is Igor
	storm_list[i].name = name_temp
	name_storm[i] = storm_list[i].name # will be used to know which GPS is linked to which specular point
    # initial position
	storm_list[i].x, storm_list[i].y =  m(lon_storm[0,i], lat_storm[0,i])
    # color of the storm in the visualization
        storm_list[i].red   = 0
        storm_list[i].green = 0
        storm_list[i].blue  = 0
        storm_list[i].marker_storm = 'o'



        ell = get_ellipse_coords(a=radius_for_tissot(radius_uncertainty_storm[0,i]), b=radius_for_tissot(radius_uncertainty_storm[0,i]), x=lon_storm[0,i], y=lat_storm[0,i])
        storm_list[i].point_plot = m.plot(ell[:,0], ell[:,1], color='red',linewidth=3.5)[0]

#	storm_list[i].point_plot = m.plot(storm_list[i].x, storm_list[i].y,  marker=storm_list[i].marker_storm, markersize=1,color = (storm_list[i].red, storm_list[i].green, storm_list[i].blue))[0]
#	storm_list[i].point_circle_plot = m.tissot(lon_storm[0,i], lat_storm[0,i], radius_for_tissot(radius_uncertainty_storm[0,i]), 256, facecolor='b', alpha = 0.0)

# this plots all the circles of uncertainty during the entire animation. Like this, it gives a cone (that is the superimposition of all these circles of uncertainty)
    for i in range(nb_storms):
        for j in range(0, 3*24*3600,600):
            if float(a[j].split()[3])  > 180:
                m.tissot(float(a[j].split()[3]) - 360, float(a[j].split()[2]), radius_for_tissot(float(a[j].split()[7])), 256, facecolor='b', alpha = 0.002)
            else:
                m.tissot(float(a[j].split()[3]), float(a[j].split()[2]), radius_for_tissot(float(a[j].split()[7])), 256, facecolor='b', alpha = 0.002)



######################################################
# Read the satellite file(s). ASSUMPTION: all satellites are propagated at the same epoch and for the same amount of time

## Say which satellites you want to run (indicate the name of the file)
directory_file = './storm_files/'
string_same = 'interpolated_position_LLA_ECI_ECEF_'
satellite_filenames = []
for i in range(nb_satellites):
        satellite_filenames.append(string_same + 'CYGNSS_' + str(which_sat[i])+'.txt')

nb_cygnss = nb_satellites 
## Loop over the satellites
### Read first filename to initialize the size of x, y, z
file = open(directory_file+satellite_filenames[0],'r')
a = file.readlines()
nb_lines_header = 0
nb_lines = len(a)-nb_lines_header
line = ["" for x in range( nb_lines )]
time_sat_temp = ["" for x in range( (int)( nb_lines / step_visualization )+1 )]
time_sat      = ["" for x in range( (int)( nb_lines / step_visualization )+1 )]
lon,lat = np.zeros([(int)(nb_lines / step_visualization)+1,nb_cygnss]), np.zeros([(int)(nb_lines / step_visualization)+1,nb_cygnss])
file.close()
sat_over_storm = np.zeros([(int)(nb_lines / step_visualization)+1,nb_cygnss,nb_storms])
#which_sat_see_storm_and_when = []


### Now read each file
for i in range(nb_cygnss):
	file = open(directory_file+satellite_filenames[i],'r')
	a = file.readlines()
	for line_count in range(0, nb_lines, step_visualization):
		lon_index_line = 1
		lat_index_line = 2
		sat_over_storm_index_line = 10
		line[line_count] = (a[nb_lines_header+line_count]).split(' ')
		time_sat_temp[(int)(line_count / step_visualization)] = line[line_count][0] + ' ' + line[line_count][1] 
		while (line[line_count][lon_index_line] == ''):
			lon_index_line = lon_index_line + 1
			lat_index_line = lat_index_line + 1
                        sat_over_storm_index_line = sat_over_storm_index_line + 1
		while (line[line_count][lat_index_line] == ''):
			lat_index_line = lat_index_line + 1
                        sat_over_storm_index_line = sat_over_storm_index_line + 1
		lon[(int)(line_count / step_visualization),i] = float( line[line_count][lon_index_line]  )
		if (lon[(int)(line_count / step_visualization),i] > 180.0):
			lon[(int)(line_count / step_visualization),i] = lon[(int)(line_count / step_visualization),i] - 360.0
		lat[(int)(line_count / step_visualization),i] = float( line[line_count][lat_index_line]  )
#                 for ss in range(nb_storms):
#                     while (line[line_count][sat_over_storm_index_line] == ''):
#                         sat_over_storm_index_line = sat_over_storm_index_line + 1
#                     sat_over_storm[(int)(line_count / step_visualization),i,ss] = float(line[line_count][sat_over_storm_index_line])
#                     if ( sat_over_storm[(int)(line_count / step_visualization),i,ss] == 1 ):
#                         which_sat_see_storm_and_when.append([line_count,i,ss])

	file.close()

# !!!!!!!! REMOVE THE LIE BELOW 
for c in range(nb_cygnss):
    for i in range((int)(nb_lines / step_visualization)+1 - 120):
        lon[i,c] = lon[i+120,c]
        lat[i,c] = lat[i+120,c]

# some time format conversion
for i in range(0,nb_lines,step_visualization):
    time_sat[(int)(i/step_visualization)] = time_sat_temp[(int)(i/step_visualization)][0:4]+'/'+time_sat_temp[(int)(i/step_visualization)][5:7]+'/'+time_sat_temp[(int)(i/step_visualization)][8:10]+' '+time_sat_temp[(int)(i/step_visualization)][11:19]


time_step = (datetime.datetime(*time.strptime(time_sat[1], "%Y/%m/%d %H:%M:%S")[0:6]) - datetime.datetime(*time.strptime(time_sat[0], "%Y/%m/%d %H:%M:%S")[0:6])).seconds # note: the time_step for the storms and the satellites is always the same (result of the propagator)

#nb_times_storm_is_seen = len(which_sat_see_storm_and_when)

if ( show_storm == 0 ):
    time_storm = time_sat
    nb_lines_storm = nb_lines

# from here to (*) don't really worry about. This is just in case the files for the satellites and the ones for the storms don't start at the same time
delta_time_between_stop_storm_and_stop_constellation = datetime.datetime(*time.strptime(time_storm[(int)(nb_lines_storm / step_visualization)-1], "%Y/%m/%d %H:%M:%S")[0:6]) - datetime.datetime(*time.strptime(time_sat[(int)(nb_lines/step_visualization)-1], "%Y/%m/%d %H:%M:%S")[0:6])
delta_time_between_stop_storm_and_stop_constellation_in_seconds = delta_time_between_stop_storm_and_stop_constellation.days*24*3600 + delta_time_between_stop_storm_and_stop_constellation.seconds
nb_time_steps_between_stop_storm_and_stop_constellation = delta_time_between_stop_storm_and_stop_constellation_in_seconds / time_step


start_array = [ datetime.datetime(*time.strptime(time_storm[0], "%Y/%m/%d %H:%M:%S")[0:6]), datetime.datetime(*time.strptime(time_sat[0], "%Y/%m/%d %H:%M:%S")[0:6])]
stop_array = [ datetime.datetime(*time.strptime(time_storm[(int)(nb_lines_storm/step_visualization)-1], "%Y/%m/%d %H:%M:%S")[0:6]), datetime.datetime(*time.strptime(time_sat[(int)(nb_lines/step_visualization)-1], "%Y/%m/%d %H:%M:%S")[0:6]) ]
start_visualization = max( start_array)# we want to start the simulation at max( time_sat[0], time_storm[0] ) (because we need data for both the storms and the satellites)
stop_visualization = min( stop_array ) # we want to stop the simulation at min( time_sat[0], time_storm[0] ) (because we need data for both the storms and the satellites)

delta_time_between_starts =  start_visualization -start_array[start_array.index(min(start_array))] 
delta_time_between_starts_in_seconds = delta_time_between_starts.days * 24 * 3600 + delta_time_between_starts.seconds
nb_time_steps_between_starts = delta_time_between_starts_in_seconds / time_step

delta_time_between_stops = stop_array[stop_array.index(max(stop_array))]  - stop_visualization
delta_time_between_stops_in_seconds = delta_time_between_stops.days * 24 * 3600 + delta_time_between_stops.seconds
nb_time_steps_between_stops = delta_time_between_stops_in_seconds / time_step

if ( start_array.index(max(start_array)) == 0 ):
    who_started_last = 'storms'
else:
    who_started_last = 'satellites'

if ( stop_array.index(min(stop_array)) == 0 ):
    who_stopped_first = 'storms'
else:
    who_stopped_first = 'satellites'
# (*)


######################################################
######################################################
# TUPLES
# Build the tuples for the visualization of the satellites
point = namedtuple('point', ['x', 'y'])
color = namedtuple('color', 'red green blue')
spacecraft_list = []
name_satellite = ["" for x in range(nb_cygnss)]
for i in range(nb_cygnss):
	spacecraft = namedtuple('spacecraft',('name',) +  point._fields + ('type',) + color._fields + ('point_plot',) + ('marker_spacecraft',))
	spacecraft_list.append(spacecraft)
    # name
	name_temp = satellite_filenames[i].replace(".txt","")
	spacecraft_list[i].name = name_temp
	name_satellite[i] = spacecraft_list[i].name # will be used to know which GPS is linked to which specular point
    # initial position
	spacecraft_list[i].x, spacecraft_list[i].y =  m(lon[0,i], lat[0,i])
    # type of spacecraft (CYGNSS or GPS)
	if name_temp[len(name_temp)-3:len(name_temp)] == "GPS":
		spacecraft_list[i].type = "GPS"
	elif name_temp[len(name_temp)-8:len(name_temp)-2] == "CYGNSS":
		spacecraft_list[i].type = "CYGNSS"
	else:
		print "Spacecraft "+str(i)+" is neither a CYGNSS nor a GPS satellite. The program will stop."
		quit()
    # color of the spacecraft in the visualization
	if (spacecraft_list[i].type == "CYGNSS"): # make CYGNSS a pink diamond
		spacecraft_list[i].red   = 1
		spacecraft_list[i].green = 0
		spacecraft_list[i].blue  = 1
		spacecraft_list[i].marker_spacecraft = 'o'
    # point on the plot
	if (spacecraft_list[i].type == "GPS"): # the color of each GPS is the same as the color of the specular it is linked to (and plain circles)
		spacecraft_list[i].red   = color_array[i-nb_cygnss,0] / 255.
		spacecraft_list[i].green = color_array[i-nb_cygnss,1] / 255.
		spacecraft_list[i].blue  = color_array[i-nb_cygnss,2] / 255.
		spacecraft_list[i].marker_spacecraft = 'o'
	spacecraft_list[i].point_plot = m.plot(spacecraft_list[i].x, spacecraft_list[i].y,  marker=spacecraft_list[i].marker_spacecraft, markersize=15,color = (spacecraft_list[i].red, spacecraft_list[i].green, spacecraft_list[i].blue))[0]

if (show_specular_points == 1):
# # # Specular analysis: 
    nb_gps = 32
    color_array = np.zeros([nb_gps, 3])
    color_array[0,:3] = [0, 0, 255]; color_array[11,:3] = [0, 255, 0]; color_array[21,:3] = [255, 0, 0]; 
    color_array[1,:3] = [0, 255, 255]; color_array[12,:3] = [255, 0, 255]; color_array[22,:3] = [255, 255, 0]; 
    color_array[2,:3] = [0, 0, 128]; color_array[13,:3] = [0, 128, 0]; color_array[23,:3] = [128, 0, 0]; 
    color_array[3,:3] = [0, 128, 128]; color_array[14,:3] = [128, 0, 128]; color_array[24,:3] = [128, 128, 0]; 
    color_array[4,:3] = [0, 0, 64]; color_array[15,:3] = [0, 64, 0]; color_array[25,:3] = [64, 0, 0]; 
    color_array[5,:3] = [0, 64, 64]; color_array[16,:3] = [64, 0, 64]; color_array[26,:3] = [64, 64, 0]; 
    color_array[6,:3] = [0, 0, 192]; color_array[17,:3] = [0, 192, 0]; color_array[27,:3] = [192, 0, 0]; 
    color_array[7,:3] = [0, 192, 192]; color_array[18,:3] = [192, 0, 192]; color_array[28,:3] = [192, 192, 0]; 
    color_array[8,:3] = [0, 0, 96]; color_array[19,:3] = [0, 96, 0]; color_array[29,:3] = [96, 0, 0]; 
    color_array[9,:3] = [0, 96, 96]; color_array[20,:3] = [96, 0, 96]; color_array[30,:3] = [96, 96, 0]; 
    color_array[10,:3] = [0, 48, 48];

    duration_propagation = 24L # in hours                                                                         
    duration_specular = 24L # in hours                                                                            
    duration_simu = min([duration_propagation, duration_specular])
    time_max = 1438L * 60# int(np.floor( duration_simu*3600L / 1.0 )) # the time step for the specular points is always 0
    specular_names = [] # name of the GPS to which the specular point is linked

## Loop over the satellites
    nb_lines_spec  = time_max+1
    time_spec      = ["" for x in range( nb_lines_spec )]
    time_spec_temp = ["" for x in range( nb_lines_spec )]
    lon_spec,lat_spec, gain_spec, distance_spec_to_storm, specular_over_storm = np.zeros([nb_lines_spec,nb_cygnss,nb_gps]), np.zeros([nb_lines_spec,nb_cygnss,nb_gps]), np.zeros([nb_lines_spec,nb_cygnss,nb_gps]), np.zeros([nb_lines_spec,nb_cygnss,nb_gps, nb_storms]), np.zeros([nb_lines_spec,nb_cygnss,nb_gps, nb_storms])
cygnss_spec_filenames = ['results_specular_see_storm_CYGNSS_']
which_specular_point_see_storm_and_when = []

if (show_specular_points == 1):
    for j in range(nb_cygnss):
        line_spec = ["" for x in range( ( nb_lines_spec * step_visualization+ 4 ) * nb_gps ) ]
	file = open(directory_file+cygnss_spec_filenames[0]+str(which_sat[j])+'.txt','r')
	a = file.readlines()
	for i in range(nb_gps):
		specular_names.append(a[ ( i * ( nb_lines_spec * step_visualization + 4 ) + 2 ) -2].split(' ')[2] + '_' + a[ ( i * ( nb_lines_spec + 4 ) + 2 ) -2].split(' ')[3])
		for line_count in range(0,nb_lines_spec*step_visualization,step_visualization):

			lon_index_line = 1
			lat_index_line = 2
			gain_index_line = 3
                        distance_index_line = 4
                        specular_over_storm_index_line = 5
			line_spec[line_count] = (a[ ( i * ( nb_lines_spec * step_visualization + 4 ) + 2 ) + line_count]).split(' ')
                        time_spec_temp[(int)(line_count / step_visualization)] = line_spec[line_count][0]
			while (line_spec[line_count][lon_index_line] == ''):
				lon_index_line = lon_index_line + 1
				lat_index_line = lat_index_line + 1
				gain_index_line = gain_index_line + 1
                                distance_index_line = distance_index_line + 1
                                specular_over_storm_index_line = specular_over_storm_index_line + 1
			while (line_spec[line_count][lat_index_line] == ''):
				lat_index_line = lat_index_line + 1
				gain_index_line = gain_index_line + 1
                                distance_index_line = distance_index_line + 1
                                specular_over_storm_index_line = specular_over_storm_index_line + 1
			while (line_spec[line_count][gain_index_line] == ''):
				gain_index_line = gain_index_line + 1
                                distance_index_line = distance_index_line + 1
                                specular_over_storm_index_line = specular_over_storm_index_line + 1
			lon_spec[(int)(line_count / step_visualization),j,i] = float( line_spec[line_count][lon_index_line]  )
			if (lon_spec[(int)(line_count / step_visualization),j,i] > 180.0):
				lon_spec[(int)(line_count / step_visualization),j,i] = lon_spec[(int)(line_count / step_visualization),j,i] - 360.0
			lat_spec[(int)(line_count / step_visualization),j,i] = float( line_spec[line_count][lat_index_line]  )
			gain_spec[(int)(line_count / step_visualization),j,i] = float( line_spec[line_count][gain_index_line]  )
                        for ss in range(nb_storms):
                            while (line_spec[line_count][distance_index_line] == ''):
                                distance_index_line = distance_index_line + 1
                                specular_over_storm_index_line = specular_over_storm_index_line + 1
                            while (line_spec[line_count][specular_over_storm_index_line] == ''):
                                specular_over_storm_index_line = specular_over_storm_index_line + 1
                            distance_spec_to_storm[(int)(line_count / step_visualization),j,i,ss] = float( line_spec[line_count][distance_index_line]  )
                            specular_over_storm[(int)(line_count / step_visualization),j,i,ss] = float( line_spec[line_count][specular_over_storm_index_line]  )
                            if ( specular_over_storm[(int)(line_count / step_visualization),j,i,ss] == 1 ):
                                which_specular_point_see_storm_and_when.append([line_count,j,i,ss])

	file.close()

    nb_times_storm_is_seen = len(which_specular_point_see_storm_and_when)
## Build the tuples for the visualization of the satellites
    specular_list = []
# !!! need to change if more than one CYGNSS
    for j in range(nb_cygnss):
        for i in range(nb_gps):
	    specular = namedtuple('specular',('name',) +  point._fields  + color._fields + ('point_plot',))
	    specular_list.append(specular)
    # name 
	    specular_list[i+j*nb_gps].name = specular_names[i]
    # initial position                                            
	    specular_list[i+j*nb_gps].x, specular_list[i+j*nb_gps].y =  m(lon_spec[0,j,i], lat_spec[0,j,i]) 
    # color of specular point = color of the GPS satellite it is linked to
# 	for k in range(nb_cygnss):
# 		if name_satellite[k] == specular_list[i+j*nb_gps].name.strip():
# 			specular_list[i+j*nb_gps].red   = 0#spacecraft_list[k].red
# 			specular_list[i+j*nb_gps].green = 0#spacecraft_list[k].green
# 			specular_list[i+j*nb_gps].blue  = 0#spacecraft_list[k].blue
    # point on the plot
	    specular_list[i+j*nb_gps].point_plot = m.plot(specular_list[i+j*nb_gps].x, specular_list[i+j*nb_gps].y, marker='o', markersize = 15, color = 'chartreuse', fillstyle = 'none', mew = 2)[0]#fillstyle='none',


# # ######################################################
# # ######################################################
when_start_visu = time_spec_temp.index(when_start_visu)
when_start_visu = (int)( when_start_visu / step_visualization ) 
when_stop_visu = time_spec_temp.index(when_stop_visu)
when_stop_visu = (int)( when_stop_visu / step_visualization ) 

# COVERAGE
times_when_storm_is_seen_temp = np.zeros([nb_times_storm_is_seen])
for i in range(nb_times_storm_is_seen):   
    times_when_storm_is_seen_temp[i] = (int)( which_specular_point_see_storm_and_when[i][0] / step_visualization)
times_when_storm_is_seen = np.sort(times_when_storm_is_seen_temp)

which_specular_gain_not_zero_and_above_storm = np.zeros([(int)(nb_lines_spec / step_visualization)+1, nb_cygnss, 4])
prob_spec = np.zeros([(int)(nb_lines_spec / step_visualization)+1, nb_cygnss, nb_gps])
nb_specular_points_gain_not_zero_and_above_storm_at_a_given_time = np.zeros([(int)(nb_lines_spec / step_visualization)+1, nb_cygnss])
sve_all_prop = []
if (show_storm == 1):
    for i in times_when_storm_is_seen:
        if (who_started_last == 'storms'):
            sigma_gaussian = radius_uncertainty_storm[(int)( i / step_visualization) +nb_time_steps_between_starts,0]
        else:
            sigma_gaussian = radius_uncertainty_storm[(int)( i / step_visualization) ,0]
        for c in range(nb_cygnss):
            nb_specular_points_gain_not_zero_and_above_storm_at_a_given_time[(int)( i / step_visualization),c] =  len(np.where( specular_over_storm[(int)( i / step_visualization),c,:,0] == 1 )[0])
            for s in range(min(4, len(np.where( specular_over_storm[(int)( i / step_visualization),c,:,0] == 1 )[0]))):
                which_specular_gain_not_zero_and_above_storm[(int)( i / step_visualization),c,s] = np.where( specular_over_storm[(int)( i / step_visualization),c,:,0] == 1)[0][s]
                distance_to_center = distance_spec_to_storm[(int)( i / step_visualization),c,which_specular_gain_not_zero_and_above_storm[(int)( i / step_visualization),c,s],0]
                prob_spec_temp = np.exp( np.log(2./3.) * distance_to_center**2 / (sigma_gaussian**2) ) 
 #( 1 / np.sqrt(2 * np.pi * sigma_gaussian**2) ) * np.exp( - distance_to_center**2 / (2*sigma_gaussian**2) )
                prob_spec[(int)( i / step_visualization), c, which_specular_gain_not_zero_and_above_storm[(int)( i / step_visualization),c,s]] = prob_spec_temp
                sve_all_prop.append(prob_spec_temp)
            
#  # sigma_gaussian = radius of uncertainty (so that if the specular point is a distance equal to the radius of uncertainty, the probability of being above the storm is 0.683 (in NOAA they define the circle of uncertainty so that if we are at a distance equal to the radius of uncertainty then the probability of being above the storm is 2/3 (so this is close to 0.683))

# # PLOT flights over
# ## TIMES WHEN STORM IS SEEN
#list_ax = [ax1, ax2]#, ax3, ax4, ax5, ax6, ax7, ax8] 
# For each specular point of each satellite, find the consecutive times when the specular point in above the storm (with a non-zero gain). Put these times in a list called "list_consecutive_times_when_specular_above_storm". The beginning and the end of each serie of times is stored in time_start_spec_above_storm and time_stop_spec_above_storm. For instance, time_start_spec_above_storm is the list of beginning time for each specular point of each satellite (it is a list (satellites) of list (specular points) of list (beginning times))
time_start_spec_above_storm = []
time_stop_spec_above_storm = []
for c in range(nb_cygnss):
    sub_list_start = []
    sub_list_stop = []
    for g in range(nb_gps):
        when_spec_is_above_storm_and_gain_not_zero = np.where(prob_spec[:,c,g] > 0)[0]
        list_consecutive_times_when_specular_above_storm = []
        for k, l in groupby(enumerate(when_spec_is_above_storm_and_gain_not_zero), lambda (i, x): i-x): # find the consecutive times when this specular points is above the storm (with a non zero gain)
            list_consecutive_times_when_specular_above_storm.append(map(itemgetter(1), l))
        sub_sub_list_start = []
        sub_sub_list_stop = []
        for e in range(len(list_consecutive_times_when_specular_above_storm)):
            sub_sub_list_start.append(list_consecutive_times_when_specular_above_storm[e][0])
            sub_sub_list_stop.append(list_consecutive_times_when_specular_above_storm[e][-1])
        sub_list_start.append(sub_sub_list_start)
        sub_list_stop.append(sub_sub_list_stop)
     #   print list_consecutive_times_when_specular_above_storm, len(list_consecutive_times_when_specular_above_storm)
    time_start_spec_above_storm.append(sub_list_start)
    time_stop_spec_above_storm.append(sub_list_stop)

# Make a list of sub-list that gather the moments when specular points fly over the storm within a same orbit
## list of starts of fly over
list_start_spec_fly_over_storm = []
all_start_for_given_sat_sorted = []
for c in range(nb_cygnss):
    all_start_for_given_sat = []
    sub_list_start_spec_fly_over_storm = []
    for g in range(nb_gps):
        for e in range(len( time_start_spec_above_storm[c][g])):
            all_start_for_given_sat.append(time_start_spec_above_storm[c][g][e])
    all_start_for_given_sat_sorted.append(np.sort(all_start_for_given_sat))
    skip_spec2 = 0
    h = 0
    while (h < len(all_start_for_given_sat_sorted[c]) ):
        hh = h + skip_spec2
        skip_spec = 0
        if (hh < len(all_start_for_given_sat_sorted[c])-1):
            sub_sub_list_start_spec_fly_over_storm = []
            sub_sub_list_start_spec_fly_over_storm.append(all_start_for_given_sat_sorted[c][hh])
            while (  (hh+1+skip_spec < len(all_start_for_given_sat_sorted[c]) - 1) & ( (all_start_for_given_sat_sorted[c][hh+1+skip_spec]-all_start_for_given_sat_sorted[c][hh+skip_spec])  < 20*60/step_visualization ) ) : 
                sub_sub_list_start_spec_fly_over_storm.append(all_start_for_given_sat_sorted[c][hh+1+skip_spec])
                skip_spec = skip_spec + 1
                skip_spec2 = skip_spec2 + 1
            sub_list_start_spec_fly_over_storm.append(sub_sub_list_start_spec_fly_over_storm)
        if (hh == len(all_start_for_given_sat_sorted[c])-1):
            if ( all_start_for_given_sat_sorted[c][hh] - all_start_for_given_sat_sorted[c][hh-1] < 20*60/step_visualization ):
                sub_list_start_spec_fly_over_storm[-1].append(all_start_for_given_sat_sorted[c][hh])
                h = len(all_start_for_given_sat_sorted[c]) + 1
            else: 
                sub_sub_list_start_spec_fly_over_storm = []
                sub_sub_list_start_spec_fly_over_storm.append(all_start_for_given_sat_sorted[c][hh])
                sub_list_start_spec_fly_over_storm.append(sub_sub_list_start_spec_fly_over_storm)
                h = len(all_start_for_given_sat_sorted[c]) + 1
        h = h + 1
    list_start_spec_fly_over_storm.append(sub_list_start_spec_fly_over_storm)

## list of ends of fly over
list_stop_spec_fly_over_storm = []
all_stop_for_given_sat_sorted = []
for c in range(nb_cygnss):
    all_stop_for_given_sat = []
    sub_list_stop_spec_fly_over_storm = []
    for g in range(nb_gps):
        for e in range(len( time_stop_spec_above_storm[c][g])):
            all_stop_for_given_sat.append(time_stop_spec_above_storm[c][g][e])
    all_stop_for_given_sat_sorted.append(np.sort(all_stop_for_given_sat))
    skip_spec2 = 0
    h = 0
    while (h < len(all_stop_for_given_sat_sorted[c]) ):
        hh = h + skip_spec2
        skip_spec = 0
        if (hh < len(all_stop_for_given_sat_sorted[c])-1):
            sub_sub_list_stop_spec_fly_over_storm = []
            sub_sub_list_stop_spec_fly_over_storm.append(all_stop_for_given_sat_sorted[c][hh])
            while (  (hh+1+skip_spec < len(all_stop_for_given_sat_sorted[c]) - 1) & ( (all_stop_for_given_sat_sorted[c][hh+1+skip_spec]-all_stop_for_given_sat_sorted[c][hh+skip_spec])  < 20*60/step_visualization ) ) : 
# if two specular points fly over storm one after each other by less than 20 minutes then we put them in the same sub-list. If > 20 minutes, then it means that the second specular point belongs to the next orbit, so we put it in a different sub-list. 20 minutes is an arbitrary choice, it could be any number smaller than the orbit of the cygnss satellite
                sub_sub_list_stop_spec_fly_over_storm.append(all_stop_for_given_sat_sorted[c][hh+1+skip_spec])
                skip_spec = skip_spec + 1
                skip_spec2 = skip_spec2 + 1
            sub_list_stop_spec_fly_over_storm.append(sub_sub_list_stop_spec_fly_over_storm)
        if (hh == len(all_stop_for_given_sat_sorted[c])-1):
            if ( all_stop_for_given_sat_sorted[c][hh] - all_stop_for_given_sat_sorted[c][hh-1] < 20*60/step_visualization ):
                sub_list_stop_spec_fly_over_storm[-1].append(all_stop_for_given_sat_sorted[c][hh])
                h = len(all_stop_for_given_sat_sorted[c]) + 1
            else: 
                sub_sub_list_stop_spec_fly_over_storm = []
                sub_sub_list_stop_spec_fly_over_storm.append(all_stop_for_given_sat_sorted[c][hh])
                sub_list_stop_spec_fly_over_storm.append(sub_sub_list_stop_spec_fly_over_storm)
                h = len(all_stop_for_given_sat_sorted[c]) + 1
        h = h + 1
    list_stop_spec_fly_over_storm.append(sub_list_stop_spec_fly_over_storm)

# if a spec flies over the strom for 1 second then we ignore it cause it makes think crash
for c in range(nb_cygnss): 
    nb_fly_over_during_simu = len(list_start_spec_fly_over_storm[c])
    for f in range(nb_fly_over_during_simu):
        if list_start_spec_fly_over_storm[c][f][0] == list_stop_spec_fly_over_storm[c][f][-1]:
            list_start_spec_fly_over_storm[c].remove(list_start_spec_fly_over_storm[c][f])
            list_stop_spec_fly_over_storm[c].remove(list_stop_spec_fly_over_storm[c][f])


## IF YOU WANT TO PLOT THE COVERAGE OF ALL SATELLITES THEN UNCOMMENT BELOW
# fig_coverage = plt.figure()#plt.figure(figsize=(10, 8))
# fig_coverage.suptitle('COVERAGE FOR CYGNSS 5, 6, 7, AND 8 FROM '+time_spec_temp[0]+' TO '+time_spec_temp[-1],fontsize = 15, weight = 'bold', color = 'w')
# background_color = (0/555,76./255,153/255.)
# fig_coverage.set_facecolor(background_color)
# list_of_axes = []
# min_yplot = 0.6
# for c in range(nb_cygnss): # !!!!!!!!!! nb_cygnss
#     sub_list_of_axes = []# this is the sub-list containing each time there is a fly over of specular points within a same orbit    
#     nb_fly_over_during_simu = len(list_start_spec_fly_over_storm[c])
#     horizontal_space = 0.022
#     vertical_space = 0.05
#     width_each_axis = (1-horizontal_space*(nb_fly_over_during_simu+1))/nb_fly_over_during_simu
#     height_each_axis = (1-vertical_space*5)/4
#     for f in range(nb_fly_over_during_simu):
#             new_ax = plt.axes([horizontal_space + f * ( width_each_axis + horizontal_space), vertical_space + c * ( height_each_axis + vertical_space), width_each_axis, height_each_axis])
#             x_plot = range(list_start_spec_fly_over_storm[c][f][0], list_stop_spec_fly_over_storm[c][f][-1], step_visualization) 
#             new_ax.text(( x_plot[0]+x_plot[-1] ) / 2.0, min_yplot + (1-min_yplot)/60, time_spec_temp[list_start_spec_fly_over_storm[c][f][0]][5:19],horizontalalignment ='center', weight = 'bold')
#             how_long_spec_seconds_temp = list_stop_spec_fly_over_storm[c][f][-1] - list_start_spec_fly_over_storm[c][f][0]
#             how_long_spec_min = (int)(how_long_spec_seconds_temp / 60.)
#             how_long_spec_seconds = how_long_spec_seconds_temp % 60
#             how_long_spec = str(how_long_spec_min)+"'"+str(how_long_spec_seconds)+"''"
#             new_ax.text(( x_plot[0]+x_plot[-1] ) / 2.0, min_yplot - (1-min_yplot)/15, how_long_spec,horizontalalignment ='center', weight = 'bold', color = 'w')
#             new_ax.set_xlim([x_plot[0], x_plot[-1]])
#             new_ax.set_ylim([min_yplot,1])
#             new_ax.xaxis.set_ticklabels([])
#             if f > 0:
#                 new_ax.yaxis.set_ticklabels([])
#             for g in range(nb_gps):
#                 color_red = color_array[g,0]/255L  
#                 color_green = color_array[g,1]/255L  
#                 color_blue = color_array[g,2]/255L          
#                 y_plot = prob_spec[list_start_spec_fly_over_storm[c][f][0] : list_stop_spec_fly_over_storm[c][f][-1], c, g]
#                 new_ax.scatter(x_plot, y_plot, color = (color_red, color_green,color_blue),marker = '.' )
#             sub_list_of_axes.append(new_ax)
# #    list_of_axes.append(sub_list_of_axes)
# #fig_coverage.savefig('coverage_specular_joaquin_cygnss_1_to_4_test.png', facecolor=fig.get_facecolor(), edgecolor='none')
# raise Exception
## END OF IF YOU WANT TO PLOT THE COVERAGE OF ALL SATELLITES THEN UNCOMMENT BELOW

list_of_axes = []
min_yplot = 0.65
sat_coverage_choice = [[1,2]] # !!! each line has to have the same numbre of columns (select the same number of flights over for each satellite)
#                       [2,4]]

for i in range(len(sat_coverage_choice)):
    for j in range(len(sat_coverage_choice[0])):
        sat_coverage_choice[i][j] = sat_coverage_choice[i][j]-1

count_vertical = 0
for cc in range(len(sat_coverage_choice)): # !!!!!!!!!! nb_cygnss
    c = sat_coverage_choice[cc][0]
    sub_list_of_axes = []# this is the sub-list containing each time there is a fly over of specular points within a same orbit    
    nb_fly_over_during_simu = len(sat_coverage_choice[0])-1
    horizontal_space = 0.035
    vertical_space = 0.05
    width_each_axis = (1-horizontal_space*(nb_fly_over_during_simu+1))/nb_fly_over_during_simu/2
    height_each_axis = (1-vertical_space*(len(sat_coverage_choice)+1))/len(sat_coverage_choice)
    count_horizontal = 0
    for ff in range(1,len(sat_coverage_choice[0])):
        f = sat_coverage_choice[cc][ff]
        new_ax = plt.axes([horizontal_space + count_horizontal * ( width_each_axis + horizontal_space), vertical_space + count_vertical * ( height_each_axis + vertical_space), width_each_axis, height_each_axis])
#        new_ax = plt.axes([0.25,0.25 , width_each_axis, height_each_axis])
        x_plot = range(list_start_spec_fly_over_storm[c][f][0], list_stop_spec_fly_over_storm[c][f][-1], step_visualization) 
        how_long_spec_seconds_temp = list_stop_spec_fly_over_storm[c][f][-1] - list_start_spec_fly_over_storm[c][f][0]
        how_long_spec_min = (int)(how_long_spec_seconds_temp / 60.)
        how_long_spec_seconds = how_long_spec_seconds_temp % 60
        how_long_spec = str(how_long_spec_min)+"'"+str(how_long_spec_seconds)+"''"
        new_ax.text(x_plot[0], min_yplot - (1-min_yplot)/15, time_spec_temp[list_start_spec_fly_over_storm[c][f][0]][14:19], weight = 'bold', color = 'w')
        new_ax.text(x_plot[-1], min_yplot - (1-min_yplot)/15, time_spec_temp[list_stop_spec_fly_over_storm[c][f][-1]][14:19], horizontalalignment = 'right',weight = 'bold', color = 'w')
        new_ax.set_xlim([x_plot[0], x_plot[-1]])
        new_ax.set_ylim([min_yplot,1])
        new_ax.xaxis.set_ticklabels([])        
        if cc == 0:
            new_ax.text(x_plot[0] + (x_plot[-1]-x_plot[0])/2., min_yplot - (1-min_yplot)/15, 'Real Time (mm:ss)', horizontalalignment = 'center',weight = 'bold', color = 'w')
        if count_horizontal > 0:
            new_ax.yaxis.set_ticklabels([])
        for g in range(nb_gps):
            color_red = color_array[g,0]/255L  
            color_green = color_array[g,1]/255L  
            color_blue = color_array[g,2]/255L          
            y_plot = prob_spec[list_start_spec_fly_over_storm[c][f][0] : list_stop_spec_fly_over_storm[c][f][-1], c, g]
            new_ax.scatter(x_plot, y_plot, color = (color_red, color_green,color_blue),marker = '.' )
        new_ax.yaxis.label.set_color('w')
        sub_list_of_axes.append(new_ax)
        count_horizontal = count_horizontal + 1
    list_of_axes.append(sub_list_of_axes)
    count_vertical = count_vertical +1 
#fig_coverage.savefig('coverage_specular_joaquin_cygnss_1_to_4_test.png', facecolor=fig.get_facecolor(), edgecolor='none')


# ## ALL TIMES
# # fig = plt.figure()
# # plt.title('Joaquin seen by the CYGNSS constellation')
# coverage_start = when_start_visu
# coverage_stop = when_stop_visu#nb_lines
# coverage_start = (int)( coverage_start / step_visualization )
# coverage_stop = (int)( coverage_stop / step_visualization )
# x_plot = range(coverage_start,coverage_stop)
# for c in range(nb_cygnss): #range(nb_cygnss): 

#     for s in range(nb_gps):
# #        s = 4
#         color_red = color_array[s,0]/255L  
#         color_green = color_array[s,1]/255L  
#         color_blue = color_array[s,2]/255L          
#         ax2.scatter(x_plot,prob_spec[ coverage_start:coverage_stop,c, s] +c, color = (color_red, color_green,color_blue),marker = '.' ) 
#         print max(prob_spec[ coverage_start:coverage_stop,c, s] * + c)
# #         s = 12
# #         color_red = color_array[s,0]/255L  
# #         color_green = color_array[s,1]/255L  
# #         color_blue = color_array[s,2]/255L          
# #         plt.scatter(x_plot,prob_spec[ coverage_start:coverage_stop,c, s] + c, color = (color_red, color_green,color_blue),marker = '.' ) 
# #         print max(prob_spec[ coverage_start:coverage_stop,c, s] + c)


# # #  3 spec: 52658:52711 # specular_point_see_storm[0,34581:34665,0],'k.'
# ax2.axis([min(x_plot), max(x_plot), 0.6, 1]) 
# # plt.xlabel('Time (seconds)')
# # plt.ylabel('CYGNSS # - specular point probabilty')
# # plt.show()


# #which_specular_gain_not_zero_and_above_storm[coverage_start:coverage_stop,c,s]
# #fig.savefig('coverage_specular_joaquin.png')

# #end of PLOT flights over

# # ######################################################
# # ######################################################

# Animation
if ( zoom == 0 ):
    pos_time_text_ax0_x = max_lon[zoom] - ( max_lon[zoom] - min_lon[zoom] ) / 7.
    pos_time_text_ax0_y = min_lat[zoom] + ( max_lat[zoom] -  min_lat[zoom] ) / 40.
else :
    pos_time_text_ax0_x = max_lon[zoom] - ( max_lon[zoom] - min_lon[zoom] ) / 6.0
    pos_time_text_ax0_y = min_lat[zoom] + ( max_lat[zoom] -  min_lat[zoom] ) / 40.

time_text_ax0=ax0.text(pos_time_text_ax0_x, pos_time_text_ax0_y,'',horizontalalignment ='center', weight = 'bold',fontsize = 15, color = 'b')

pos_time_text_ax_coverage_x = 0
pos_time_text_ax_coverage_y = 0
# ax_coverage = list_of_axes[0][0]
#time_text_ax_coverage=list_of_axes[0][0].text(pos_time_text_ax_coverage_x, pos_time_text_ax_coverage_y,'',horizontalalignment ='center', weight = 'bold')
# time_now, = ax_coverage.plot([], [], lw=2)
list_time_now = []
for cc in range(len(sat_coverage_choice)): # !!!!!!!!!! nb_cygnss
    sub_list_time_now = []
    for ff in range(len(sat_coverage_choice[0])-1):
        time_now, = list_of_axes[cc][ff].plot([], [], lw=2)
#                    list_of_axes[cc][ff].plot([], [], lw=2)[0].set_data([], [])
        sub_list_time_now.append(time_now)
    list_time_now.append(sub_list_time_now)

def init():
	tuple_spacecraft_point_to_plot = (spacecraft_list[0].point_plot,)
        tuple_specular_point_to_plot = ()
        if (show_specular_points == 1):
            tuple_specular_point_to_plot = (specular_list[0].point_plot,)
	for i in range(nb_cygnss):
		spacecraft_list[i].point_plot.set_data([], [])
		tuple_spacecraft_point_to_plot = tuple_spacecraft_point_to_plot + (spacecraft_list[i].point_plot,)
                if (show_specular_points == 1):
                    for j in range(nb_gps):
                        specular_list[j+i*nb_gps].point_plot.set_data([], [])
                        tuple_specular_point_to_plot = tuple_specular_point_to_plot + (specular_list[j+i*nb_gps].point_plot,)

#	tuple_storm_point_to_plot = (storm_list[0].point_circle_plot,)
        tuple_storm_point_to_plot = ()
        if (show_storm == 1):
            tuple_storm_point_to_plot = (storm_list[0].point_plot,) 
            for i in range(nb_storms):
		storm_list[i].point_plot.set_data([], []) 
#		tuple_storm_point_to_plot = tuple_storm_point_to_plot + (storm_list[i].point_circle_plot,)
		tuple_storm_point_to_plot = tuple_storm_point_to_plot + (storm_list[i].point_plot,) 

        time_text_ax0.set_text('')
#        time_text_ax_coverage.set_text('')
        # vertical line to represent the current time on the coverage graph
#        tuple_time_now = (list_of_axes[0][0].plot([], [], lw=2)[0],)#[0]
        tuple_time_now = ()
        for cc in range(len(sat_coverage_choice)): # !!!!!!!!!! nb_cygnss
            sub_list_time_now = []
            for ff in range(len(sat_coverage_choice[0])-1):
                list_time_now[cc][ff].set_data([], [])
                tuple_time_now = tuple_time_now + (list_time_now[cc][ff],)

	tuple_spacecraft_and_storm_and_specular = tuple_spacecraft_point_to_plot + tuple_storm_point_to_plot + tuple_specular_point_to_plot + (time_text_ax0,) + tuple_time_now   # + (time_text_ax_coverage,) #+ tuple_time_now 
#+ (time_now,)


	return tuple_spacecraft_and_storm_and_specular


time_template = 'time = %.1fs'
print time_sat[when_start_visu],'-->',time_sat[when_stop_visu]
def animate(iii):
#       	global iii
        global tuple_spacecraft_point_to_plot
        global tuple_storm_point_to_plot
        global tuple_specular_point_to_plot
        global list_of_tuple_storm_point_to_plot
        global tuple_time_now
        global time_now
	if (iii < when_stop_visu - (when_start_visu)): # useless, should be removed
                tuple_specular_point_to_plot = ()
                tuple_storm_point_to_plot = ()
                tuple_time_now = ()
                list_of_tuple_storm_point_to_plot = []
                list_of_tuple_specular_point_to_plot = []
                color_storm = 'b'

                for cc in range(len(sat_coverage_choice)): # !!!!!!!!!! nb_cygnss
                    for ff in range(len(sat_coverage_choice[0])-1):
                        if (who_started_last == 'storms'):
                            list_time_now[cc][ff].set_data([iii + when_start_visu+nb_time_steps_between_starts, iii + when_start_visu+nb_time_steps_between_starts], [0,8])
                        else:
                            time_now, = list_of_axes[cc][ff].plot([], [], lw=2)
                            time_now.set_data([iii + when_start_visu+nb_time_steps_between_starts, iii + when_start_visu+nb_time_steps_between_starts], [0,8])

                        tuple_time_now = tuple_time_now + (list_time_now[cc][ff],)

                # CYGNSS 
		tuple_spacecraft_point_to_plot = (spacecraft_list[0].point_plot,)
		for i in range(nb_cygnss):
                        # the if and else you see in the animation are always just because I had to deal with the fact that the satellites files and the storms files could possibly not start at the same time. So don't really care about it, just keep in mind that if you add a line in the if, you need to add the same line in the else and change the index (take the same index each time as the one in lon[] and lat[])
                        if (who_started_last == 'storms'):
                            lon_a, lat_a = lon[iii + when_start_visu+nb_time_steps_between_starts,i], lat[iii + when_start_visu+nb_time_steps_between_starts,i]
# SPECULAR POINTS
                            if (show_specular_points == 1):
                                gain_array = gain_spec[iii+ when_start_visu+nb_time_steps_between_starts,i,:]
                                index_gain_not_zero = np.where(gain_array !=0)
                                index_gain_zero = np.where(gain_array ==0)
                                for ii in range(nb_gps):
                                    if ( (ii in index_gain_not_zero[0]) == True ):

                                        lon_spec_a, lat_spec_a = lon_spec[iii + when_start_visu+nb_time_steps_between_starts,i,ii], lat_spec[iii + when_start_visu+nb_time_steps_between_starts,i,ii]
                                        specular_list[ii + i*nb_gps].x,  specular_list[ii + i*nb_gps].y = m(lon_spec_a, lat_spec_a)
                                        specular_list[ii+i*nb_gps].point_plot = m.plot(specular_list[ii + i*nb_gps].x, specular_list[ii + i*nb_gps].y, marker='o', markersize = 15, color = 'chartreuse', fillstyle = 'none', mew = 2)[0]#fillstyle='none',
                                        list_of_tuple_specular_point_to_plot.append((specular_list[ii+i*nb_gps].point_plot,))


                                        if (show_storm == 1):
                                            if (specular_over_storm[iii + when_start_visu+nb_time_steps_between_starts,i, ii,0] == 1 ):
                                                color_storm = 'r'
                                            for sss in range(nb_storms):
                                                ell = get_ellipse_coords(a=radius_for_tissot(radius_uncertainty_storm[iii + when_start_visu+nb_time_steps_between_starts,sss]), b=radius_for_tissot(radius_uncertainty_storm[iii + when_start_visu+nb_time_steps_between_starts,sss]), x=lon_storm[iii + when_start_visu+nb_time_steps_between_starts,sss], y=lat_storm[iii + when_start_visu+nb_time_steps_between_starts,sss])
                                                storm_list[sss].point_plot = m.plot(ell[:,0], ell[:,1], color=color_storm,linewidth=3.5)[0]
                                                list_of_tuple_storm_point_to_plot.append((storm_list[sss].point_plot,))

                                    

                        else:
                            lon_a, lat_a = lon[iii + when_start_visu,i], lat[iii + when_start_visu,i]
                            if (show_specular_points == 1):
                                gain_array = gain_spec[iii + when_start_visu,i,:]
                                index_gain_not_zero = np.where(gain_array !=0)
                            
                                for ii in index_gain_not_zero[0]:
                                    if ii == index_gain_not_zero[0][0]:
                                        tuple_specular_point_to_plot = (specular_list[ii + i*nb_gps].point_plot,)
                                    lon_spec_a, lat_spec_a = lon_spec[iii + when_start_visu,i,ii], lat_spec[iii + when_start_visu,i,ii]

                                    specular_list[ii + i*nb_gps].x,  specular_list[ii + i*nb_gps].y = m(lon_spec_a, lat_spec_a)
                                    specular_list[ii+ i*nb_gps].point_plot.set_data(specular_list[ii+ i*nb_gps].x, specular_list[ii+ i*nb_gps].y)
                                    tuple_specular_point_to_plot = tuple_specular_point_to_plot + (specular_list[ii+ i*nb_gps].point_plot,)        

                            if (show_storm == 1):
                                if ( sat_over_storm[iii + when_start_visu,i,0] == 1):
                                    color_storm = 'r'
                                    tuple_storm_point_to_plot = (storm_list[0].point_plot,) 
                                    for sss in range(nb_storms):
                                        ell = get_ellipse_coords(a=radius_for_tissot(radius_uncertainty_storm[iii + when_start_visu,sss]), b=radius_for_tissot(radius_uncertainty_storm[iii + when_start_visu,sss]), x=lon_storm[iii + when_start_visu,sss], y=lat_storm[iii + when_start_visu,sss])
                                        storm_list[sss].point_plot.set_data(ell[:,0],  ell[:,1])
                                        list_of_tuple_storm_point_to_plot.append(tuple_storm_point_to_plot)



			spacecraft_list[i].x,  spacecraft_list[i].y = m(lon_a, lat_a)
			spacecraft_list[i].point_plot.set_data(spacecraft_list[i].x, spacecraft_list[i].y)
			tuple_spacecraft_point_to_plot = tuple_spacecraft_point_to_plot + (spacecraft_list[i].point_plot,)        

        if (show_storm == 1):
            for i in range(len(list_of_tuple_storm_point_to_plot)):
                tuple_storm_point_to_plot = tuple_storm_point_to_plot + list_of_tuple_storm_point_to_plot[i]

        if (show_specular_points == 1):
            for i in range(len(list_of_tuple_specular_point_to_plot)):
                tuple_specular_point_to_plot = tuple_specular_point_to_plot + list_of_tuple_specular_point_to_plot[i]

        time_text_ax0.set_text(time_sat[iii + when_start_visu+nb_time_steps_between_starts])

        tuple_spacecraft_and_storm_and_specular = tuple_spacecraft_point_to_plot + tuple_storm_point_to_plot + tuple_specular_point_to_plot + (time_text_ax0,) + tuple_time_now #+ (time_text_ax_coverage,) #+ tuple_time_now  #+ (time_now,)

        fig.savefig('./video6/specular'+str(iii), facecolor=fig.get_facecolor(), edgecolor='none')

        for i in range(nb_cygnss):
            if (show_specular_points == 1):
                gain_array = gain_spec[iii+ when_start_visu+nb_time_steps_between_starts,i,:]
                index_gain_not_zero = np.where(gain_array !=0)
                for ii in range(nb_gps):
                    if ( (ii in index_gain_not_zero[0]) == True ):
                        specular_list[ii+i*nb_gps].point_plot.remove() 
#                         if (specular_over_storm[iii + when_start_visu+nb_time_steps_between_starts,i, ii,0] == 1 ):                        
#                             storm_list[120sss].point_plot = m.plot(ell[:,0], ell[:,1], color=color_storm,linewidth=3.5)[0]
#                             list_of_tuple_storm_point_to_plot.append((storm_list[sss].point_plot,))
                    

	return tuple_spacecraft_and_storm_and_specular 


# call the animator.  blit=True means only re-draw the parts that have changed.
fps_we_want = 20L
anim = animation.FuncAnimation(fig, animate,init_func=init, 
                               frames = when_stop_visu - when_start_visu,
                               interval=1000L/fps_we_want, blit=True, repeat = False)#1000L/fps_we_want
#anim.save('the_movie.mp4', writer = 'ffmpeg', codec="libx264", fps = fps_we_want)


mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
plt.show()

# faire pareil quand le satellite commence avant le storm
# peut etre failt il faire un else de si on n est pas au dessus du storm alors m.plot(rien) (pareil avec les specular)





