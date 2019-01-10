import matplotlib as mpl
mpl.use('WXAgg')
from matplotlib.figure import Figure
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
from datetime import timedelta
from itertools import groupby
from operator import itemgetter
from scipy.interpolate import spline
from math import *

from scipy.interpolate import *

from bisect import *

import urllib
import gdal
from gdalconst import *

import matplotlib.gridspec as gridspec
from get_ellipse_coords import *
from radius_for_tissot import *
from formatMapString import *

from read_input_file import *

from forecast import Forecast

#kenny imports

from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

try:
    import wx
except ImportError:
    raise ImportError, "The wxPython module is required to run this application."


class SatelliteManager(object):
    curSample = 0
    curTracker = []
    speedFactor = 1 #Integer that determines speed of animation
    firstTime = ''
    lastTime = ''
    satColors = ['lime', 'blue', 'indigo', 'purple', 'dodgerblue', 'steelblue', 'seagreen', 'limegreen']
    colorByIntensity = {'DB':'skyblue', 'TD':'aqua', 'TS':'orange', 'TY':'darkcyan', 'ST':'darkslategray', 'TC':'salmon',
                                    'HU':'red', 'SD':'lime', 'SS':'g', 'EX':'purple', 'PT':'limegreen',
                                    'IN':'indigo', 'DS':'darkorchid', 'LO': 'slateblue', 'WV': 'mediumslateblue', 'ET':'plum',
                                    'XX':'slategray'}
    #This dictionary is meant to hold a 'lighter' version of the main colors representing each storm type for forecasts in future
    colorByIntensity_faded = {'DB':'skyblue', 'TD':'paleturquoise', 'TS':'navajowhite', 'TY':'darkcyan', 'ST':'darkslategray', 'TC':'salmon',
                                    'HU':'pink', 'SD':'lime', 'SS':'g', 'EX':'purple', 'PT':'limegreen',
                                    'IN':'indigo', 'DS':'darkorchid', 'LO': 'slateblue', 'WV': 'mediumslateblue', 'ET':'plum',
                                    'XX':'slategray'}


    #User dictated file locations

    driver_filepath = 'C:/Users/kcarlsen/Desktop/EasyFly/input/main_input/CYGNSS_9_30_16/cygnss_9_30_16.txt'
    storm_dir = 'C:/Users/kcarlsen/Desktop/EasyFly/input/storm_forecasts'
    gsReportLoc = 'C:/Users/kcarlsen/Desktop/EasyFly/cygnss_gs_reports'
    specReportLoc = 'C:/Users/kcarlsen/Desktop/EasyFly/cygnss_spec_reports'
    sat_filepath = "C:/Users/kcarlsen/Desktop/EasyFly/input/main_input/CYGNSS_9_30_16/interpolated_position_LLA_ECEF_"
    specular_filepath = "C:/Users/kcarlsen/Desktop/EasyFly/input/main_input/CYGNSS_9_30_16/coverage_specular_"
    gsInteraction_filepath = "C:/Users/kcarlsen/Desktop/EasyFly/input/main_input/CYGNSS_9_30_16/report_coverage_all_stations_by_"
    specInteraction_filepath = "C:/Users/kcarlsen/Desktop/EasyFly/input/main_input/CYGNSS_9_30_16/coverage_specular_"

    def __init__(self, m, maxlonIn, minlonIn, latIn, heightIn, widthIn, axMapIn, inZoom):

        self.step_visualization = 1

        self.nb_storms = 0

        self.show_graph = 0 #you can't show graph if you don't plot coverage. but you can plot coverage and not show graph

        self.stormsShowed = False

        self.inZoom = inZoom

        #These parameters are used only for placing clock on map, copied from map instantiation
        self.max_lon = maxlonIn
        self.min_lon = minlonIn
        self.min_lat = latIn
        self.height_fig_with_map = heightIn
        self.width_fig_with_map = widthIn
        self.ax_map = axMapIn

        ###Accumulate satellite data###
        self.satNames = ['SAT1', 'SAT2', 'SAT3', 'SAT4', 'SAT5', 'SAT6', 'SAT7', 'SAT8']
        self.specNames = ['SPEC1', 'SPEC2', 'SPEC3', 'SPEC4']

        #Read propagator input files
        self.input_variables, self.order_input_variables = read_input_file(self.driver_filepath)
        self.date_start = self.input_variables[0]
        self.date_stop = self.input_variables[1]
        self.dt = self.input_variables[2]
        self.nb_steps = self.input_variables[3]
        self.nb_satellites = self.input_variables[4]
        self.gps_name = self.input_variables[5]
        self.output_file_path_propagator = self.input_variables[6]
        self.output_filename_propagator = self.input_variables[7]
        #self.nb_storms = self.input_variables[9]
        self.storm_name = self.input_variables[10]

        #Read propagator output files
        self.nb_spec_pts = 4
        self.color_spec = ['c', 'r','b', 'k','g','m','y','w']
        self.interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
        self.nb_steps_interpolation = (int)((self.nb_steps-1) * self.dt / self.interpolation_step) +1

        #Declare necessary variables
        self.lon_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
        self.lat_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
        self.ecef_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation, 3])
        self.time_sat = []
        self.lon_spec, self.lat_spec, self.gain_spec, self.distance_spec_to_storm, self.specular_over_storm = np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation])
        self.ecef_spec = np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation, 3])
        self.name_spec = []
        self.time_spec = []
        self.date = []
        self.list_output_variables_to_read = ["longitude","latitude"]
        self.point = namedtuple('point', ['x', 'y'])
        self.color = namedtuple('color', 'red green blue')
        self.spacecraft_list = []
        self.specular_list = []
        self.init_index = 3 #27, 5, WHO KNOWS

        #Build out satellite data

        #SAT TRACE LINES - used to show ground tracks
        self.satTraceX = [[] for i in range(self.nb_satellites)]
        self.satTraceY = [[] for i in range(self.nb_satellites)]
        self.specTraceX = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.specTraceY = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.satTraceRefs = {}
        self.specTraceRefs = {}

        #Ground station interaction parameters
        self.GSinteractionIdx = None


        

        for i in range(self.nb_satellites):
            self.time_spec_sublist = []
            self.name_spec_between_list_and_sublist = []

            ### SATELLITES ###
            output_file = open(self.sat_filepath + self.output_filename_propagator[i], "r")
            read_output_file = output_file.readlines()
            nb_lines_header_output_file_sat = 0
            
            for j in range(self.nb_steps_interpolation):
                self.lon_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[1])
                if (self.lon_sat[i,j] > 180):
                    self.lon_sat[i,j] = self.lon_sat[i,j] - 360.
                self.lat_sat[i,j] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[2])
                self.ecef_sat[i,j,0] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[4])
                self.ecef_sat[i,j,1] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[5])
                self.ecef_sat[i,j,2] = np.float(read_output_file[j+nb_lines_header_output_file_sat].split()[6])
                self.time_sat.append(read_output_file[j+nb_lines_header_output_file_sat].split()[0])
                

            # Build the tuples for the visualization of the satellites
            spacecraft = namedtuple('spacecraft',('name',) +  self.point._fields + ('point_plot',) + ('marker_spacecraft',))
            self.spacecraft_list.append(spacecraft)
            name_temp = self.output_filename_propagator[i].replace(".txt","")
            self.spacecraft_list[i].name = name_temp
            # initial position
            self.spacecraft_list[i].x, self.spacecraft_list[i].y =  m(self.lon_sat[i,self.init_index], self.lat_sat[i,self.init_index])
            self.spacecraft_list[i].marker_spacecraft = 'v'
            # point on the plot
            self.spacecraft_list[i].point_plot = m.plot([],[],  marker=self.spacecraft_list[i].marker_spacecraft, markersize=10,color = self.satColors[i])[0]

            ###Add in specular points
            self.spec_dir = ""
            for j in range(len(self.output_file_path_propagator[i].split('/'))-2):
                if (j > 0):
                    self.spec_dir = self.spec_dir + "/" + self.output_file_path_propagator[i].split('/')[j]

            file_specular = open(self.specular_filepath + self.output_filename_propagator[i], "r")
            #file_specular = open(self.spec_dir + "/coverage/storm/coverage_specular_" + self.output_filename_propagator[i], "r")
            read_file_specular  = file_specular.readlines()
            # Nb of lines in the spec file header
            if (i == 0):
                nb_lines_header_output_file_spec = 0
                while (read_file_specular[nb_lines_header_output_file_spec].split()[0] != "#START"):
                    nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
                nb_lines_header_output_file_spec = nb_lines_header_output_file_spec + 1
            ispec_save = 0
            ## the output of find_specular_points.c does not start at the same time as the propagation 
            j = 0
            while ( datetime.strptime(read_file_specular[nb_lines_header_output_file_spec].split()[0], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(self.time_sat[j], "%Y-%m-%dT%H:%M:%S") ):
                j = j +1
                self.time_spec_sublist.append('')
                self.name_spec_between_list_and_sublist.append('')
            j = j-1
            while (ispec_save < len(read_file_specular)-1-self.nb_spec_pts):
                j = j + 1
                self.time_spec_sublist_temp_ini = read_file_specular[ispec_save+nb_lines_header_output_file_spec].split()[0] 
                self.time_spec_sublist_temp_ini = datetime.strptime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
                self.time_spec_sublist.append(datetime.strftime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
                self.name_spec_sublist = []
                self.lon_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[1]) #CHANGED
                if self.lon_spec[0,i,j] > 180:
                    self.lon_spec[0,i,j] = self.lon_spec[0,i,j] - 360.
                self.lat_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[2]) #CHANGED
                self.gain_spec[0,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[3]) #CHANGED
                self.name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save].split()[4]) #CHANGED
                ispec = 1
                while (datetime.strptime(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == self.time_spec_sublist_temp_ini):
                    self.lon_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
                    if self.lon_spec[ispec,i,j] > 180:
                        self.lon_spec[ispec,i,j] = self.lon_spec[ispec,i,j] - 360.
                    self.lat_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
                    self.gain_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
                    self.name_spec_sublist.append(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[4])
                    ispec = ispec + 1
                    if (nb_lines_header_output_file_spec+ispec_save+ispec == len(read_file_specular) - 1):
                        break
                ispec_save = ispec + ispec_save
                self.name_spec_between_list_and_sublist.append(self.name_spec_sublist)
            # if up to here we still ahve read the entire spec file
            if ( datetime.strptime(self.time_spec_sublist[-1], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(read_file_specular[len(read_file_specular)-2].split()[0], "%Y-%m-%dT%H:%M:%S") ):
                j = j + 1
                first_spec_of_last_time_step = +ispec_save
                self.time_spec_sublist_temp_ini = read_file_specular[first_spec_of_last_time_step].split()[0]
                self.time_spec_sublist_temp_ini = datetime.strptime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S")
                self.time_spec_sublist.append(datetime.strftime(self.time_spec_sublist_temp_ini, "%Y-%m-%dT%H:%M:%S"))
                self.name_spec_sublist = []
                self.lon_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[1])
                if self.lon_spec[0,i,j] > 180:
                    self.lon_spec[0,i,j] = self.lon_spec[0,i,j] - 360.
                self.lat_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[2])
                self.gain_spec[0,i,j] = np.float(read_file_specular[first_spec_of_last_time_step].split()[3])
                ispec = 1
                while (datetime.strptime(read_file_specular[first_spec_of_last_time_step+ispec].split()[0], "%Y-%m-%dT%H:%M:%S")  == self.time_spec_sublist_temp_ini):
                    self.lon_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[1])
                    if self.lon_spec[ispec,i,j] > 180:
                        self.lon_spec[ispec,i,j] = self.lon_spec[ispec,i,j] - 360.
                    self.lat_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[2])
                    self.gain_spec[ispec,i,j] = np.float(read_file_specular[nb_lines_header_output_file_spec+ispec_save+ispec].split()[3])
                    if (first_spec_of_last_time_step+ispec < len(read_file_specular) - 1):
                        ispec = ispec + 1
                    else: 
                        break
                self.name_spec_between_list_and_sublist.append(self.name_spec_sublist)

            ## the output of find_specular_points.c does not end at the same time as the propagation 
            j_end = 0
            while ( datetime.strptime(self.time_spec_sublist[-1-j_end], "%Y-%m-%dT%H:%M:%S") != datetime.strptime(self.time_sat[-1-j_end], "%Y-%m-%dT%H:%M:%S") ):
                j_end = j_end +1
                self.time_spec_sublist.append('')
                self.name_spec_between_list_and_sublist.append('')
            self.time_spec.append(self.time_spec_sublist)
        # Build the tuples for the visualization of the specular points
            for k in range(self.nb_spec_pts):
                self.specular = namedtuple('specular',('name',) +  self.point._fields  + self.color._fields + ('point_plot',))
                self.specular_list.append(self.specular)
            # name 
            #       self.specular_list[i+j*nb_gps].name = self.specular_names[i]
            # initial position                    
                self.specular_list[k+i*self.nb_spec_pts].x, self.specular_list[k+i*self.nb_spec_pts].y =  m(self.lon_spec[k,i,self.init_index], self.lat_spec[k,i,self.init_index])
                self.specular_list[k+i*self.nb_spec_pts].point_plot = m.plot([],[], marker='o', markersize = 2, color = self.satColors[i], fillstyle = 'none', mew = 2)[0]

            self.name_spec.append(self.name_spec_between_list_and_sublist)


        ###Calculate Default Ground traces###
        self.calcTraces(0, self.nb_steps_interpolation)
        ###Build Ground Stations
        self.gs_name = ['Hawaii', 'Chile', 'Austrailia'] #Note: Order is important for indexing purposes
        self.gs_lat = [21, -33, -31]
        self.gs_lon = [-157, -70, 115]
        self.gsRefs = {}

        self.nb_gs = len(self.gs_name)

        ###Build Storms

        self.forecastTimes = []
        self.forecastLats = []
        self.forecastLons = []
        self.forecastColors = []
        self.forecastRadii = []
        
        self.forecastTimeIdxs = []
        self.forecastTimestamps = []
        self.mapTextRefs = [] #For use in plot storms, probably shouldn't be declared here
        self.activeStorms = []

        self.build_raw_storms()
        self.generate_interpolators()
        self.interpolate_storms()

        self.forecastArtistlist = [[] for x in range(self.nb_storms)]
        self.trackerArtistlist = [[] for x in range(self.nb_storms)]
        self.trajArtistlist = []

        ###Build Map Time display###
        self.mapTextCorrectionFactor = 2
        #CurTime clock + Display start stop location parameters
        self.pos_display_text_ax_map_x = self.max_lon[self.inZoom] - self.mapTextCorrectionFactor
        self.pos_display_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.display_text_ax_map = self.ax_map.text(self.pos_display_text_ax_map_x, self.pos_display_text_ax_map_y,'',horizontalalignment ='right', weight = 'bold',fontsize = 15, color = 'r')
        #Warning time (forecast time) display
        self.pos_warning_text_ax_map_x = self.min_lon[self.inZoom] + self.mapTextCorrectionFactor
        self.pos_warning_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.time_warning_text_ax_map = self.ax_map.text(self.pos_warning_text_ax_map_x, self.pos_warning_text_ax_map_y,'',horizontalalignment ='left', weight = 'bold',fontsize = 15, color = 'r')
        

        ###Build Sat/GS interaction data structures###
        self.calcGSinteraction()

        ###Build Spec/Storm interaction data structures###
        #self.findSpecStormPasses()

        ###Call initializing function###
        self.init_data(m)

    def build_raw_storms(self):
        #This list holds tuples of raw forecast data
        self.rawForecasts = []
        #This list holds tuples of (stormName, warningTime)
        self.stormInfo = []

        #Need start/end times when deciding indices
        startTime = datetime.strptime(self.time_sat[0], "%Y-%m-%dT%H:%M:%S")
        endTime = datetime.strptime(self.time_sat[-1], "%Y-%m-%dT%H:%M:%S")

        #Get all forecast filenames in ../input/storm_forecasts directory
        stormFiles = [f for f in os.listdir(self.storm_dir) if os.path.isfile(os.path.join(self.storm_dir, f))]

        for file in range(len(stormFiles)):
            STORM_FILEPATH = self.storm_dir + '/' + stormFiles[file]
            stormFile = open(STORM_FILEPATH, 'r')
            print self.storm_dir + '/' + stormFiles[file]
            forecastLines = stormFile.readlines()

            #tmpForecasts allows us to append on a per-storm basis
            tmpForecasts = []
            #tmpTimes allows us to only look at 34 kt radii for each storm
            tmpTimes = []
            #tmpLons lets us address duplicate longitudes
            tmpLons = []
            if(STORM_FILEPATH[-3:] == 'txt'):
                #Need first and last cast to determine if storm is in interval
                firstCast = [x.strip() for x in forecastLines[0].split(',')]
                firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H") + timedelta(hours=int(lastCast[3]))

                #Check if storm falls at all in display interval
                if firstCastTime >= startTime or lastCastTime >= startTime:
                    #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                    self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                    #For each line in forecast file
                    for iStorm in range(len(forecastLines)):
                        noaa = [x.strip() for x in forecastLines[iStorm].split(',')]
                        #Build basic storm parameters
                        if noaa[3] not in tmpTimes: #So we only read forecast one time 
                            #Timestamp
                            ts = datetime.strptime(noaa[2], "%Y%m%d%H") + timedelta(hours=int(noaa[3])) 

                            #Index
                            if ts < startTime:
                                sIdx = -2 #Arbitrary flag 
                            elif ts > endTime:
                                sIdx = -1 #Arbitrary flag
                            else:
                                sIdx = int(self.time_sat.index(ts.isoformat()))

                            #Lat
                            lat = float(noaa[4])

                            #Lon
                            lon = float(noaa[5])
                            while lon in tmpLons:
                                print "Found lon in tmp lons, incrementing."
                                lon -= 0.3
                            tmpLons.append(lon)

                            #Storm type
                            typ = noaa[9]

                            #Size
                            windRadii = [int(noaa[12]), int(noaa[13]), int(noaa[14]), int(noaa[15])]
                            size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding
                            
                            tmpTimes.append(noaa[3]) 

                            newForecast = (ts, sIdx, lat, lon, typ, size)
                            
                            tmpForecasts.append(newForecast)
                    if tmpForecasts:
                        #Create list of lists of tuples (one list for each active storm)
                        self.rawForecasts.append(tmpForecasts)
            
            elif(STORM_FILEPATH[-3:] == 'fst'):
                print("IN ELIF")
                firstCast = [x.strip() for x in forecastLines[0].split(',')]
                firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H") + timedelta(hours=int(lastCast[5]))

                #Check if storm falls at all in display interval
                if firstCastTime >= startTime or lastCastTime >= startTime:
                    #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                    self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                    #For each line in forecast file
                    for iStorm in range(len(forecastLines)):
                        jtwc = [x.strip() for x in forecastLines[iStorm].split(',')]
                        #Build basic storm parameters
                        if jtwc[5] not in tmpTimes: #So we only read forecast one time 
                            #Timestamp
                            ts = datetime.strptime(jtwc[2], "%Y%m%d%H") + timedelta(hours=int(jtwc[5])) 

                            #Index
                            if ts < startTime:
                                sIdx = -2 #Arbitrary flag 
                            elif ts > endTime:
                                sIdx = -1 #Arbitrary flag
                            else:
                                sIdx = int(self.time_sat.index(ts.isoformat()))

                            #Lat
                            lat = jtwc[6]
                            if lat[-1:] == 'N':
                                lat = float(lat[0:3])/10.0
                            else:
                                lat = -1*float(lat[0:3])/10.0
                            print lat

                            #Lon
                            lon = jtwc[7]
                            if lon[-1:] =='E':
                                lon = float(lon[0:3])/10.0
                            else:
                                lon = -1*float(lon[0:3])/10.0

                            print lon

                            while lon in tmpLons:
                                print "Found lon in tmp lons, incrementing."
                                lon -= 0.3
                            tmpLons.append(lon)

                            #Storm type
                            typ = jtwc[10]

                            #Size
                            windRadii = [int(jtwc[13]), int(jtwc[14]), int(jtwc[15]), int(jtwc[16])]
                            size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding
                            
                            tmpTimes.append(jtwc[5]) 

                            newForecast = (ts, sIdx, lat, lon, typ, size)
                            
                            tmpForecasts.append(newForecast)
                    if tmpForecasts:
                        #Create list of lists of tuples (one list for each active storm)
                        self.rawForecasts.append(tmpForecasts)

        self.nb_storms = len(self.rawForecasts)

    #Slightly alters data to force strictly increasing requirement
    #inList must be presorted
    def fix_duplicates(self, inList):
        if len(inList) < 1:
            return inList

        cur = inList[0]
        counter = 0
        for idx in range(1, len(inList)):
            if inList[idx] == cur:
                counter += 1
                inList[idx] = inList[idx] + (0.01*counter)
            else:
                counter = 0
                cur = inList[idx]

        return inList


    def generate_interpolators(self):
        self.interpolators = []

        for system in self.rawForecasts:
            interpLats = [x[2] for x in system]
            interpLons = [x[3] for x in system]
            #Must remove any duplicate lons to satisfy spline condition
            combined = zip(interpLons, interpLats)
            #pruned = list(OrderedDict(combined[::-1]).items())[::-1]
            combined = sorted(combined, key=itemgetter(0))
            interpLons, interpLats = zip(*combined)
            interpLons = list(interpLons)
            interpLats = list(interpLats)
            

            self.interpolators.append(InterpolatedUnivariateSpline(interpLons, interpLats, ext=0))

        #print "Interpolators length: " + str(len(interpolators))

    def interpolate_storms(self):
        finalTimestamps = [[] for x in range(self.nb_storms)]
        finalIdxs = [[] for x in range(self.nb_storms)]
        finalLats = [[] for x in range(self.nb_storms)]
        finalLons = [[] for x in range(self.nb_storms)]
        finalTypes = [[] for x in range(self.nb_storms)]
        finalSizes = [[] for x in range(self.nb_storms)]
        
        timeSkip = 720 #in seconds
        startTime = datetime.strptime(self.time_sat[0], "%Y-%m-%dT%H:%M:%S")
        endTime = datetime.strptime(self.time_sat[-1], "%Y-%m-%dT%H:%M:%S")

        systemIdx = 0
        for system in self.rawForecasts:
            if len(system) > 1:
                for i in range(len(system) - 1):
                    #Figure out number of steps in this interval
                    intervals = int(divmod((system[i+1][0] - system[i][0]).total_seconds(), timeSkip)[0])
                    for j in range(intervals):
                        #Build new timestamps
                        fts = system[i][0] + j*timedelta(seconds=timeSkip)
                        finalTimestamps[systemIdx].append(fts)

                        #Retrieve new indices
                        if fts < startTime:
                            finalIdxs[systemIdx].append(-2) #Arbitrary flag 
                        elif fts > endTime:
                            finalIdxs[systemIdx].append(-1) #Arbitrary flag
                        else:
                            finalIdxs[systemIdx].append(int(self.time_sat.index(fts.isoformat())))

                    #Generate new lons
                    
                    newLons = np.linspace(system[i][3], system[i+1][3], intervals, endpoint=False)
                    newLons = [round(x, 2) for x in newLons]

                    finalLons[systemIdx].extend(newLons)

                    #Generate new types
                    finalTypes[systemIdx].extend([system[i][4]] * intervals)

                    #Generate new sizes
                    finalSizes[systemIdx].extend(np.linspace(system[i][5], system[i+1][5], intervals, endpoint=False))
            
            #Generate new lats from interpolator
            finalLats[systemIdx].extend(self.interpolators[systemIdx](finalLons[systemIdx]))

            systemIdx += 1

        #Find start/end indices for each storm
        self.stormStart = []
        self.stormEnd = []
        self.curStorm = []

        for i in range(self.nb_storms):
            self.stormStart.append(bisect(finalIdxs[i], -2)) #Find first value not before sim start
            self.curStorm.append(self.stormStart[i])
            try:
                self.stormEnd.append(finalIdxs[i].index(-1)) #Find first value past end of interval
            except ValueError:
                self.stormEnd.append(len(finalIdxs[i]) - 1)

        self.finalForecasts = []

        for i in range(self.nb_storms):
            self.finalForecasts.append(zip(finalTimestamps[i], finalIdxs[i], finalLats[i], finalLons[i], finalTypes[i], finalSizes[i]))

    ### Data/Simulation managing functions ###
    def init_data(self, m):
        #SATELLITES AND SPECULAR POINTS
        self.tuple_spacecraft_point_to_plot = ()
        self.tuple_specular_point_to_plot = ()
        for isat in range(self.nb_satellites):
            self.spacecraft_list[isat].point_plot.set_data([], []) # set_data must be applied to a line object
            self.tuple_spacecraft_point_to_plot = self.tuple_spacecraft_point_to_plot + (self.spacecraft_list[isat].point_plot,)
            for ispec in range(self.nb_spec_pts):
                self.specular_list[ispec+isat*self.nb_spec_pts].point_plot.set_data([], [])
                self.tuple_specular_point_to_plot = self.tuple_specular_point_to_plot + (self.specular_list[ispec+isat*self.nb_spec_pts].point_plot,)

        #Prepare for animation
        self.firstTime = self.time_sat[0] #First and last time maintained separately to maintain absolute bounds
        self.lastTime = self.time_sat[-1]
        #These values are strings
        self.start_visu = self.firstTime
        self.stop_visu = self.lastTime
        #These values are indices
        self.when_start_visu = self.time_sat.index(self.start_visu)
        self.when_stop_visu = self.time_sat.index(self.stop_visu)
        self.curSample = self.when_start_visu
        self.delay_show_prob_before = 0 # in seconds # the probability distributino functions of the spec on the left of the map will showp up before the spec enters 
        self.delay_show_prob_after = 0 # in seconds # the probability distributino functions of the spec on the left of the map will stay after the spec leaves
        self.delay_show_prob_before = self.delay_show_prob_before / self.interpolation_step
        self.delay_show_prob_after = self.delay_show_prob_after / self.interpolation_step
        
        print "Analysis Interval: " + self.time_sat[self.when_start_visu] + ' -> ' + self.time_sat[self.when_stop_visu]

        self.tuple_spacecraft_and_specular = self.tuple_spacecraft_point_to_plot + self.tuple_specular_point_to_plot

        #TIME
        self.display_text_ax_map.set_text('')
        self.time_warning_text_ax_map.set_text("Reported Warning Times \n ---\n" + formatWarningTime(self.stormInfo))
        self.tuple_spec_text_ax_map = ()

        #GS_INTERACTION
        self.GSinteractionIdx = [0] * self.nb_satellites #We need 8 indices to mark our progress, one for each sat

    def advSim(self, plotMap):
        self.advSimtuple_spacecraft_point_to_plot = ()
        self.advSimtuple_specular_point_to_plot = ()
        
        #This is where the sats are updated every iteration during animation
        #To not show, only run inner code block for satellites marked as 'true'
        for isat in range(self.nb_satellites):
            # SPACECRAFT
            lon_a, lat_a = self.lon_sat[isat,self.curSample], self.lat_sat[isat,self.curSample]
            self.spacecraft_list[isat].x,  self.spacecraft_list[isat].y = plotMap(lon_a, lat_a)
            self.spacecraft_list[isat].point_plot.set_data(self.spacecraft_list[isat].x, self.spacecraft_list[isat].y)       
            # SPECULAR POINTS
            kspec = 0
            # if there is less specular points than self.nb_spec_pts for for this iteration and satellite, then the old specular point stays on the image. So for loop blow is to make sure that does not happen (there are smarter ways to do that...)
            for ispec in range(self.nb_spec_pts):
                self.specular_list[ispec + isat*self.nb_spec_pts].point_plot.set_data([], []) 
            for ispec in range(len(self.name_spec[isat][self.curSample])):
                lon_spec_a, lat_spec_a = self.lon_spec[ispec,isat,self.curSample], self.lat_spec[ispec,isat,self.curSample]
                self.specular_list[ispec + isat*self.nb_spec_pts].x,  self.specular_list[ispec + isat*self.nb_spec_pts].y = plotMap(lon_spec_a, lat_spec_a)
                self.specular_list[ispec+isat*self.nb_spec_pts].point_plot.set_data(self.specular_list[ispec+isat*self.nb_spec_pts].x, self.specular_list[ispec+isat*self.nb_spec_pts].y)
                self.advSimtuple_specular_point_to_plot = self.advSimtuple_specular_point_to_plot + (self.specular_list[ispec+isat*self.nb_spec_pts].point_plot,)

        #TIME
        self.display_text_ax_map.set_text("Start: " + formatDisplayTime(self.time_sat[self.when_start_visu]) + "\n"  +  "Current: "  + formatDisplayTime(self.time_sat[self.curSample]) + "\n" + "End:   " + formatDisplayTime(self.time_sat[self.when_stop_visu]))

        #GS_INTERACTION
        for i in range(self.nb_satellites):
            #I don't love this try block, maybe reconsider - also integrate this expression into above block?
            try:
                if datetime.strptime(self.time_sat[self.curSample], '%Y-%m-%dT%H:%M:%S') > self.overGSstart[i][self.GSinteractionIdx[i]]:
                    self.GSinteractionIdx[i] += 1
                else:
                    pass

            except IndexError:
                pass

        #STORMS
        if self.stormsShowed:
            for storm in range(self.nb_storms):            
                if self.curStorm[storm] != self.stormEnd[storm] and self.curSample >= self.finalForecasts[storm][self.curStorm[storm] + 1][1]:
                    self.curStorm[storm] += 1
                    self.removeFirstStorm(storm)
                    self.updateFirstStorm(storm)
 
    def incrementIndex(self):

        if(self.curSample == self.when_stop_visu):
            self.curSample = self.when_start_visu
            self.forwardLoopReset()

        if(self.curSample + self.speedFactor < self.when_stop_visu):
            self.curSample += self.speedFactor
        else:
            self.curSample = self.when_stop_visu

    def decrementIndex(self):
        if(self.curSample == self.when_start_visu):
            self.curSample = self.when_stop_visu
            self.findGSinteractionIdx()

        if(self.curSample - self.speedFactor > self.when_start_visu):
            self.curSample -= self.speedFactor
        else:
            self.curSample = self.when_start_visu

    def updateCurStorm(self):

        for storm in range(self.nb_storms):
            print "idxs: "
            idxs = [x[1] for x in self.finalForecasts[storm]]
            print idxs

            self.curStorm[storm] = idxs.index(self.find_le(idxs, self.curSample))
            print self.curStorm[storm]

    #TODO: Create a reverseLoopReset()

    def forwardLoopReset(self):
        self.findGSinteractionIdx()
        if self.stormsShowed:
            self.removeStorms()
            self.updateCurStorm()
            self.plotStorms(self.plotMap)


    ### Map Managing functions ###
    def rebuildTimeDisplay(self, zoom):

        self.inZoom = zoom
        #Remove old text
        self.display_text_ax_map.remove()
        self.time_warning_text_ax_map.remove()

        self.pos_display_text_ax_map_x = self.max_lon[self.inZoom] - self.mapTextCorrectionFactor
        self.pos_display_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.display_text_ax_map = self.ax_map.text(self.pos_display_text_ax_map_x, self.pos_display_text_ax_map_y,'',horizontalalignment ='right', weight = 'bold',fontsize = 15, color = 'r')

        #Interval display
        self.pos_warning_text_ax_map_x = self.min_lon[self.inZoom] + self.mapTextCorrectionFactor
        self.pos_warning_text_ax_map_y = self.min_lat[self.inZoom] + ( self.height_fig_with_map - self.width_fig_with_map ) / 40. + self.mapTextCorrectionFactor
        self.time_warning_text_ax_map = self.ax_map.text(self.pos_warning_text_ax_map_x, self.pos_warning_text_ax_map_y,'',horizontalalignment ='left', weight = 'bold',fontsize = 15, color = 'r')
        #Set warning text here since it doesn't change over the course of the animation
        self.time_warning_text_ax_map.set_text("Reported Warning Times \n ---\n" + formatWarningTime(self.stormInfo))

    def plotGroundStation(self, plotMap, idx):
        #self.gsRefs[self.gs_name[idx]] = plotMap.plot(self.gs_lon[idx], self.gs_lat[idx], color='black', marker='o', markersize=40, fillstyle='none', label=self.gs_name[idx])
        self.gsRefs[self.gs_name[idx]] = plotMap.tissot(self.gs_lon[idx], self.gs_lat[idx], radius_for_tissot(2000), 256, facecolor='none', linestyle='dashdot', linewidth=3)
        
    def removeGroundStation(self, idx):
        doomedGS = self.gsRefs.pop(self.gs_name[idx]) #This pop refers to the dictionary holding map objects
        doomedGS.remove() #This pop refers to the artist managing the map object

    def find_le(self, inList, x):
        #'Find rightmost value less than or equal to x'
        i = bisect_right(inList, x)
        if i:
            return inList[i-1]
        raise ValueError

    def find_gt(self, a, x):
        #'Find leftmost value greater than x'
        i = bisect_right(a, x)
        if i != len(a):
            return a[i]
        raise ValueError

    def plotStorms(self, plotMap):
        self.activeStorms = []
        #This save and then repass in forwardloopreset seems unnecessary - improve?
        self.plotMap = plotMap
        first_in_past = True

        #Perhaps I should move the updateCurStorm call to here, instead of in forwardloopreset
        #This way, no matter when you plot the storms, curStorm will automatically be fast forwarded
        #To the appropriate location
        
        for storm in range(self.nb_storms):
            for idx in range(self.curStorm[storm], len(self.finalForecasts[storm])):
                self.forecastArtistlist[storm].append(plotMap.tissot(self.finalForecasts[storm][idx][3], self.finalForecasts[storm][idx][2], radius_for_tissot(self.finalForecasts[storm][idx][5]), 256, facecolor=self.colorByIntensity_faded[self.finalForecasts[storm][idx][4]], edgecolor='none'))
            
            self.updateFirstStorm(storm)
        self.stormsShowed = True
             

    def removeStorms(self):
        for storm in range(self.nb_storms):
            while len(self.forecastArtistlist[storm]) > 0:
                self.forecastArtistlist[storm][len(self.forecastArtistlist[storm]) - 1].remove()
                self.forecastArtistlist[storm].pop()
            #self.removeForecastText()
            # for storm in range(self.nb_storms):
            #     self.removeStormTracker(storm)
            self.stormsShowed = False

    def removeFirstStorm(self, storm):
        if len(self.forecastArtistlist[storm]) > 0:
            self.forecastArtistlist[storm][0].remove()
            self.forecastArtistlist[storm].pop(0)

    def updateFirstStorm(self, storm):
        if len(self.forecastArtistlist[storm]) > 0:
            print "In if in updateFirstStorm"
            print storm
            self.forecastArtistlist[storm][0].set_zorder(2)
            self.forecastArtistlist[storm][0].set_facecolor(self.colorByIntensity[self.finalForecasts[storm][self.curStorm[storm]][4]])

    def plotStormTracker(self, storm, idx, plotMap):
        self.removeStormTracker(storm) #Remove old tracker
        self.trackerArtistlist[storm].append(plotMap.tissot(self.newLons[storm][idx], self.newLats[storm][idx], radius_for_tissot(self.newRad[storm][idx]), 256, facecolor='none', edgecolor='black', linewidth=2))

    def removeStormTracker(self, storm):
        while len(self.trackerArtistlist[storm]) > 0:
            self.trackerArtistlist[storm][len(self.trackerArtistlist[storm]) - 1].remove()
            self.trackerArtistlist[storm].pop()

    def plotTrajectory(self, plotMap):
        for storm in range(self.nb_storms):
            stormLons = [x[3] for x in self.finalForecasts[storm]]
            stormLats = [x[2] for x in self.finalForecasts[storm]]
            self.trajArtistlist.append(plotMap.plot(stormLons, stormLats, zorder=5, color='black', linewidth=1, linestyle='dashed', alpha=0.5))

    def removeTrajectory(self):
        print self.nb_storms
        for storm in range(self.nb_storms):
            doomedTraj = self.trajArtistlist.pop(0)
            doomedTraj.pop(0).remove()

    def plotForecastText(self):
        for idx in self.activeStorms:
            self.mapTextRefs.append(self.ax_map.text(self.forecastLons[idx] + 5, self.forecastLats[idx] + 8, '', rotation=45, horizontalalignment ='right', fontsize = 8, color = 'black').set_text(self.forecastTimestamps[idx]))

    def removeForecastText(self):
        while len(self.mapTextRefs) > 0:
            self.mapTextRefs[len(self.mapTextRefs) - 1].remove(0)
            self.mapTextRefs.pop()

    def plotSatTrace(self, plotMap, satNum):
        self.satTraceRefs[self.satNames[satNum]] = plotMap.plot(self.satTraceX[satNum], self.satTraceY[satNum], 'red', linestyle='', marker='.', markersize=1)

    def removeSatTrace(self, satNum):
        doomedTrace = self.satTraceRefs.pop(self.satNames[satNum])
        doomedTrace.pop(0).remove()

    def plotSpecTrace(self, plotMap, satNum):

        self.specTraceRefs[self.satNames[satNum]] = []
        for idx in range(len(self.specNames)):
            #self.specTraceRefs[self.satNames[satNum]].append(plotMap.scatter(self.specTraceX[satNum][idx], self.specTraceY[satNum][idx], s=self.specTraceGain[satNum][idx])) This allows variable markersize but is slow
            self.specTraceRefs[self.satNames[satNum]].append(plotMap.plot(self.specTraceX[satNum][idx], self.specTraceY[satNum][idx], color=self.satColors[satNum], linestyle='', marker='.', markersize=1))

    def removeSpecTrace(self, satNum):
        doomedTraces = self.specTraceRefs[self.satNames[satNum]]
        for trace in doomedTraces:
            trace.pop(0).remove()

    def calcTraces(self, startSample, endSample):
        self.satTraceX = [[] for i in range(self.nb_satellites)]
        self.satTraceY = [[] for i in range(self.nb_satellites)]
        self.specTraceX = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.specTraceY = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]
        self.specTraceGain = [[[] for i in range(self.nb_spec_pts)] for i in range(self.nb_satellites)]

        for i in range(self.nb_satellites):
            for j in range(startSample, endSample):
                self.satTraceX[i].append(self.lon_sat[i, j])
                self.satTraceY[i].append(self.lat_sat[i, j])

            for k in range(self.nb_spec_pts):
                for step in range(startSample, endSample-1): #Not sure about this -1, but otherwise there is an extra 0
                    self.specTraceX[i][k].append(self.lon_spec[k,i,step])
                    self.specTraceY[i][k].append(self.lat_spec[k,i,step])
                    if self.gain_spec[k,i,step] < 10:
                        self.specTraceGain[i][k].append(.01)
                    elif self.gain_spec[k,i,step] < 30:
                        self.specTraceGain[i][k].append(4)
                    else:
                        self.specTraceGain[i][k].append(40)



    def calcGSinteraction(self):
        self.overGSstart = []
        self.overGSend = []
        self.overGSdur = []
        self.overGSwhich = []

        for i in range(self.nb_satellites):
            #Read in file
            data_filename = self.gsInteraction_filepath + self.output_filename_propagator[i]
            data_file = open(data_filename, 'r')
            data_lines = data_file.readlines()
            #Pop off Header and blank line
            del data_lines[0]
            del data_lines[0]

            temp_start = []
            temp_end = []
            temp_dur = []
            temp_which = []

            for j in range(len(data_lines)):
                line = data_lines[j].split()

                start_time = line[1] + "T" + line[2]
                temp_start.append(datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S'))

                end_time = line[4] + "T" + line[5]
                temp_end.append(datetime.strptime(end_time, '%Y-%m-%dT%H:%M:%S'))

                temp_dur.append(temp_end[j] - temp_start[j]) #Appends as timedelta object

                temp_which.append(line[12])

            self.overGSstart.append(temp_start)
            self.overGSend.append(temp_end)
            self.overGSdur.append(temp_dur)
            self.overGSwhich.append(temp_which)

    def findSpecStormPasses(self):
        self.storm_names = []
        self.specPassStart = []
        self.specPassGps = []
        self.specPassDist = []
        self.specPassIn = []
        
        for i in range(self.nb_satellites):
            data_filename = self.specInteraction_filepath + self.output_filename_propagator[i]
            data_file = open(data_filename, 'r')
            data_lines = data_file.readlines()
            #Pop off Header line
            del data_lines[0]
            
            if len(self.storm_names) == 0:
                line = data_lines[0].split()
                for m in range(self.nb_storms):
                    name = line[5+2*m].split('_')[1]
                    self.storm_names.append(name)

            del data_lines[0]
            del data_lines[0]

            temp_start = [[] for x in range(self.nb_storms)]
            temp_gps = [[] for x in range(self.nb_storms)]
            temp_dist = [[] for x in range(self.nb_storms)]
            temp_in = [[] for x in range(self.nb_storms)]

            for j in range(len(data_lines) - 1):
                line = data_lines[j].split()
                
                if datetime.strptime(line[0], '%Y-%m-%dT%H:%M:%S') > datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S') and datetime.strptime(line[0], '%Y-%m-%dT%H:%M:%S') < datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S'):            
                    for k in range(self.nb_storms):
                        if line[6+2*k] == "1": #Checks that spec is in storm
                            temp_start[k].append(line[0])
                            temp_gps[k].append(line[4])
                            temp_dist[k].append(line[5+2*k])
                            temp_in[k].append(line[6+2*k])

            self.specPassStart.append(temp_start)
            self.specPassGps.append(temp_gps)
            self.specPassDist.append(temp_dist)
            self.specPassIn.append(temp_in)


    def generateSpecReport(self):
        file_loc = self.specReportLoc + "/" + str(datetime.now().strftime('%Y-%m-%dT%H-%M-%S')) + ".txt"
        outfile = open(file_loc, 'w')
        startTime = datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S')
        endTime = datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S')
        self.findSpecStormPasses()

        #Write header
        outfile.write("CYGNSS Specular Flyover Report\n")
        outfile.write("-----\n")
        outfile.write("Generated on: " + str(datetime.now()) + "\n")
        outfile.write("For time range: " + self.start_visu + " to " + self.stop_visu + "\n")
        outfile.write("Number of Satellites in this report: " + str(self.nb_satellites) + "\n")
        outfile.write("-----\n\n")

        for sat in range(self.nb_satellites):
            outfile.write("###CYGNSS " + str(sat + 1) + "###\n") #Sat+1 converts from 0 to 1 indexing
            for storm in range(len(self.storm_names)):
                outfile.write(self.storm_names[storm])
                outfile.write("\n")
                for spec in range(len(self.specPassStart[sat][storm])):
                    outfile.write(self.specPassStart[sat][storm][spec] + "   " + self.specPassDist[sat][storm][spec] + "   " + self.specPassGps[sat][storm][spec])
                    outfile.write("\n")
                outfile.write("\n")
            outfile.write("---")
            outfile.write("\n")

        outfile.write("-----\n")
        outfile.write("Note: This file was autogenerated by SIFT.\n")

        outfile.close()
        print "Specular/Storm report generated."

        os.startfile(file_loc)

    #Sets appropriate gsinteraction to current indices
    def findGSinteractionIdx(self):
        for i in range(self.nb_satellites):
            self.GSinteractionIdx[i] = 0 #Reset value to base
            tmpDT = datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S')
            while tmpDT > self.overGSstart[i][self.GSinteractionIdx[i]]:
                self.GSinteractionIdx[i] += 1

    def generateGSreport(self):
        file_loc = self.gsReportLoc + "/" + str(datetime.now().strftime('%Y-%m-%dT%H-%M-%S')) + ".txt"
        outfile = open(file_loc, 'w')
        startTime = datetime.strptime(self.start_visu, '%Y-%m-%dT%H:%M:%S')
        endTime = datetime.strptime(self.stop_visu, '%Y-%m-%dT%H:%M:%S')

        #Write header
        outfile.write("CYGNSS Ground Station Flyover Report\n")
        outfile.write("-----\n")
        outfile.write("Generated on: " + str(datetime.now()) + "\n")
        outfile.write("For time range: " + self.start_visu + " to " + self.stop_visu + "\n")
        outfile.write("Number of Satellites in this report: " + str(self.nb_satellites) + "\n")
        outfile.write("-----\n")

        #Find first GSInteractionIdx of interest for all sats
        firstGSIdx = self.GSinteractionIdx

        #Find last GSInteractionIdx of intererst for all sats
        lastGSIdx = []
        for i in range(self.nb_satellites):
            tmpIdx = 0
            while tmpIdx < len(self.overGSend[i]) and endTime > self.overGSend[i][tmpIdx]:
                tmpIdx += 1

            lastGSIdx.append(tmpIdx)

        #Write information to file
        for i in range(self.nb_satellites):
            outfile.write("CYGNSS " + str(i + 1) + "\n")
            for interactionIdx in range(firstGSIdx[i], lastGSIdx[i]):
                outfile.write("\tGround Station: " + self.overGSwhich[i][interactionIdx] + "  Connection Start: " + str(self.overGSstart[i][interactionIdx]) + 
                    "  Connection End: " + str(self.overGSend[i][interactionIdx]) + "  Duration: " + str(self.overGSdur[i][interactionIdx]) + "\n")
                interactionIdx += 1

            outfile.write("\n")

        outfile.write("-----\n")
        outfile.write("Note: This file was autogenerated by SIFT.  It only includes Ground Station interactions that COMPLETELY fall within the analysis interval given at the top of this file\n")

        outfile.close()
        print "Satellite/Ground Station report generated."

        os.startfile(file_loc)

    

##############################################


class CygnssFrame(wx.Frame):
    """ The main application frame
    """
    title = 'CYGNSS Forecaster'
    zoom = 0
    show_graph = 0
    mapRef = None
    localGSInteractionIdxs = None
    meridians = {}
    parallels = {}

    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)

        
        self.initialize_plot_map()
        self.satManager = SatelliteManager(self.plotMap, self.max_lon, self.min_lon, self.min_lat, self.height_fig_with_map, self.width_fig_with_map, self.ax_map, self.zoom)
        self.localGSInteractionIdxs = [-1] * self.satManager.nb_satellites
        self.create_main_panel()

        self.paused = True
        self.revPaused = True
        self.create_menu()
        

        self.redraw_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_redraw_timer, self.redraw_timer)
        self.redraw_timer.Start(10) #Controls how often the UI is redrawn


    def create_menu(self):
        self.menubar = wx.MenuBar()

        menu_file = wx.Menu()
        m_exit = menu_file.Append(-1, "&Exit\t Ctrl-X", "Exit")
        self.Bind(wx.EVT_CLOSE, self.on_exit, m_exit)

        menu_view = wx.Menu()
        m_fullGlobe = menu_view.Append(-1, "&Full Globe", "Full Globe")
        self.Bind(wx.EVT_MENU, self.on_full_globe, m_fullGlobe)
        nwAtlantic = menu_view.Append(-1, "&Northwest Atlantic", "NorthWest Atlantic")
        self.Bind(wx.EVT_MENU, self.on_northwest_atl, nwAtlantic)
        swAtlantic = menu_view.Append(-1, "&Southwest Atlantic", "SouthWest Atlantic")
        self.Bind(wx.EVT_MENU, self.on_southwest_atl, swAtlantic)
        nePacific = menu_view.Append(-1, "&Northeast Pacific", "NorthEast Pacific")
        self.Bind(wx.EVT_MENU, self.on_northeast_pac, nePacific)
        sePacific = menu_view.Append(-1, "&Southeast Pacific", "SouthEast Pacific")
        self.Bind(wx.EVT_MENU, self.on_southeast_pac, sePacific)
        asiaPacific = menu_view.Append(-1, "&Asia Pacific", "Asia Pacific")
        self.Bind(wx.EVT_MENU, self.on_asia_pac, asiaPacific)
        allAtlantic = menu_view.Append(-1, "&All Atlantic", "All Atlantic")
        self.Bind(wx.EVT_MENU, self.on_all_atl, allAtlantic)
        gulf = menu_view.Append(-1, "&Gulf", "Gulf")
        self.Bind(wx.EVT_MENU, self.on_gulf, gulf)
        indian = menu_view.Append(-1, "&Indian Ocean", "Indian Ocean")
        self.Bind(wx.EVT_MENU, self.on_indian, indian)

        menu_generate = wx.Menu()
        stormReport = menu_generate.Append(-1, "&Specular/Storm Report", "Specular/Storm Report")
        self.Bind(wx.EVT_MENU, self.on_generate_spec_report, stormReport)
        gsReport = menu_generate.Append(-1, "&Satellite/Ground Station Report", "Satellite/Ground Station Report")
        self.Bind(wx.EVT_MENU, self.on_generate_gs_report, gsReport)


        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_view, "&View")
        self.menubar.Append(menu_generate, "&Generate")
        self.SetMenuBar(self.menubar)


    def initialize_plot_map(self):
        ### MAP ###

        self.height_fig_with_map = 14
        self.width_fig_with_map = 8
        background_color = (0/555,76./255,153/255.)
        self.fig_with_map = plt.figure(figsize=(self.height_fig_with_map, self.width_fig_with_map))
        self.fig_with_map.set_facecolor(background_color)
        # Axes
        self.gs_with_map = gridspec.GridSpec(1, self.show_graph+1)
        self.gs_with_map.update(left=0.04, right=0.99, top = 0.98,bottom = 0.04,hspace=0.08,wspace = 0.06)
        self.ax_map = self.fig_with_map.add_subplot(self.gs_with_map[0, self.show_graph])


        #This is where we set map params for different zoom settings
        #Order is: Globe, NWAtlantic, SWAtlantic, NEPacific, SEPacific, AllAtlantic, Asia Pacific, Gulf, Indian
        self.min_lon = [-180, -90, -80, -180, -180, -90, 60, -110, 10] #0
        self.max_lon = [180, -10, -0, -100, -100, 30, 180, -30, 140] #60
        self.min_lat = [-90, 0, -50, 0, -50, -40, -40, 0, -70] #-10
        self.max_lat = [90,50, 0, 50, 0, 40, 40, 50, 40]#40
        self.step_lon = [45,2, 2, 2, 2, 2, 2, 2, 5]#2
        self.step_lat = [30, 2, 2, 2, 2, 2, 2, 2, 5]#2

        self.build_plot_map(self.mapRef)

    def build_plot_map(self, mapRef):
        #Del old keys to remove old map lines; must use del here not clear
        for key in self.meridians.keys():
            del self.meridians[key]
        for key in self.parallels.keys():
            del self.parallels[key]

        self.array_lon = [ str(ss) for ss in np.arange(self.min_lon[self.zoom], self.max_lon[self.zoom],self.step_lon[self.zoom]) ]
        self.array_lat = [ str(ss) for ss in np.arange(self.min_lat[self.zoom], self.max_lat[self.zoom],self.step_lat[self.zoom]) ]
        self.ax_map.xaxis.set_major_locator(FixedLocator(np.arange(self.min_lon[self.zoom], self.max_lon[self.zoom]+1, self.step_lon[self.zoom])))
        self.ax_map.yaxis.set_major_locator(FixedLocator(np.arange(self.min_lat[self.zoom], self.max_lat[self.zoom]+1, self.step_lat[self.zoom])))
        self.ax_map.set_xticklabels(self.array_lon, color = 'w')
        self.ax_map.set_yticklabels(self.array_lat, color = 'w')

        self.plotMap = Basemap( projection       = 'cyl',
                 llcrnrlon        = self.min_lon[self.zoom] , #Lower Left  CoRNeR Longitude
                 urcrnrlon        = self.max_lon[self.zoom]  , #Upper Right CoRNeR Longitude
                 llcrnrlat        = self.min_lat[self.zoom]  , #Lower Left  CoRNeR Latitude
                 urcrnrlat        = self.max_lat[self.zoom],   #Upper Right CoRNeR Latitude
                 resolution       = 'l'  ,
                 suppress_ticks   = False,
                 ax = self.ax_map,
                 )
        self.meridians = self.plotMap.drawmeridians(np.arange(self.min_lon[self.zoom], self.max_lon[self.zoom],self.step_lon[self.zoom]))
        self.parallels = self.plotMap.drawparallels(np.arange(self.min_lat[self.zoom], self.max_lat[self.zoom],self.step_lat[self.zoom]))
        mapRef = self.plotMap.drawcoastlines(linewidth=0.7, color='blue')

    def create_main_panel(self):
        self.panel = wx.Panel(self)
        self.canvas = FigCanvas(self.panel, wx.ID_ANY, self.fig_with_map)

        #Create infrastructure for staticbox outline effect (this must be done first)
        self.paramSb = wx.StaticBox(self.panel, -1, 'Parameters: ')
        self.graphicSb = wx.StaticBox(self.panel, -1, 'Visualization: ')
        self.readoutSb = wx.StaticBox(self.panel, -1, 'Output: ')

        #Create buttons
        self.ppBtn = wx.Button(self.panel, -1, "Play")
        self.Bind(wx.EVT_BUTTON, self.on_pause_button, self.ppBtn)
        self.revBtn = wx.Button(self.panel, -1, "Play Reverse")
        self.Bind(wx.EVT_BUTTON, self.on_rev_button, self.revBtn)

        self.resetBtn = wx.Button(self.panel, -1, "Reset")
        self.Bind(wx.EVT_BUTTON, self.on_reset_button, self.resetBtn)

        self.resetIntervalBtn = wx.Button(self.panel, -1, "Reset Interval")
        self.Bind(wx.EVT_BUTTON, self.on_reset_interval_button, self.resetIntervalBtn)

        self.startValBtn = wx.Button(self.panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.on_startVal_set, self.startValBtn)

        self.endValBtn = wx.Button(self.panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.on_endVal_set, self.endValBtn)

        self.curValBtn = wx.Button(self.panel, -1, "Set")
        self.Bind(wx.EVT_BUTTON, self.on_curVal_set, self.curValBtn)

        self.curPThreeBtn = wx.Button(self.panel, -1, "Set: Current + 3h", name="curPThree")
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.curPThreeBtn)
        self.curPTwelveBtn = wx.Button(self.panel, -1, "Set: Current + 12h", name="curPTwelve")
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.curPTwelveBtn)
        self.nowPThreeBtn = wx.Button(self.panel, -1, "Set: Now + 3h", name="nowPThree")
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.nowPThreeBtn)
        self.nowPTwelveBtn = wx.Button(self.panel, -1, "Set: Now + 12h", name="nowPTwelve")
        self.Bind(wx.EVT_BUTTON, self.on_hotBtn, self.nowPTwelveBtn)

       

        #Create static text labels
        self.satLabel = wx.StaticText(self.panel, label="Show/Hide Ground Tracks")
        self.gsLabel = wx.StaticText(self.panel, label="Show/Hide Ground Stations")
        self.speedLabel = wx.StaticText(self.panel, label="Animation Speed")
        self.startLabel = wx.StaticText(self.panel, label="Visualization Start:")
        self.endLabel = wx.StaticText(self.panel, label="Visualization End:  ")
        self.jumpLabel = wx.StaticText(self.panel, label="Jump to:                ") #Extra spaces to align textboxes
        self.gsinteractionLabel = wx.StaticText(self.panel, label="Ground Station Interactions ")
        self.stormColorKey = wx.StaticText(self.panel, label="Storm Color Key")

        #For output section
        self.satCol = wx.StaticText(self.panel, label="Satellite")
        self.nextGSCol = wx.StaticText(self.panel, label="Next Ground Station")
        self.nextPassCol = wx.StaticText(self.panel, label="Time until next pass")
        self.durationCol = wx.StaticText(self.panel, label="Duration")

        self.cyg1lab = wx.StaticText(self.panel, label="CYGNSS 1")
        self.cyg2lab = wx.StaticText(self.panel, label="CYGNSS 2")
        self.cyg3lab = wx.StaticText(self.panel, label="CYGNSS 3")
        self.cyg4lab = wx.StaticText(self.panel, label="CYGNSS 4")
        self.cyg5lab = wx.StaticText(self.panel, label="CYGNSS 5")
        self.cyg6lab = wx.StaticText(self.panel, label="CYGNSS 6")
        self.cyg7lab = wx.StaticText(self.panel, label="CYGNSS 7")
        self.cyg8lab = wx.StaticText(self.panel, label="CYGNSS 8")
        self.cyg1lab.SetForegroundColour(wx.Colour(0, 255, 0)) #lime
        self.cyg2lab.SetForegroundColour(wx.Colour(0, 0, 255)) #blue
        self.cyg3lab.SetForegroundColour(wx.Colour(75, 0, 130)) #indigo
        self.cyg4lab.SetForegroundColour(wx.Colour(128, 0, 128)) #purple
        self.cyg5lab.SetForegroundColour(wx.Colour(30, 144, 255)) #dodgerblue
        self.cyg6lab.SetForegroundColour(wx.Colour(70, 130, 180)) #steelblue
        self.cyg7lab.SetForegroundColour(wx.Colour(60, 179, 113)) #seagreen
        self.cyg8lab.SetForegroundColour(wx.Colour(50, 205, 50)) #limegreen

        
        self.cyg1NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg2NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg3NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg4NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg5NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg6NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg7NextTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg8NextTC = wx.StaticText(self.panel, size=(100, 20))

        #For Storm Color Key
        self.HUlab = wx.StaticText(self.panel, label="Hurricane")
        self.TSlab = wx.StaticText(self.panel, label="Tropical Storm")
        self.TDlab = wx.StaticText(self.panel, label="Tropical Depression")
        self.HUlab.SetForegroundColour(wx.Colour(255, 0, 0)) #red
        self.TSlab.SetForegroundColour(wx.Colour(255, 165, 0)) #orange
        self.TDlab.SetForegroundColour(wx.Colour(0, 255, 255)) #aqua
      

        self.cygNextTCs = [self.cyg1NextTC, self.cyg2NextTC, self.cyg3NextTC, self.cyg4NextTC, self.cyg5NextTC, self.cyg6NextTC, self.cyg7NextTC, self.cyg8NextTC]

        self.cyg1TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg2TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg3TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg4TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg5TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg6TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg7TimeTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg8TimeTC = wx.StaticText(self.panel, size=(100, 20))

        self.cygTimeTCs = [self.cyg1TimeTC, self.cyg2TimeTC, self.cyg3TimeTC, self.cyg4TimeTC, self.cyg5TimeTC, self.cyg6TimeTC, self.cyg7TimeTC, self.cyg8TimeTC]

        self.cyg1DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg2DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg3DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg4DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg5DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg6DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg7DurTC = wx.StaticText(self.panel, size=(100, 20))
        self.cyg8DurTC = wx.StaticText(self.panel, size=(100, 20))

        self.cygDurTCs = [self.cyg1DurTC, self.cyg2DurTC, self.cyg3DurTC, self.cyg4DurTC, self.cyg5DurTC, self.cyg6DurTC, self.cyg7DurTC, self.cyg8DurTC]

        #Call to set up readout correctly
        self.updateGSNextAndDur()
        self.updateGSTime()

        self.vizParamLab = wx.StaticText(self.panel, size=(150,20), label="Visualization Parameters")
        self.startParamLab = wx.StaticText(self.panel, size=(120,20), label="Current Analysis Start: ")
        self.endParamLab = wx.StaticText(self.panel, size=(120,20), label="Current Analysis End: ")
        self.maxStartParamLab = wx.StaticText(self.panel, size=(100,20), label="Max Analysis Start: ")
        self.maxEndParamLab = wx.StaticText(self.panel, size=(100,20), label = "Max Analysis End: ")

        self.startParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.start_visu)
        self.endParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.stop_visu)
        self.maxStartParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.firstTime)
        self.maxEndParam = wx.StaticText(self.panel, size=(110,20), label=self.satManager.lastTime)

        #Create checkboxes
        self.cyg1cb = wx.CheckBox(self.panel, label='CYGNSS 1')
        self.cyg1cb.SetForegroundColour(wx.Colour(0, 255, 0)) #lime
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg1cb)

        self.cyg2cb = wx.CheckBox(self.panel, label='CYGNSS 2')
        self.cyg2cb.SetForegroundColour(wx.Colour(0, 0, 255)) #blue
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg2cb)

        self.cyg3cb = wx.CheckBox(self.panel, label='CYGNSS 3')
        self.cyg3cb.SetForegroundColour(wx.Colour(75, 0, 130)) #indigo
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg3cb)

        self.cyg4cb = wx.CheckBox(self.panel, label='CYGNSS 4')
        self.cyg4cb.SetForegroundColour(wx.Colour(128, 0, 128)) #purple
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg4cb)

        self.cyg5cb = wx.CheckBox(self.panel, label='CYGNSS 5')
        self.cyg5cb.SetForegroundColour(wx.Colour(30, 144, 255)) #dodgerblue
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg5cb)

        self.cyg6cb = wx.CheckBox(self.panel, label='CYGNSS 6')
        self.cyg6cb.SetForegroundColour(wx.Colour(70, 130, 180)) #steelblue
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg6cb)

        self.cyg7cb = wx.CheckBox(self.panel, label='CYGNSS 7')
        self.cyg7cb.SetForegroundColour(wx.Colour(60, 179, 113)) #seagreen
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg7cb)

        self.cyg8cb = wx.CheckBox(self.panel, label='CYGNSS 8')
        self.cyg8cb.SetForegroundColour(wx.Colour(50, 205, 50)) #limegreen
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg, self.cyg8cb)

        self.cyg1specCB = wx.CheckBox(self.panel, label='Specs 1')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg1specCB)
        self.cyg2specCB = wx.CheckBox(self.panel, label='Specs 2')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg2specCB)
        self.cyg3specCB = wx.CheckBox(self.panel, label='Specs 3')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg3specCB)
        self.cyg4specCB = wx.CheckBox(self.panel, label='Specs 4')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg4specCB)
        self.cyg5specCB = wx.CheckBox(self.panel, label='Specs 5')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg5specCB)
        self.cyg6specCB = wx.CheckBox(self.panel, label='Specs 6')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg6specCB)
        self.cyg7specCB = wx.CheckBox(self.panel, label='Specs 7')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg7specCB)
        self.cyg8specCB = wx.CheckBox(self.panel, label='Specs 8')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_cyg_spec, self.cyg8specCB)

        self.gs1cb = wx.CheckBox(self.panel, label='Hawaii')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gsHawaii, self.gs1cb)
        self.gs2cb = wx.CheckBox(self.panel, label='Chile')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gsSantiago, self.gs2cb)
        self.gs3cb = wx.CheckBox(self.panel, label='Australia')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gsPerth, self.gs3cb)

        self.stormcb = wx.CheckBox(self.panel, label='Show/Hide Storms')
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_storms, self.stormcb)

        self.cygCheckBoxes = [self.cyg1cb, self.cyg2cb, self.cyg3cb, self.cyg4cb, self.cyg5cb, self.cyg6cb, self.cyg7cb, self.cyg8cb]
        self.gsCheckBoxes = [self.gs1cb, self.gs2cb, self.gs3cb]
        self.specCheckBoxes = [self.cyg1specCB, self.cyg2specCB, self.cyg3specCB, self.cyg4specCB, self.cyg5specCB, self.cyg6specCB, self.cyg7specCB, self.cyg8specCB]

        #Create sliders
        self.speedSld = wx.Slider(self.panel, value=self.satManager.speedFactor, minValue=self.satManager.speedFactor, maxValue=100, style=wx.SL_HORIZONTAL)
        self.speedSld.Bind(wx.EVT_SCROLL, self.on_speed_scroll)

        #Create text entry fields
        self.startTC = wx.TextCtrl(self.panel, value="YYYY-MM-DDTHH:MM:SS", size=(130,20), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_startVal_set, self.startTC)
        self.endTC = wx.TextCtrl(self.panel, value="YYYY-MM-DDTHH:MM:SS", size=(130,20), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_endVal_set, self.endTC)
        self.curTC = wx.TextCtrl(self.panel, value="YYYY-MM-DDTHH:MM:SS", size=(130,20), style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT_ENTER, self.on_curVal_set, self.curTC)


        #Main Sizers
        self.windowSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.controlsSizer = wx.StaticBoxSizer(self.paramSb, wx.VERTICAL)
        self.outputSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.graphicSizer = wx.StaticBoxSizer(self.graphicSb, wx.VERTICAL)
        self.textoutSizer = wx.StaticBoxSizer(self.readoutSb, wx.VERTICAL)

        #Row sizers for controls section
        self.cRow1 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow2 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow3 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow4 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow5 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow6 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow7 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow8 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow9 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow10 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow11 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow12 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow13 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow14 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow15 = wx.BoxSizer(wx.HORIZONTAL)
        self.cRow16 = wx.BoxSizer(wx.HORIZONTAL)

        #Row sizers for output section
        self.oRow0 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow1 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow2 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow3 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow4 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow5 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow6 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow7 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow8 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow9 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow10 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow11 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow12 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow13 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow14 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow15 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow16 = wx.BoxSizer(wx.HORIZONTAL)
        self.oRow17 = wx.BoxSizer(wx.HORIZONTAL)

        #Fill sizers in order from smallest to largest
        #Control rows
        self.cRow1.Add(self.cyg1cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow1.Add(self.cyg1specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow1.Add(self.cyg5cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow1.Add(self.cyg5specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow2.Add(self.cyg2cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow2.Add(self.cyg2specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow2.Add(self.cyg6cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow2.Add(self.cyg6specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow3.Add(self.cyg3cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow3.Add(self.cyg3specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow3.Add(self.cyg7cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow3.Add(self.cyg7specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow4.Add(self.cyg4cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow4.Add(self.cyg4specCB, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow4.Add(self.cyg8cb, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow4.Add(self.cyg8specCB, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow5.Add(self.gs1cb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow6.Add(self.gs2cb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow7.Add(self.gs3cb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow8.Add(self.stormcb, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow9.Add(self.speedLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow9.Add(self.speedSld, 0, wx.ALL | wx.EXPAND, border = 8)

        self.cRow10.Add(self.startLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow10.Add(self.startTC, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow10.Add(self.startValBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow11.Add(self.endLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow11.Add(self.endTC, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow11.Add(self.endValBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow12.Add(self.jumpLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow12.Add(self.curTC, 0, wx.ALL | wx.EXPAND, border = 8)
        self.cRow12.Add(self.curValBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow13.Add(self.ppBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow13.Add(self.revBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow16.Add(self.resetBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow16.Add(self.resetIntervalBtn, 0, wx.ALL | wx.EXPAND, border = 4)

        self.cRow14.Add(self.curPThreeBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow14.Add(self.curPTwelveBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow15.Add(self.nowPThreeBtn, 0, wx.ALL | wx.EXPAND, border = 4)
        self.cRow15.Add(self.nowPTwelveBtn, 0, wx.ALL | wx.EXPAND, border = 4)


        #Output rows
        self.oRow0.Add(self.satCol, 0, wx.ALL | wx.EXPAND, border= 4)
        self.oRow0.Add(self.nextGSCol, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow0.Add(self.nextPassCol, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow0.Add(self.durationCol, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow1.Add(self.cyg1lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow1.Add(self.cyg1NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow1.Add(self.cyg1TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow1.Add(self.cyg1DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow2.Add(self.cyg2lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow2.Add(self.cyg2NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow2.Add(self.cyg2TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow2.Add(self.cyg2DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow3.Add(self.cyg3lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow3.Add(self.cyg3NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow3.Add(self.cyg3TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow3.Add(self.cyg3DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow4.Add(self.cyg4lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow4.Add(self.cyg4NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow4.Add(self.cyg4TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow4.Add(self.cyg4DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow5.Add(self.cyg5lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow5.Add(self.cyg5NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow5.Add(self.cyg5TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow5.Add(self.cyg5DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow6.Add(self.cyg6lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow6.Add(self.cyg6NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow6.Add(self.cyg6TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow6.Add(self.cyg6DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow7.Add(self.cyg7lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow7.Add(self.cyg7NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow7.Add(self.cyg7TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow7.Add(self.cyg7DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow8.Add(self.cyg8lab, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow8.Add(self.cyg8NextTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow8.Add(self.cyg8TimeTC, 0, wx.ALL | wx.EXPAND, border=4)
        self.oRow8.Add(self.cyg8DurTC, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow9.Add(self.vizParamLab, 0, wx.ALL | wx.EXPAND, border=4)

        self.oRow10.Add(self.startParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow10.Add(self.startParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow11.Add(self.endParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow11.Add(self.endParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow12.Add(self.maxStartParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow12.Add(self.maxStartParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow13.Add(self.maxEndParamLab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow13.Add(self.maxEndParam, 0, wx.ALL | wx.EXPAND, border=2)

        self.oRow14.Add(self.stormColorKey, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow15.Add(self.HUlab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow16.Add(self.TSlab, 0, wx.ALL | wx.EXPAND, border=2)
        self.oRow17.Add(self.TDlab, 0, wx.ALL | wx.EXPAND, border=2)


        self.controlsSizer.Add(self.satLabel, 0, wx.ALL, border = 10)
        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow1, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow2, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow3, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow4, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.gsLabel, 0, wx.ALL, border = 10)
        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow5, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow6, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow7, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow8, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.controlsSizer.Add(self.cRow9, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow10, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow11, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow12, 0, wx.ALL | wx.EXPAND, border = 2)
        self.controlsSizer.Add(self.cRow13, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow16, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow14, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)
        self.controlsSizer.Add(self.cRow15, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)


        self.textoutSizer.Add(self.gsinteractionLabel, 0, wx.ALL, border=5)
        self.textoutSizer.Add(self.oRow0, 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow1, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow2, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow3, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow4, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow5, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow6, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow7, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow8, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow9, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow10, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow11, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow12, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow13, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        self.textoutSizer.Add(self.oRow14, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow15, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow16, 0, wx.ALL | wx.EXPAND, border=2)
        self.textoutSizer.Add(self.oRow17, 0, wx.ALL | wx.EXPAND, border=2)

        self.graphicSizer.Add(self.canvas, 1, wx.ALL | wx.EXPAND, border = 15)
        #self.graphicSizer.Add(wx.StaticLine(self.panel,), 0, wx.ALL | wx.EXPAND, border=5)
        #self.graphicSizer.Add(self.outputLabel, 0, wx.ALL | wx.EXPAND, border = 8)
        #self.graphicSizer.Add(self.gRow0, 0, wx.TOP | wx.BOTTOM | wx.LEFT | wx.EXPAND, border = 8)

        self.outputSizer.Add(self.graphicSizer, 0, wx.ALL | wx.EXPAND, border = 5)
        self.outputSizer.Add(self.textoutSizer, 0, wx.ALL | wx.EXPAND, border = 5)

        self.windowSizer.Add(self.controlsSizer, 0, wx.ALL | wx.EXPAND, border = 5)
        self.windowSizer.Add(self.outputSizer, 0, wx.ALL | wx.EXPAND, border = 5)

        self.panel.SetSizer(self.windowSizer)
        self.windowSizer.Fit(self)

    #Helper functions ###################################

    def draw_plot(self):
        self.satManager.advSim(self.plotMap)
        self.updateGSNextAndDur()
        self.updateGSTime()

        self.canvas.draw()

    def on_redraw_timer(self, event):
        if not self.paused:
            self.satManager.incrementIndex()
        elif not self.revPaused:
            self.satManager.decrementIndex()

        self.draw_plot()


    #GUI Control Functionality ###########################    

    def on_pause_button(self, event):
        self.paused = not self.paused

        label = "Play" if self.paused else "Pause"
        self.ppBtn.SetLabel(label)

        if self.paused:
            self.revBtn.Enable()
        else:
            self.revBtn.Disable()

    def on_rev_button(self, event):
        self.revPaused = not self.revPaused

        label = "Play Reversed" if self.revPaused else "Pause"
        self.revBtn.SetLabel(label)

        if self.revPaused:
            self.ppBtn.Enable()
        else:
            self.ppBtn.Disable()

    def on_reset_button(self, event):
        #We choose to automatically pause when resetting
        self.paused = True
        self.ppBtn.SetLabel("Play")

        self.satManager.curSample = self.satManager.when_start_visu

        #Perform updates to GS indexing to ensure readout is correct
        self.satManager.findGSinteractionIdx()
        self.updateGSNextAndDur()
        self.updateGSTime()

    def on_reset_interval_button(self, event):
        self.satManager.start_visu = self.satManager.firstTime
        self.satManager.stop_visu = self.satManager.lastTime
        #These values are indices
        self.satManager.when_start_visu = self.satManager.time_sat.index(self.satManager.start_visu)
        self.satManager.when_stop_visu = self.satManager.time_sat.index(self.satManager.stop_visu)

        #Update visualization parameter readout
        self.startParam.SetLabel(self.satManager.firstTime)
        self.endParam.SetLabel(self.satManager.lastTime)

        print "New Analysis Interval: " + self.satManager.time_sat[self.satManager.when_start_visu] + ' -> ' + self.satManager.time_sat[self.satManager.when_stop_visu]

    def on_hotBtn(self, event):
        caller = event.GetEventObject().GetName()
        
        if caller == "curPThree":
            newStartTime = datetime.strptime(self.satManager.time_sat[self.satManager.curSample],'%Y-%m-%dT%H:%M:%S')
            newEndTime = newStartTime + timedelta(hours=3)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))

        elif caller == "curPTwelve":
            newStartTime = datetime.strptime(self.satManager.time_sat[self.satManager.curSample],'%Y-%m-%dT%H:%M:%S')
            newEndTime = newStartTime + timedelta(hours=12)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))

        elif caller == "nowPThree":
            newStartTime = datetime.utcnow()
            newEndTime = newStartTime + timedelta(hours=3)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))


        elif caller == "nowPTwelve":
            newStartTime = datetime.utcnow()
            newEndTime = newStartTime + timedelta(hours=12)
            self.time_set_helper(newStartTime.strftime('%Y-%m-%dT%H:%M:%S'), newEndTime.strftime('%Y-%m-%dT%H:%M:%S'))

    def time_set_helper(self, newStartTime, newEndTime):
        try:
            indexStartTime = self.satManager.time_sat.index(newStartTime)
            indexEndTime = self.satManager.time_sat.index(newEndTime)
            #If there is no exception in the above statement, set values
            #Note: start_visu is format "YYYY-MM-DDTHH:MM:SS"
            #      when_start_visu is index corresponding to that timestamp
            self.satManager.start_visu = newStartTime
            self.satManager.when_start_visu = indexStartTime
            self.satManager.curSample = indexStartTime

            self.satManager.stop_visu = newEndTime
            self.satManager.when_stop_visu = indexEndTime

            #Recalculate ground traces
            self.satManager.calcTraces(self.satManager.when_start_visu, self.satManager.when_stop_visu)
            self.reset_traces()

            #Recalculate GS interaction indices
            self.satManager.findGSinteractionIdx()

            #Update visualization param readout
            self.startParam.SetLabel(newStartTime)
            self.endParam.SetLabel(newEndTime)
            
        except ValueError:
            print "OPERATION FAILED: " + newStartTime + " -> " + newEndTime + " is outside allowable range; Analysis start time has not changed."
            return

        print "New Analysis Interval: " + self.satManager.time_sat[self.satManager.when_start_visu] + ' -> ' + self.satManager.time_sat[self.satManager.when_stop_visu]


    def on_cb_cyg(self, event):
        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        senderLabel = sender.GetLabel()

        cbId = [int(s) for s in senderLabel.split() if s.isdigit()]


        if(isChecked):
            self.satManager.plotSatTrace(self.plotMap, cbId[0] - 1)
        else:
            self.satManager.removeSatTrace(cbId[0] - 1)

    def on_cb_cyg_spec(self, event):
        sender = event.GetEventObject()
        isChecked=sender.GetValue()

        senderLabel = sender.GetLabel()

        cbId = [int(s) for s in senderLabel.split() if s.isdigit()]

        if(isChecked):
            self.satManager.plotSpecTrace(self.plotMap, cbId[0] - 1)
        else:
            self.satManager.removeSpecTrace(cbId[0] - 1)

    def on_cb_gsHawaii(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotGroundStation(self.plotMap, 0)
        else:
            self.satManager.removeGroundStation(0)

        self.draw_plot()

    def on_cb_gsSantiago(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotGroundStation(self.plotMap, 1)
        else:
            self.satManager.removeGroundStation(1)

        self.draw_plot()

    def on_cb_gsPerth(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotGroundStation(self.plotMap, 2)
        else:
            self.satManager.removeGroundStation(2)

    def on_cb_storms(self, event):

        sender = event.GetEventObject()
        isChecked = sender.GetValue()

        if(isChecked):
            self.satManager.plotTrajectory(self.plotMap)
            self.satManager.plotStorms(self.plotMap)
        else:
            self.satManager.removeTrajectory()
            self.satManager.removeStorms()

    def on_speed_scroll(self, event):

        sender = event.GetEventObject()
        newFactor = sender.GetValue()

        self.satManager.speedFactor = newFactor

    def on_startVal_set(self, event):
        newTime = self.startTC.GetValue()

        #Check if input time is a valid timestamp
        try:
            indexTime = self.satManager.time_sat.index(newTime)
            #If there is no exception in the above statement, set values
            #Note: start_visu is format "YYYY-MM-DDTHH:MM:SS"
            #      when_start_visu is index corresponding to that timestamp
            self.satManager.start_visu = newTime
            self.satManager.when_start_visu = indexTime
            self.satManager.curSample = indexTime

            #Recalculate ground traces
            self.satManager.calcTraces(self.satManager.when_start_visu, self.satManager.when_stop_visu)
            self.reset_traces()

            #Recalculate GS interaction indices
            self.satManager.findGSinteractionIdx()

            #Recalculate spec/storm interaction
            #self.satManager.findSpecStormPasses()

            #Update visualization param readout
            self.startParam.SetLabel(newTime)

            print "New Analysis Interval: " + self.satManager.time_sat[self.satManager.when_start_visu] + ' -> ' + self.satManager.time_sat[self.satManager.when_stop_visu]
    
        except ValueError:
            print "OPERATION FAILED: Start time entered is outside allowable range; Analysis start time has not changed."

        
            
    def on_endVal_set(self, event):
        newTime = self.endTC.GetValue()

        #Check if input time is a valid timestamp
        try:
            indexTime = self.satManager.time_sat.index(newTime)
            #If there is no exception in the above statement, set values
            #Note: stop_visu is format "YYYY-MM-DDTHH:MM:SS"
            #      when_stop_visu is index corresponding to that timestamp
            self.satManager.stop_visu = newTime
            self.satManager.when_stop_visu = indexTime

            #Recalculate ground traces
            self.satManager.calcTraces(self.satManager.when_start_visu, self.satManager.when_stop_visu)
            self.reset_traces()

            #Recalculate spec/storm interaction
            #self.satManager.findSpecStormPasses()
            #print self.satManager.specPassDist

            self.endParam.SetLabel(newTime)

            print "New Analysis Interval: " + self.satManager.time_sat[self.satManager.when_start_visu] + ' -> ' + self.satManager.time_sat[self.satManager.when_stop_visu] 
    
        except ValueError:
            print "OPERATION FAILED: End time entered is outside allowable range; Analysis end time has not changed."

            #Update visualization param readout
            self.endParam.SetLabel(newTime)

        

    def on_curVal_set(self, event):
        newTime = self.curTC.GetValue()

        #First check to make sure input is a valid value at all
        try:
            curTimeIdx = self.satManager.time_sat.index(newTime)
            
        except ValueError:
            print "OPERATION FAILED: Invalid timestamp entered; Please enter a timestamp between " + self.satManager.start_visu + " and " + self.satManager.stop_visu + "."
            return

        #Now test if input is between CURRENT start and stop times
        if curTimeIdx <= self.satManager.when_stop_visu and curTimeIdx >= self.satManager.when_start_visu:
            self.satManager.curSample = curTimeIdx
        else:
            print "OPERATION FAILED: Invalid timestamp entered; Please enter a timestamp between " + self.satManager.start_visu + " and " + self.satManager.stop_visu + "."


    def on_exit(self, event):
        self.Destroy()
        exit(0)

    def on_full_globe(self, event):
        self.zoom = 0
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)
        
    def on_northwest_atl(self, event):
        self.zoom = 1
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_southwest_atl(self, event):
        self.zoom = 2
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_all_atl(self, event):
        self.zoom = 5
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_northeast_pac(self, event):
        self.zoom = 3
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_southeast_pac(self, event):
        self.zoom = 4
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_asia_pac(self, event):
        self.zoom = 6
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_gulf(self, event):
        self.zoom = 7
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)

    def on_indian(self, event):
        self.zoom = 8
        self.build_plot_map(self.mapRef)
        self.satManager.rebuildTimeDisplay(self.zoom)
        
    def on_generate_gs_report(self, event):
        self.satManager.generateGSreport()

    def on_generate_spec_report(self, event):
        self.satManager.generateSpecReport()

    def update_map(self):
        self.initialize_plot_map()
        #self.satManager = SatelliteManager(self.plotMap, self.max_lon, self.min_lat, self.height_fig_with_map, self.width_fig_with_map, self.ax_map)
        self.canvas = FigCanvas(self.panel, wx.ID_ANY, self.fig_with_map)

    def reset_sat_cb(self):
        for i in range(len(self.cygCheckBoxes)):
            self.cygCheckBoxes[i].SetValue(False)

    def reset_traces(self):
        for idx in range(len(self.cygCheckBoxes)):
            #print self.cygCheckBoxes[idx].IsChecked()
            if self.cygCheckBoxes[idx].IsChecked():
                self.satManager.removeSatTrace(idx)
                self.satManager.plotSatTrace(self.plotMap, idx)

            if self.specCheckBoxes[idx].IsChecked():
                self.satManager.removeSpecTrace(idx)
                self.satManager.plotSpecTrace(self.plotMap, idx)

    #GUI Info Readout Functionality #######################

    def updateGSNextAndDur(self):
        for i in range(self.satManager.nb_satellites):
            if(self.satManager.GSinteractionIdx[i] != self.localGSInteractionIdxs[i]):
                self.localGSInteractionIdxs[i] = self.satManager.GSinteractionIdx[i]
                self.cygNextTCs[i].SetLabel(self.satManager.overGSwhich[i][self.satManager.GSinteractionIdx[i]])
                self.cygDurTCs[i].SetLabel(str(self.satManager.overGSdur[i][self.satManager.GSinteractionIdx[i]]))
        

    def updateGSTime(self):
        for i in range(self.satManager.nb_satellites):
            time = (self.satManager.overGSstart[i][self.satManager.GSinteractionIdx[i]]) - datetime.strptime(self.satManager.time_sat[self.satManager.curSample], '%Y-%m-%dT%H:%M:%S')
            self.cygTimeTCs[i].SetLabel(str(time))


if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = CygnssFrame()
    app.frame.Show()
    app.MainLoop()
