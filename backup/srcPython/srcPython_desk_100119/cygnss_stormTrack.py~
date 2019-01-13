# -*- coding: utf-8 -*-
#import pdb
debug = 0
import gzip
import ipdb
import os
import sys
import subprocess
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from spock_cygnss_spec_parallel_for_sift import *
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import FixedLocator
from mpl_toolkits.basemap import Basemap, shiftgrid
from spock_main_input import *
import numpy as np
from collections import *
import time
import datetime
from datetime import timedelta
from itertools import groupby
from operator import itemgetter
from math import *
from scipy.interpolate import *
from bisect import *
import urllib
import matplotlib.gridspec as gridspec
import logging
import ctypes
import wx.lib.scrolledpanel
from struct import *
#Non-Library function imports
from get_ellipse_coords import *
from radius_for_tissot import *
from formatMapString import *
from read_input_file import *
from find_in_read_input_order_variables import *
from storm_track_downloader import DownloadStorms
from readZoom import readZoom
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from os import listdir
from os.path import isfile, join
from scipy.spatial import distance
import matplotlib.patches as mpatches



save_ani = 1
show_ani = 0
show_spec = 1
if show_ani == 1:
    plt.ion()

class interact(object):
    def __init__(self):
        self.userZoom_filepath = 'user_zoom.txt'
        self.mapRef = None            
        self.meridians = {}
        self.parallels = {}
        self.curSample = 0
        self.speedFactor = (int)(3*3600)
        #Assumes the year you are targeting is the year you are looking at
        self.forecastStart = '2017-01-01T00:00:00'# !!!!datetime.today().isoformat()[:11] + '00:00:00'
        self.forecastEnd = '2020-12-31T23:59:59'#!!!!!!(datetime.today() + timedelta(days = 100)).isoformat()[:11] + '00:00:00'
        self.targetYear = '2017'# self.forecastStart[:4]    #'2018'#self.forecastStart[:4]
        self.inputFilepath = 'data'        
        # ask for password to download storm files 
        # v1.2
        self.username_jtwc = 'bussyvirat.umich.edu'
        self.password_jtwc = '[CygNss]2017[D474]'



        # Download storm forecasts for this simulation
        # v1.2
        #ipdb.set_trace()
        #DownloadStorms(self.inputFilepath, 'nhc_observations', self.targetYear, self.username_jtwc, self.password_jtwc)
        #sys.exit()
        # ipdb.set_trace()
        # DownloadStorms(self.inputFilepath, 'nhc', self.targetYear, self.username_jtwc, self.password_jtwc)
        
        # DownloadStorms(self.inputFilepath, 'jtwc', self.targetYear, self.username_jtwc, self.password_jtwc) #Note: requires login, to suppress login request simply comment out


        self.list_storm = []
        self.list_storm.append(sys.argv[1])
        print self.list_storm

        # self.list_storm_old = [f for f in listdir('data_previous') if isfile(join('data_previous', f))]
        # self.list_storm_new = [f for f in listdir('data') if isfile(join('data', f))]
        # self.list_storm = []
        # for ifilenew in self.list_storm_new:
        #     if ((ifilenew in self.list_storm_old) == False):
        #         self.list_storm.append(ifilenew)
        # #os.system("cp -p data/* data_previous") # next time we call this script data will be compared to the previous data
        # if len(self.list_storm) == 0: # no new storm since last time the script was run
        #     print "***! No new storm. The program will stop. !***"; sys.exit()




        self.nb_storms = 0
        self.build_raw_storms()

        startDate_date = self.all_date_storm_arrange[0]
        endDate_date = self.all_date_storm_arrange[-1]
        self.time_sat = []
        istep = 0
        new_date = startDate_date


        self.startDate = startDate_date.isoformat()
        self.endDate = (endDate_date+timedelta(seconds = 2*60)).isoformat() # need to go 2 time steps beyond in SpOCK to make sure that last time step of the storm is included in the output of SpOCK
        
        # Run SpOCK from start of first storm to end of last storm
        list_dont_run = ['spock_spec_start_2017-04-16T06_00_00_end_2017-04-22T18_02_00',
                         'spock_spec_start_2017-06-18T18_00_00_end_2017-06-20T09_02_00',
                         'spock_spec_start_2018-05-20T18_00_00_end_2018-05-29T06_02_00',
                         'spock_spec_start_2018-07-02T00_00_00_end_2018-07-16T06_02_00',
                         'spock_spec_start_2017-06-19T18_00_00_end_2017-06-24T06_02_00',
                         'spock_spec_start_2018-07-04T12_00_00_end_2018-07-12T12_02_00',
                         'spock_spec_start_2017-07-05T12_00_00_end_2017-07-07T12_02_00',
                         'spock_spec_start_2018-08-02T18_00_00_end_2018-08-09T18_02_00',
                         'spock_spec_start_2017-07-17T00_00_00_end_2017-07-18T12_02_00',
                         'spock_spec_start_2018-08-11T12_00_00_end_2018-08-18T06_02_00',
                         'spock_spec_start_2017-07-30T18_00_00_end_2017-08-02T00_02_00',
                         'spock_spec_start_2018-08-29T06_00_00_end_2018-09-18T12_02_00',
                         'spock_spec_start_2017-08-06T18_00_00_end_2017-08-10T12_02_00',
                         'spock_spec_start_2018-09-01T12_00_00_end_2018-09-08T00_02_00',
                         'spock_spec_start_2017-08-12T00_00_00_end_2017-08-18T18_02_00',
                         'spock_spec_start_2018-09-06T00_00_00_end_2018-09-16T12_02_00',
                         'spock_spec_start_2017-08-16T06_00_00_end_2017-09-02T12_02_00'
        ]
        if show_spec == 1:
            self.driverFile = "spock_spec_start_" + self.startDate.replace(":","_") + "_end_" + self.endDate.replace(":","_")
            if (self.driverFile in list_dont_run) == False:
                spock_cygnss_spec_parallel_for_sift(self.startDate, self.endDate, 'spec')

            self.driver_filepath = self.driverFile +'.txt'
            ###Accumulate satellite data###
            self.satNames = ['CYGFM05', 'CYGFM04', 'CYGFM02', 'CYGFM01', 'CYGFM08', 'CYGFM06', 'CYGFM07', 'CYGFM03']
            self.satNameInt = []
            for iii in range(len(self.satNames)):
                self.satNameInt.append((int)((self.satNames[iii]).replace("CYGFM", ""))-1)
            self.satNameInt_back = [3, 2, 7, 1, 0, 5, 6, 4]
            self.specNames = ['SPEC1', 'SPEC2', 'SPEC3', 'SPEC4']

            self.input_variables, self.order_input_variables = read_input_file(self.driver_filepath)
            self.date_start = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'date_start')]
            self.date_stop = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'date_stop')]
            self.dt = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'dt_output')] # used to be dt but cbv changed on june 3 2018 to dt_output
            self.nb_steps = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'nb_steps')]
            self.nb_satellites = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'nb_sc')]
            self.gps_name = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'gps_name')]
            self.output_file_path_propagator = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'output_file_path_list')]
            self.output_filename_propagator = self.input_variables[find_in_read_input_order_variables(self.order_input_variables, 'output_file_name_list')]

            #Read propagator output files
            self.nb_spec_pts = 4
            self.interpolation_step = 1 # in second, interpolation step of find_specular_points.c (1 s usually)
            self.nb_steps_interpolation = (int)((self.nb_steps-1) * self.dt / self.interpolation_step) #+1 cbv commented the +1 on june 3 2018 because find_specular_points.c stops one second before final_epoch (the way the interpolation works), before wwe had to go to final_epoch (so keep the +1) because satellite files went all the way to final_epoch. But now the satellite files are not read anyoire (sat positions are part of the spec files)

            #Declare necessary variables
            self.lon_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
            self.lat_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
            #self.ecef_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation, 3])
            self.heading_sat = np.zeros([self.nb_satellites, self.nb_steps_interpolation])
            self.time_sat = []
            self.time_sat_datetime = []
            self.lon_spec, self.lat_spec, self.gain_spec, self.distance_spec_to_storm, self.specular_over_storm = np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation]), np.zeros([self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation])

            self.lon_spec_order_report,self.lat_spec_order_report = np.zeros([self.nb_satellites, self.nb_steps_interpolation, self.nb_spec_pts ]), np.zeros([self.nb_satellites, self.nb_steps_interpolation, self.nb_spec_pts ]) # the report is ordered by sat then time.
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


            ### SPECULAR POINTS ### and satellites fter june 3 2018
            for i in range(self.nb_satellites):
                print i, self.nb_satellites-1
                self.time_spec_sublist = []
                self.name_spec_between_list_and_sublist = []

                self.spec_dir = ""
                for j in range(len(self.output_file_path_propagator[i].split('/'))-2):
                    if (j > 0):
                        self.spec_dir = self.spec_dir + "/" + self.output_file_path_propagator[i].split('/')[j]

                self.specular_filepath = self.output_file_path_propagator[i] + 'specular_' + self.output_filename_propagator[i].replace(".txt", ".bin")
                spec_file = open(self.specular_filepath, "rb")
                iGps_temp = spec_file.read(4)
                count = 0
                j = 0
                while iGps_temp != "":
                    iGps = unpack('i', iGps_temp)[0] 
                    iPt = unpack('i', spec_file.read(4))[0]  # 0
                    iPtInner  = ( iPt % (int) (self.dt) );
                    if count == 0:
                        iPtInner_previous = iPtInner
                        lon_sat_this_time = []; lat_sat_this_time = []; heading_sat_this_time = []; lon_spec_this_time = []; lat_spec_this_time = []; gain_spec_this_time = []; prn_this_time= []

                    count = 1

                    if iPtInner != iPtInner_previous: # new time step
                        self.lon_sat[i,j] = lon_sat_this_time[0]
                        if (self.lon_sat[i,j] > 180):
                            self.lon_sat[i,j] = self.lon_sat[i,j] - 360.
                        self.lat_sat[i,j] = lat_sat_this_time[0]
                        self.heading_sat[i,j] = heading_sat_this_time[0]
                        if (i == 0):
                            self.time_sat.append(datetime.strftime(date_now, "%Y-%m-%dT%H:%M:%S"))
                            self.time_sat_datetime.append(date_now)
                        for ispec_now in range(len(lon_spec_this_time)):
                            #print ispec_now,i,j, '|', self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation
                            self.lon_spec[ispec_now,i,j] = lon_spec_this_time[ispec_now]
                            if self.lon_spec[ispec_now,i,j] > 180:
                                self.lon_spec[ispec_now,i,j] = self.lon_spec[ispec_now,i,j] - 360.
                            self.lat_spec[ispec_now,i,j] = lat_spec_this_time[ispec_now]
                            self.gain_spec[ispec_now,i,j] = gain_spec_this_time[ispec_now]
                        self.name_spec_between_list_and_sublist.append(prn_this_time)
                        lon_sat_this_time = []; lat_sat_this_time = []; heading_sat_this_time = []; lon_spec_this_time = []; lat_spec_this_time = []; gain_spec_this_time = []; prn_this_time= []
                        j = j + 1
                    if (iPtInner  == 0):
                        time_ymdhmsm  = np.zeros([7])
                        for sss in range(7):
                            time_ymdhmsm[sss]  = unpack('f', spec_file.read(4))[0]
                        date_now_str_temp  = str((int)(time_ymdhmsm[0] )) + '-' + str((int)(time_ymdhmsm[1])) + '-' + str((int)(time_ymdhmsm[2])) + 'T' + str((int)(time_ymdhmsm[3])) + ':' + str((int)(time_ymdhmsm[4])) + ':' + str((int)(time_ymdhmsm[5]))
                        date_now  = datetime.strptime(date_now_str_temp, "%Y-%m-%dT%H:%M:%S")
                        date_now_save  = date_now
                    else:
                        date_now  = date_now_save + timedelta(seconds  = iPtInner)

                    lon_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) 
                    lat_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) 
                    heading_sat_this_time.append(unpack('f', spec_file.read(4))[0] ) 
                    #ipdb.set_trace()
                    prn_this_time.append((int)(self.gps_name[iGps]))
                    lon_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # 197.093384
                    lat_spec_this_time.append(unpack('f', spec_file.read(4))[0] ) # -28.814657
                    gain_spec_this_time.append(unpack('B', spec_file.read(1))[0] ) # 8
                    iGps_temp = spec_file.read(4)
                    iPtInner_previous = iPtInner
                    #print  date_now, lon_spec, lat_spec, gain_spec, 'PRN_' + gps_name[iGps], 

                # foro the last time step
                self.lon_sat[i,j] = lon_sat_this_time[0]
                if (self.lon_sat[i,j] > 180):
                    self.lon_sat[i,j] = self.lon_sat[i,j] - 360.
                self.lat_sat[i,j] = lat_sat_this_time[0]
                self.heading_sat[i,j] = heading_sat_this_time[0]
                if (i == 0):
                    self.time_sat.append(datetime.strftime(date_now, "%Y-%m-%dT%H:%M:%S"))
                    self.time_sat_datetime.append(date_now)

                for ispec_now in range(len(lon_spec_this_time)):
                    #print ispec_now,i,j, '|', self.nb_spec_pts, self.nb_satellites, self.nb_steps_interpolation
                    self.lon_spec[ispec_now,i,j] = lon_spec_this_time[ispec_now]
                    if self.lon_spec[ispec_now,i,j] > 180:
                        self.lon_spec[ispec_now,i,j] = self.lon_spec[ispec_now,i,j] - 360.
                    self.lat_spec[ispec_now,i,j] = lat_spec_this_time[ispec_now]
                    self.lon_spec_order_report[self.satNameInt[i],j,ispec_now] = self.lon_spec[ispec_now,i,j]
                    self.lat_spec_order_report[self.satNameInt[i],j,ispec_now] = self.lat_spec[ispec_now,i,j]
                    self.gain_spec[ispec_now,i,j] = gain_spec_this_time[ispec_now]
                self.name_spec_between_list_and_sublist.append(prn_this_time)

                spec_file.close()

            # Build the tuples for the visualization of the satellites  - was moved here n june 3 2018 (when specular files are read as nimary not ascii)
                spacecraft = namedtuple('spacecraft',('name',) +  self.point._fields + ('point_plot',) + ('marker_spacecraft',))
                self.spacecraft_list.append(spacecraft)
                name_temp = self.output_filename_propagator[i].replace(".txt","")
                self.spacecraft_list[i].name = name_temp
                # initial position
                #self.spacecraft_list[i].x, self.spacecraft_list[i].y =  m(self.lon_sat[i,self.init_index], self.lat_sat[i,self.init_index])
                self.spacecraft_list[i].marker_spacecraft = 'v'
                # point on the plot
                #self.spacecraft_list[i].point_plot = m.plot([],[],  marker=self.spacecraft_list[i].marker_spacecraft, markersize=10,color = self.satColors[i])[0]


            # Build the tuples for the visualization of the specular points
                for k in range(self.nb_spec_pts):
                    self.specular = namedtuple('specular',('name',) +  self.point._fields  + self.color._fields + ('point_plot',))
                    self.specular_list.append(self.specular)
                    # initial position                    
                    #self.specular_list[k+i*self.nb_spec_pts].x, self.specular_list[k+i*self.nb_spec_pts].y =  m(self.lon_spec[k,i,self.init_index], self.lat_spec[k,i,self.init_index])
                    #self.specular_list[k+i*self.nb_spec_pts].point_plot = m.plot([],[], marker='o', markersize = 2, color = self.satColors[i], fillstyle = 'none', mew = 2)[0]

                self.name_spec.append(self.name_spec_between_list_and_sublist)
                
        else:
            self.dt_sat = 1 # in seconds !!!! has to be the same as dt_output in spock_cygnss_spec_parallel_for_sift
            while new_date <= endDate_date:
                self.time_sat.append(new_date.isoformat())
                new_date = new_date + timedelta(seconds = self.dt_sat)


        # Set filepath for rundir of this simulation        
        self.runFilepath = (self.startDate + "_end_" + self.endDate).replace(":", "_")

        #ipdb.set_trace();
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


        if debug == 1:
            ipdb.set_trace()

        # v1.2
        #        generate_interpolators()
        #ipdb.set_trace()        
        self.interpolate_storms()

#         if debug == 1:
#             ipdb.set_trace()



        self.forecastArtistlist = [[] for x in range(self.nb_storms)]
        self.trackerArtistlist = [[] for x in range(self.nb_storms)]
        self.trajArtistlist = []

        print 'Number of storms:', self.nb_storms
        ###Build Map Time display###
        self.mapTextCorrectionFactor = 2

        self.satColors = ['black', 'blue', 'red', 'mediumorchid', 'dodgerblue', 'magenta', 'darkgreen', 'limegreen'] # !!!! new after june 4 2018
        self.label_arr_conversion = [3, 2, 7, 1, 0, 5, 6, 4]
        self.label_arr = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
#         self.colorByIntensity = {'DB':'deepskyblue', 'TD':'darkturquoise', 'TS':'orange', 'TY':'darkcyan', 'ST':'darkslategray', 'TC':'tomato',
#                                     'HU':'red', 'SD':'chartreuse', 'SS':'forestgreen', 'EX':'purple', 'PT':'limegreen',
#                                     'IN':'blueviolet', 'DS':'darkorchid', 'LO': 'slateblue', 'WV': 'mediumslateblue', 'ET':'blueviolet',
#                                     'XX':'slategray'}
#     #This dictionary is meant to hold a 'lighter' version of the main colors representing each storm type for forecasts in future
#         self.colorByIntensity_faded = {'DB':'skyblue', 'TD':'paleturquoise', 'TS':'navajowhite', 'TY':'lightcyan', 'ST':'slategray', 'TC':'lightsalmon',
#                                     'HU':'pink', 'SD':'palegreen', 'SS':'lightgreen', 'EX':'orchid', 'PT':'yellowgreen',
#                                     'IN':'thistle', 'DS':'orchid', 'LO': 'azure', 'WV': 'lightblue', 'ET':'plum',
#                                     'XX':'whitesmoke'}
        self.colorByIntensity = {'DB':'deepskyblue', 'TD':'darkturquoise', 'TS':'orange', 'TY':'darkcyan', 'ST':'darkslategray', 'TC':'tomato',
                                    'HU':'red', 'SD':'chartreuse', 'SS':'forestgreen', 'EX':'purple', 'PT':'limegreen',
                                    'IN':'blueviolet', 'DS':'darkorchid', 'LO': 'slateblue', 'WV': 'mediumslateblue', 'ET':'blueviolet',
                                    'XX':'slategray'}


        if debug == 1:
            ipdb.set_trace()

#        for self.istorm in range(self.nb_storms):
        self.istorm = 0
        self.bla = 0
        while self.istorm < self.nb_storms:
            # if self.bla == 0:
            #     self.initialize_plot_map()
            self.curSample = self.finalForecasts[self.istorm][0][1]
            self.colors = ['grey', 'dodgerblue', 'blue', 'limegreen', 'yellow', 'red' ,'darkred']#cm.rainbow(np.linspace(0, 1, len(self.finalForecasts[self.istorm])))
            while (datetime.strptime(self.time_sat[self.curSample], "%Y-%m-%dT%H:%M:%S") <= self.finalForecasts[self.istorm][-1][0]):
                self.initialize_plot_map()
                self.plotStorms(self.plotMap)
                self.curSample += self.speedFactor
                if self.curSample >= len(self.time_sat):
                    break
                # if ( self.curSample - self.finalForecasts[self.istorm][0][1] ) / self.speedFactor == 3:
                #     plt.close('all'); sys.exit()
            ani_name = self.rawForecasts[self.istorm][0][0] + self.time_sat[0][:4] + '.mp4' 
            if save_ani == 1:
                os.system('/opt/local/bin/ffmpeg -y -r 3 -i /Users/cbv/cygnss/stormTrackCygnssWebsite/animation_new/temp' + self.rawForecasts[self.istorm][0][0] + self.time_sat[0][:4] + '_%d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p /Users/cbv/cygnss/stormTrackCygnssWebsite/animation_new/' + ani_name)                
                #ipdb.set_trace()
            self.bla = self.bla + 1
            self.istorm = self.istorm + 1
            if (self.istorm == self.nb_storms):
                sys.exit()
    def initialize_plot_map(self):
        ### MAP ###
        dpi = 80.0
        self.height_fig_with_map = 14
        self.width_fig_with_map = 8

        self.fontsize_plot = 15
#v1.2
#        self.fig_with_map = plt.figure(figsize=(self.height_fig_with_map, self.width_fig_with_map))
        self.fig_with_map = plt.figure(num=None, figsize=(self.height_fig_with_map, self.width_fig_with_map),dpi = dpi, facecolor='w', edgecolor='k')
        plt.rc('font', weight='normal') ## make the labels of the ticks in normal
        # Axes
        self.gs_with_map = gridspec.GridSpec(1, 1)
        self.gs_with_map.update(left=0.08, right=0.98, top = 0.94,bottom = 0.1)
        self.ax_map = self.fig_with_map.add_subplot(self.gs_with_map[0, 0])

        ax_title = self.rawForecasts[self.istorm][-1][8].title() + ' ('+ self.rawForecasts[self.istorm][0][0]+ ') - ' + self.time_sat[self.curSample].replace("-","/").replace("T", " ") + ' UTC'
        y_label = 'Latitude '+ u'(\N{DEGREE SIGN})'
        x_label = 'Longitude '+ u'(\N{DEGREE SIGN})'         # the extra \n\n are so that the storm category lggend doesn't get cropped
        self.ax_map.set_title(ax_title, weight = 'normal', fontsize  = (int)(self.fontsize_plot*1.1), y = 1.)
        self.ax_map.set_ylabel(y_label, weight = 'normal', fontsize  = self.fontsize_plot)
        self.ax_map.set_xlabel(x_label, weight = 'normal', fontsize  = self.fontsize_plot)

        [i.set_linewidth(1) for i in self.ax_map.spines.itervalues()] # change the width of the frame of the figure
        self.ax_map.tick_params(axis='both', which='major', labelsize=self.fontsize_plot, size = 7, width = 1, pad = 7) 
        plt.rc('font', weight='normal') ## make the labels of the ticks in normal

        
        #This is where we set map params for different zoom settings
        self.zoom_names = ['Globe', 'NWAtlantic', 'SWAtlantic', 'NEPacific', 'SEPacific', 'All Atlantic', 'Asia Pacific', 'Gulf of Mexico', 'Indian Ocean']
        
        self.padding_lon = np.max([5., self.size_max_storm[self.istorm] / 110.+ 1.])
        self.padding_lat = np.max([5., self.size_max_storm[self.istorm] / 110. + 1.])
        if self.lon_min_storm[self.istorm] >= -179.9 + self.padding_lon:
            self.min_lon = self.lon_min_storm[self.istorm] - self.padding_lon
        else:
            self.min_lon = -179.9
        if self.lon_max_storm[self.istorm] <= 179.9 - self.padding_lon:
            self.max_lon = self.lon_max_storm[self.istorm] + self.padding_lon
        else:
            self.max_lon = 179.9
        if self.lat_min_storm[self.istorm] >= -89.9 + self.padding_lat:
            self.min_lat = self.lat_min_storm[self.istorm] - self.padding_lat
        else:
            self.min_lat = -179.9
        if self.lat_max_storm[self.istorm] <= 89.9 - self.padding_lat:
            self.max_lat = self.lat_max_storm[self.istorm] + self.padding_lat
        else:
            self.max_lat = 179.9


        self.step_lon = (int) ( (self.max_lon - self.min_lon)/8. * 10. ) / 10.
        self.step_lat = (int) ( (self.max_lat - self.min_lat)/8. * 10. ) / 10.

        # #Read in custom zoom levels from ../input_sift/user_zoom.txt
        # cNames, cMin_lons, cMax_lons, cMin_lats, cMax_lats, cStep_lons, cStep_lats = readZoom(self.userZoom_filepath)

        # #Append info to appropriate list above
        # self.zoom_names += cNames
        # self.min_lon += cMin_lons
        # self.max_lon += cMax_lons
        # self.min_lat += cMin_lats
        # self.max_lat += cMax_lats
        # self.step_lon += cStep_lons
        # self.step_lat += cStep_lats

        self.build_plot_map(self.mapRef)
    def build_plot_map(self, mapRef):
        #Del old keys to remove old map lines; must use del here not clear
        for key in self.meridians.keys():
            del self.meridians[key]
        for key in self.parallels.keys():
            del self.parallels[key]
        #ipdb.set_trace()
        self.array_lon_temp = [ ss for ss in np.arange(self.min_lon, self.max_lon,self.step_lon) ]
        self.array_lat_temp = [ ss for ss in np.arange(self.min_lat, self.max_lat,self.step_lat) ]
        self.array_lon = []
        self.array_lat = []
        for i in range(len(self.array_lon_temp)):
            if (self.array_lon_temp[i] < 0):
                self.array_lon.append( format(np.abs(self.array_lon_temp[i]), ".1f") + 'W' )
            else:
                self.array_lon.append( format(self.array_lon_temp[i], ".1f") + 'S' )
        for i in range(len(self.array_lat_temp)):
            if (self.array_lat_temp[i] < 0):
                self.array_lat.append( format(np.abs(self.array_lat_temp[i]), ".1f") + 'S' )
            else:
                self.array_lat.append( format(self.array_lat_temp[i], ".1f") + 'N' )

                #ipdb.set_trace()

        self.ax_map.xaxis.set_major_locator(FixedLocator(np.arange(self.min_lon, self.max_lon+1, self.step_lon)))
        self.ax_map.yaxis.set_major_locator(FixedLocator(np.arange(self.min_lat, self.max_lat+1, self.step_lat)))

        self.ax_map.set_xticklabels(self.array_lon, color = 'k')
        self.ax_map.set_yticklabels(self.array_lat, color = 'k')

        self.plotMap = Basemap( projection       = 'cyl',
                                llcrnrlon        = self.min_lon,#-97.4,#self.min_lon,#-180,#self.min_lon , #Lower Left  CoRNeR Longitude
                                urcrnrlon        = self.max_lon,#-19.8,#self.max_lon  ,#180,#self.max_lon  , #Upper Right CoRNeR Longitude
                                llcrnrlat        = self.min_lat,#8.8,#self.min_lat,#-40,#,self.min_lat  , #Lower Left  CoRNeR Latitude
                                urcrnrlat        = self.max_lat,#44.0,#self.max_lat,#40,#self.max_lat,   #Upper Right CoRNeR Latitude
                 resolution       = 'l'  ,
                 suppress_ticks   = False,
                 ax = self.ax_map,
                 )

        alpha = 0.5
        color_continents = [65,105,225,alpha*256]
        color_continents = np.array(color_continents) / 256.
        color_water  = [100,149,237,alpha*256]
        color_water = np.array(color_water) / 256.

        self.plotMap.fillcontinents(color=tuple(color_continents),lake_color=tuple(color_water))
        self.plotMap.drawmapboundary(fill_color=tuple(color_water))

        mapRef = self.plotMap.drawcoastlines(linewidth=0.7, color='blue')
        # self.ax_map.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # self.ax_map.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        #self.meridians = self.plotMap.drawmeridians(np.arange(self.min_lon, self.max_lon,self.step_lon))
        #self.parallels = self.plotMap.drawparallels(np.arange(self.min_lat, self.max_lat,self.step_lat))


    def build_raw_storms(self):
        #if debug == 1:
            #ipdb.set_trace()
        #This list holds tuples of raw forecast data
        self.rawForecasts = []
        #This list holds tuples of (stormName, warningTime) for map
        self.stormInfo = []

        #Need start/end times when deciding indices
        startTime = datetime.strptime(self.forecastStart, "%Y-%m-%dT%H:%M:%S")
        endTime = datetime.strptime(self.forecastEnd, "%Y-%m-%dT%H:%M:%S")

        stormFiles = self.list_storm
        all_date_storm = []
        self.lon_max_storm = []
        self.lon_min_storm = []
        self.lat_max_storm = []
        self.lat_min_storm = []
        self.size_max_storm = []

        for file in range(len(stormFiles)):
            commonName_good = ''
            STORM_FILEPATH =  'data/' + stormFiles[file]
            stormFile = open(STORM_FILEPATH, 'r')
            #print self.storm_dir + '/' + stormFiles[file]
            forecastLines_old = stormFile.readlines()

            forecastLines = []
            with gzip.open(STORM_FILEPATH,'r') as fin:        
                for line in fin:        
                    forecastLines.append(line)
            if len(forecastLines) > 1: #Don't include single forecasts
                #tmpForecasts allows us to append on a per-storm basis
                tmpForecasts = []
                #tmpTimes allows us to only look at 34 kt radii for each storm
                tmpTimes = []
                #tmpLons lets us address duplicate longitudes
                tmpLons = []
                if(STORM_FILEPATH[-7:] == '.dat.gz'): # observations
                    #Need first and last cast to determine if storm is in interval
                    firstCast = [x.strip() for x in forecastLines[0].split(',')]
                    firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                    lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                    lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H")

                    #TEMP FIX - Checks if any forecast for storm falls in analysis interval
                    #TODO: Fully iron out best way to weed out storms and reduce code duplication here
                    inIntervalFlag = 0
                    for i in range(len(forecastLines)):
                        cast = [x.strip() for x in forecastLines[i].split(',')]
                        castTime = datetime.strptime(cast[2], "%Y%m%d%H")
                        if castTime >= startTime and castTime < endTime:
                            inIntervalFlag = 1
                            break #Stop once we know at least one forecast is in interval

                    if inIntervalFlag:
                        logging.debug("Storm %s in analysis interval", firstCast[0])
                        #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                        self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                        lon_max = -1000.
                        lon_min = 1000.
                        lat_max = -1000.
                        lat_min = 1000.
                        size_max = -1000.
                        #For each line in forecast file
                        for iStorm in range(len(forecastLines)):
                            noaa = [x.strip() for x in forecastLines[iStorm].split(',')]
                            #Build basic storm parameters
                            if noaa[2] not in tmpTimes: #So we only read forecast one time 
                                #Timestamp
                                ts = datetime.strptime(noaa[2], "%Y%m%d%H") 

                                #Index
                                if ts < startTime:
                                    sIdx = -2 #Arbitrary flag 
                                elif ts > endTime:
                                    sIdx = -1 #Arbitrary flag
                                else:
                                    sIdx = 0#int(self.time_sat.index(ts.isoformat()))
                                lat = noaa[6]
                                if lat[-1:] == 'N':
                                    lat = lat.replace("N","")
                                    lat = float(lat)/10.0
                                else:
                                    lat = lat.replace("S","")
                                    lat = -1*float(lat)/10.0
                                

                                #Lon
                                lon = noaa[7]
                                if lon[-1:] =='E':
                                    lon = lon.replace("E","")
                                    lon = float(lon)/10.0
                                else:
                                    lon = lon.replace("W","")
                                    lon = -1*float(lon)/10.0

                               #Must be no duplicates to satisfy spline condition
                                while lon in tmpLons:
                                    logging.debug("Found %s in tmp lons, altering.", str(lon))
                                    lon -= 0.0001 # used to be 0.3 
                                tmpLons.append(lon)

                                if lat >= lat_max:
                                    lat_max = lat
                                if lat <= lat_min:
                                    lat_min = lat

                                if lon >= lon_max:
                                    lon_max = lon
                                if lon <= lon_min:
                                    lon_min = lon

                                # Sustainded wind speed
                                windSpeedValue = float(noaa[8]) # in kts

                                #Storm type
                                typ = noaa[10]

                                #Size
                                windRadii = [int(noaa[13]), int(noaa[14]), int(noaa[15]), int(noaa[16])]
                                size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding #!!!! padding whould be + 30

                                if size >= size_max:
                                    size_max = size

                                tmpTimes.append(noaa[3]) 
                                name = noaa[0]+noaa[1]
                                if len(noaa) > 26:
                                    if noaa[27] != 'INVEST':
                                        commonName_good = noaa[27]
                                commonName = commonName_good
                                newForecast = (name, ts,  lat, lon, typ, size, sIdx, windSpeedValue, commonName)
                                
                                tmpForecasts.append(newForecast)
                                all_date_storm.append(ts)
                        if tmpForecasts:
                            #Create list of lists of tuples (one list for each active storm)
                            self.rawForecasts.append(tmpForecasts)
                            self.lon_max_storm.append(lon_max)
                            self.lon_min_storm.append(lon_min)
                            self.lat_max_storm.append(lat_max)
                            self.lat_min_storm.append(lat_min)
                            self.size_max_storm.append(size_max)                
                elif(STORM_FILEPATH[-3:] == 'txt' and os.stat(STORM_FILEPATH).st_size != 0):
                    #Need first and last cast to determine if storm is in interval
                    firstCast = [x.strip() for x in forecastLines[0].split(',')]
                    firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                    lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                    lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H") + timedelta(hours=int(lastCast[3]))

                    #TEMP FIX - Checks if any forecast for storm falls in analysis interval
                    #TODO: Fully iron out best way to weed out storms and reduce code duplication here
                    inIntervalFlag = 0
                    for i in range(len(forecastLines)):
                        cast = [x.strip() for x in forecastLines[i].split(',')]
                        castTime = datetime.strptime(cast[2], "%Y%m%d%H") + timedelta(hours=int(cast[3]))
                        if castTime >= startTime and castTime < endTime:
                            inIntervalFlag = 1
                            break #Stop once we know at least one forecast is in interval

                    if inIntervalFlag:
                        logging.debug("Storm %s in analysis interval", firstCast[0])
                        #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                        self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                        lon_max = -1000.
                        lon_min = 1000.
                        lat_max = -1000.
                        lat_min = 1000.
                        size_max = -1000.
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
                                    sIdx = 0#int(self.time_sat.index(ts.isoformat()))

                                #Lat
                                lat = float(noaa[4])
                                
                                #Lon
                                lon = float(noaa[5])

                                #Must be no duplicates to satisfy spline condition
                                while lon in tmpLons:
                                    logging.debug("Found %s in tmp lons, altering.", str(lon))
                                    lon -= 0.0001 # used to be 0.3
                                tmpLons.append(lon)
                                if lat >= lat_max:
                                    lat_max = lat
                                if lat <= lat_min:
                                    lat_min = lat

                                if lon >= lon_max:
                                    lon_max = lon
                                if lon <= lon_min:
                                    lon_min = lon

                                #Storm type
                                typ = noaa[9]

                                #Size
                                windRadii = [int(noaa[12]), int(noaa[13]), int(noaa[14]), int(noaa[15])]
                                size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding #!!!! padding whould be + 30

                                if size >= size_max:
                                    size_max = size

                                tmpTimes.append(noaa[3]) 
                                name = noaa[0]+noaa[1]
                                newForecast = (name, ts,  lat, lon, typ, size, sIdx)

                                tmpForecasts.append(newForecast)
                                all_date_storm.append(ts)
                        if tmpForecasts:
                            #Create list of lists of tuples (one list for each active storm)
                            self.rawForecasts.append(tmpForecasts)
                            self.lon_max_storm.append(lon_max)
                            self.lon_min_storm.append(lon_min)
                            self.lat_max_storm.append(lat_max)
                            self.lat_min_storm.append(lat_min)
                            self.size_max_storm.append(size_max)
                elif(STORM_FILEPATH[-3:] == 'fst' and os.stat(STORM_FILEPATH).st_size != 0):
                    firstCast = [x.strip() for x in forecastLines[0].split(',')]
                    firstCastTime = datetime.strptime(firstCast[2], "%Y%m%d%H")
                    lastCast = [x.strip() for x in forecastLines[-1].split(',')]
                    lastCastTime = datetime.strptime(lastCast[2], "%Y%m%d%H") + timedelta(hours=int(lastCast[5]))
                    
                    inIntervalFlag = 0
                    for i in range(len(forecastLines)):
                        cast = [x.strip() for x in forecastLines[i].split(',')]
                        castTime = datetime.strptime(cast[2], "%Y%m%d%H") + timedelta(hours=int(cast[3]))
                        if castTime >= startTime and castTime < endTime:
                            inIntervalFlag = 1
                            break #Stop once we know at least one forecast is in interval

                    #Check if storm falls at all in display interval
                    #TEMP FIX
                    #TODO: Fully iron out best way to weed out storms and reduce code duplication here
                    if inIntervalFlag:
                        #Store name and warning time for on-map display; store warning time as string since no work will be done with it
                        self.stormInfo.append((firstCast[0]+ firstCast[1], firstCast[2]))
                        lon_max = -1000.
                        lon_min = 1000.
                        lat_max = -1000.
                        lat_min = 1000.
                        size_max = -1000.
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
                                    sIdx = int(0)#!!!!!self.startDate)

                                #Lat
                                lat = jtwc[6]
                                if lat[-1:] == 'N':
                                    lat = lat.replace("N","")
                                    lat = float(lat)/10.0
                                else:
                                    lat = lat.replace("S","")
                                    lat = -1*float(lat)/10.0
                                

                                #Lon
                                lon = jtwc[7]
                                if lon[-1:] =='E':
                                    lon = lon.replace("E","")
                                    lon = float(lon)/10.0
                                else:
                                    lon = lon.replace("W","")
                                    lon = -1*float(lon)/10.0

                               #Must be no duplicates to satisfy spline condition
                                while lon in tmpLons:
                                    logging.debug("Found %s in tmp lons, altering.", str(lon))
                                    lon -= 0.0001 # used to be 0.3 
                                tmpLons.append(lon)
                                
                                if lat >= lat_max:
                                    lat_max = lat
                                if lat <= lat_min:
                                    lat_min = lat

                                if lon >= lon_max:
                                    lon_max = lon
                                if lon <= lon_min:
                                    lon_min = lon

                                #Storm type
                                typ = jtwc[10]

                                #Size
                                windRadii = [int(jtwc[13]), int(jtwc[14]), int(jtwc[15]), int(jtwc[16])]
                                size = float(max(windRadii) *1.852 + 30) #Conversion from naut. mi. to km + padding

                                if size >= size_max:
                                    size_max = size
                                
                                tmpTimes.append(jtwc[5]) 
                                name = jtwc[0] + jtwc[1]
                                newForecast = (name, ts, lat, lon, typ, size, sIdx)
                                
                                tmpForecasts.append(newForecast)
                                all_date_storm.append(ts)
                        if tmpForecasts:
                            #Create list of lists of tuples (one list for each active storm)
                            self.rawForecasts.append(tmpForecasts)
                            self.lon_max_storm.append(lon_max)
                            self.lon_min_storm.append(lon_min)
                            self.lat_max_storm.append(lat_max)
                            self.lat_min_storm.append(lat_min)
                            self.size_max_storm.append(size_max)

        #if debug == 1:
#             ipdb.set_trace()

        self.nb_storms = len(self.rawForecasts)
#V1.2
        self.all_date_storm_arrange = np.sort(np.array(all_date_storm))

    def do_interpolate( self, system ):
#         if debug == 1:                                                                                                           
#             ipdb.set_trace()                                                                                                    

        from scipy import interpolate
        lats = [x[2] for x in system]
        lons = [x[3] for x in system]
        if len( lats ) == 3:
            tck,u=interpolate.splprep([lons,lats],s=0.0, k=2)
        elif len( lats ) == 2:
            tck,u=interpolate.splprep([lons,lats],s=0.0, k=1)
        else:
            tck,u=interpolate.splprep([lons,lats],s=0.0)
        # 40000 > num values rounded to tenths place between -180 and 180
        newLons,newLats = interpolate.splev(np.linspace(0,1,40000),tck) 
        d = {}
        for i,newLon in enumerate(newLons):
            d[newLon] = newLats[i]
#         if debug == 1:                                                                                                           
#             ipdb.set_trace()                                                                                                    

        return( d )
        
    def interpolate_storms(self):
#         if debug == 1:
        #ipdb.set_trace()

        finalTimestamps = [[] for x in range(self.nb_storms)]
        finalIdxs = [[] for x in range(self.nb_storms)]
        finalLats = [[] for x in range(self.nb_storms)]
        finalLons = [[] for x in range(self.nb_storms)]
        finalTypes = [[] for x in range(self.nb_storms)]
        finalSizes = [[] for x in range(self.nb_storms)]
        finalWindSpeeds = [[] for x in range(self.nb_storms)]
        
        timeSkip = self.speedFactor#np.max([720, self.speedFactor]) #12 minutes in seconds
        startTime = datetime.strptime(self.forecastStart, "%Y-%m-%dT%H:%M:%S")
        endTime = datetime.strptime(self.forecastEnd, "%Y-%m-%dT%H:%M:%S")

        systemIdx = 0
        if debug == 1:
            ipdb.set_trace()
        
        for system in self.rawForecasts:
            if len(system) > 1:
                for i in range(len(system) - 1):
                    #Figure out number of steps in this interval
                    intervals = int(divmod((system[i+1][1] - system[i][1]).total_seconds(), timeSkip)[0])
                    for j in range(intervals):
                        #Build new timestamps
                        fts = system[i][1] + j*timedelta(seconds=timeSkip)
#                         if str(fts) == '2017-09-10 06:00:00':
#                             ipdb.set_trace()
                        #print fts, system[i][1]
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
                    finalLons[systemIdx].extend(newLons)


                    #Generate new types
                    finalTypes[systemIdx].extend([system[i][4]] * intervals)


                    #Generate new sizes
                    finalSizes[systemIdx].extend(np.linspace(system[i][5], system[i+1][5], intervals, endpoint=False))

                    #Generate new wind speeds
                    finalWindSpeeds[systemIdx].extend(np.linspace(system[i][7], system[i+1][7], intervals, endpoint=False))


                # last time step
                #Figure out number of steps in this interval
                i = -1
                #Build new timestamps
                fts = system[i][1]
                finalTimestamps[systemIdx].append(fts)                        
                #Retrieve new indices
                if fts < startTime:
                    finalIdxs[systemIdx].append(-2) #Arbitrary flag 
                elif fts > endTime:
                    finalIdxs[systemIdx].append(-1) #Arbitrary flag
                else:
                    finalIdxs[systemIdx].append(int(self.time_sat.index(fts.isoformat())))
                #Generate new lons
                finalLons[systemIdx].extend([system[i][3]])
                #Generate new types
                finalTypes[systemIdx].extend([system[i][4]])
                #Generate new sizes
                finalSizes[systemIdx].extend([system[i][5]])
                #Generate new wind speeds
                finalWindSpeeds[systemIdx].extend([system[i][7]])
            

            #Generate new lats from interpolator
            finalLats[systemIdx] = []
            dLatLons = self.do_interpolate( system )
            for finalLon in finalLons[systemIdx]:
                lookup = dLatLons.get( finalLon, -9999 )
                if lookup == -9999:
                    arr = np.array(dLatLons.keys())
                    closest = arr[ np.abs( arr-finalLon ).argmin() ]
                    finalLats[systemIdx].append( dLatLons[closest] )
                else:
                    finalLats[systemIdx].append( dLatLons[finalLon] )
#             if debug == 1:
#                 ipdb.set_trace()
 
# Adding code to sort finalLats/Lons by euclidian distance so storms update smoothly
            combined = zip( finalLons[systemIdx],finalLats[systemIdx] )
            whatsLeft = combined[:]
            myPt = whatsLeft[0]
            final = []
            final.append( myPt)
            whatsLeft = [ i for i in whatsLeft if i!=myPt ]
            while whatsLeft:
                nextPt = whatsLeft[ distance.cdist( [myPt],np.array(whatsLeft) ).argmin()]
                final.append( nextPt )
                whatsLeft = [ i for i in whatsLeft if i!=nextPt ]
                myPt = nextPt
            newLons,newLats = zip(*final)
            finalLons[systemIdx] = list(newLons)
            finalLats[systemIdx] = list(newLats)
                
                        
            systemIdx += 1
#             if debug == 1:
#                 ipdb.set_trace()

#         if debug == 1:
#             ipdb.set_trace()

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
            self.finalForecasts.append(zip(finalTimestamps[i], finalIdxs[i], finalLats[i], finalLons[i], finalTypes[i], finalSizes[i], finalWindSpeeds[i]))

    def updateFirstStorm(self, storm):
        if len(self.forecastArtistlist[storm]) > 0:
            self.forecastArtistlist[storm][0].set_zorder(2)
            self.forecastArtistlist[storm][0].set_facecolor(self.colorByIntensity[self.finalForecasts[storm][self.curStorm[storm]][4]])
            self.forecastArtistlist[storm][0].set_edgecolor('black')

    def find_le(self, inList, x):
#         #'Find rightmost value less than or equal to x'
#         i = bisect_right(inList, x)
        if len(np.where(np.array(inList) > x)[0]) > 0:
            i = np.where(np.array(inList) > x)[0][0]
            if inList[i-1] >=0:
                return inList[i-1], 1
            else:
                return inList[0], 0
        elif (inList[-1] == x):
            return inList[-1], 1
        else:
            return inList[0], 0

    def updateCurStorm(self):
        self.validFlags = []
        idxs = [x[1] for x in self.finalForecasts[self.istorm]]
        logging.debug("idxs: %s", str(idxs))
        stormIdx, validFlag = self.find_le(idxs, self.curSample)
        self.validFlags.append(validFlag)
        self.curStorm[self.istorm] = idxs.index(stormIdx)
        logging.debug("curStorm = %s",self.curStorm[self.istorm])
        if validFlag != 1:
            ipdb.set_trace()
    def plotStorms(self, plotMap):
        self.activeStorms = []
        #This save and then repass in forwardloopreset seems unnecessary - improve?
        self.plotMap = plotMap
        first_in_past = True

        self.updateCurStorm()
        #print self.validFlags[self.istorm]
        #print self.time_sat[self.curSample], self.validFlags[self.istorm]
        #idx = self.curStorm[self.istorm]
        #print self.time_sat[self.curSample], self.validFlags[self.istorm]
        if self.validFlags[-1]:
            idx = self.curStorm[self.istorm]
            #print self.time_sat[self.curSample],  self.finalForecasts[self.istorm][idx][0].isoformat(), self.finalForecasts[self.istorm][idx][2], self.finalForecasts[self.istorm][idx][3],self.finalForecasts[self.istorm][idx][4], self.finalForecasts[self.istorm][idx][6], self.validFlags[self.istorm]

            #for ll in range(len(self.finalForecasts[self.istorm])): # previous locaiton of storm
                #self.forecastArtistlist[self.istorm].append(plotMap.tissot(self.finalForecasts[self.istorm][ll][3], self.finalForecasts[self.istorm][ll][2], radius_for_tissot(self.finalForecasts[self.istorm][ll][5]), 256, facecolor=self.colorByIntensity_faded[self.finalForecasts[self.istorm][ll][4]], edgecol®or='none'))
            for ll in range(len(self.finalForecasts[self.istorm])): # previous locaiton of storm

                windSpeedValue = np.array(self.finalForecasts[self.istorm])[ll,6]
                ## catagoery of storm according to the Saffir-Simpson scale (https://www.nhc.noaa.gov/aboutsshws.php?)
                if  (windSpeedValue <= 34 ):
                    windSpeed = -1 # Tropical depression
                elif ( (windSpeedValue >= 35 ) & (windSpeedValue <= 63)):
                    windSpeed = 0 # Tropical storm
                elif ( (windSpeedValue >= 64 ) & (windSpeedValue <= 82)):
                    windSpeed = 1 # Saffir-Simpson scale starts here (hurricane)
                elif ( (windSpeedValue >= 83) & (windSpeedValue <= 95) ):
                    windSpeed = 2
                elif ( (windSpeedValue >= 96) & (windSpeedValue <= 112) ):
                    windSpeed = 3
                elif ( (windSpeedValue >= 113) & (windSpeedValue <= 136) ):
                    windSpeed = 4
                else:
                    windSpeed = 5

                x, y =  plotMap(np.array(self.finalForecasts[self.istorm])[ll,3],np.array(self.finalForecasts[self.istorm])[ll, 2])
                plotMap.scatter(x, y,  marker='o', color = self.colors[windSpeed+1], s = 100, zorder = 5)

            windSpeedValue = np.array(self.finalForecasts[self.istorm])[idx,6]
            ## catagoery of storm according to the Saffir-Simpson scale (https://www.nhc.noaa.gov/aboutsshws.php?)
            if  (windSpeedValue <= 34 ):
                windSpeed = -1 # Tropical depression
            elif ( (windSpeedValue >= 35 ) & (windSpeedValue <= 63)):
                windSpeed = 0 # Tropical storm
            elif ( (windSpeedValue >= 64 ) & (windSpeedValue <= 82)):
                windSpeed = 1 # Saffir-Simpson scale starts here (hurricane)
            elif ( (windSpeedValue >= 83) & (windSpeedValue <= 95) ):
                windSpeed = 2
            elif ( (windSpeedValue >= 96) & (windSpeedValue <= 112) ):
                windSpeed = 3
            elif ( (windSpeedValue >= 113) & (windSpeedValue <= 136) ):
                windSpeed = 4
            else:
                windSpeed = 5

            x, y =  plotMap(np.array(self.finalForecasts[self.istorm])[idx,3],np.array(self.finalForecasts[self.istorm])[idx, 2])
            plotMap.scatter(x, y,  marker='o', color = self.colors[windSpeed+1], s = 300, zorder = 5)

#             handles_arr = [mpatches.Patch(color='lightblue', label='')]
#             legend = self.ax_map.legend( loc='center',  bbox_to_anchor=(0.5, -0.144), fontsize = self.fontsize_plot, handles=handles_arr, ncol=2, frameon=False)#    legend = ax_2d_big.legend(loc='center l
            label_storm = ['TD\n            ','TS\n            ','1\n            ','2\n            ','3\n            ','4\n            ', '5\n            ']
            handles_arr = []
            for icat in range(len(self.colors)):
                handles_arr.append(mpatches.Patch(color=self.colors[icat], label=label_storm[icat]))
            self.legend_storm = self.ax_map.legend( loc='center left',  bbox_to_anchor=(0, -0.16), fontsize = self.fontsize_plot, handles=handles_arr, ncol=len(self.colors), frameon=False, columnspacing = -4.45, title = 'Storm Category')
            self.legend_storm.get_title().set_fontsize(str(self.fontsize_plot))
            self.legend_storm.get_title().set_position((-37, 0))
            icat = -1
            for txt in self.legend_storm.get_texts():
                icat = icat + 1
                if icat == 0: # TD
                    txt.set_x(-40);txt.set_y(5); txt.set_color(self.colors[icat])
                elif icat == 1: # TS
                    txt.set_x(-38);txt.set_y(5); txt.set_color(self.colors[icat])
                else:
                    txt.set_x(-35);txt.set_y(5); txt.set_color(self.colors[icat])
            #self.ax_map.text(0.038,-0.19,'Storm Category',transform = self.ax_map.transAxes, fontsize = self.fontsize_plot, horizontalalignment = 'left');


#             self.forecastArtistlist[self.istorm].append(plotMap.tissot(self.finalForecasts[self.istorm][idx][3], self.finalForecasts[self.istorm][idx][2], radius_for_tissot(self.finalForecasts[self.istorm][idx][5]), 256, facecolor=self.colorByIntensity_faded[self.finalForecasts[self.istorm][idx][4]], edgecolor=self.colorByIntensity[self.finalForecasts[self.istorm][idx][4]] , linewidth = 4)) # current location
#             print self.finalForecasts[self.istorm][idx][3], self.finalForecasts[self.istorm][idx][2], self.finalForecasts[self.istorm][idx][0].isoformat(), self.time_sat[self.curSample]
#        self.ax_map.text(0.5,0.02,self.time_sat[self.curSample].replace("-","/").replace("T", " "), fontsize = self.fontsize_plot, transform = self.ax_map.transAxes, horizontalalignment = 'center', verticalalignment = 'center')

        if show_spec == 1:
            for isc_temp in range(self.nb_satellites):
                isc = self.label_arr_conversion[isc_temp]
                self.specular_list[isc_temp].x, self.specular_list[isc_temp].y =  plotMap(self.lon_spec[:,isc,self.curSample-3*3600:self.curSample], self.lat_spec[:,isc,self.curSample-3*3600:self.curSample])
                self.specular_list[isc_temp].point_plot = plotMap.scatter(self.specular_list[isc_temp].x, self.specular_list[isc_temp].y,  marker='o', color = self.satColors[isc], s = 10, zorder = 4,alpha = 0.05, label = self.label_arr[isc])
            for isc_temp in range(self.nb_satellites):
                isc = self.label_arr_conversion[isc_temp]
                self.specular_list[isc_temp].x, self.specular_list[isc_temp].y =  plotMap(self.lon_spec[:,isc,self.curSample], self.lat_spec[:,isc,self.curSample])
                self.specular_list[isc_temp].point_plot = plotMap.scatter(self.specular_list[isc_temp].x, self.specular_list[isc_temp].y,  marker='o', color = self.satColors[isc], s = 10, zorder = 4)
                #self.distance_sp_center_storm = np.linalg.norm

            #self.legend_sc = self.ax_map.legend(loc='center left', bbox_to_anchor=(1, 0.5), numpoints = 1,  title="", fontsize = self.fontsize_plot,handlelength=0, handletextpad=0)
                self.legend_sc = self.ax_map.legend(loc='center right', bbox_to_anchor=(1, -0.17), numpoints = 1,  title="", fontsize = self.fontsize_plot,handlelength=0, handletextpad=0, ncol = 3)
            for isc_temp in range(len(self.legend_sc.get_texts())):
                isc = self.label_arr_conversion[isc_temp]
                self.legend_sc.get_texts()[isc_temp].set_color(self.satColors[isc]) # set the label the same color as the plot
                self.legend_sc.legendHandles[isc_temp].set_visible(False) # hide the line in the label

            self.ax_map.add_artist(self.legend_storm)
#             all_legend = (self.legend_sc,) 
#             all_legend = all_legend + (self.legend_storm,)
#            self.ax_map.set_position([0.1,0.1,0.5,0.8])
        if show_ani == 1:
            self.fig_with_map.canvas.draw() # if comment then ani runs faster but when try to quit with ctrl-c it doesn't work and need to ctrl-z
            self.fig_with_map.canvas.flush_events()
        self.stormsShowed = True
        if save_ani == 1:
            self.fig_save_name = '/Users/cbv/cygnss/stormTrackCygnssWebsite/animation_new/temp' + self.rawForecasts[self.istorm][0][0]  + self.time_sat[0][:4] + '_' + str( ( self.curSample - self.finalForecasts[self.istorm][0][1] ) / self.speedFactor) + '.png'
            self.fig_with_map.subplots_adjust(bottom=0.8)
            self.fig_with_map.savefig(self.fig_save_name, facecolor= self.fig_with_map.get_facecolor(), edgecolor='none', bbox_inches='tight',bbox_extra_artists=(self.legend_sc,))
            


if __name__ == '__main__':
    interact_here = interact()
    

