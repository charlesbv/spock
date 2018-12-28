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
# this script creates a fake swpc f107 ap file (fake because we make it ourselves, not swpc) with a constant ap and f107. This was sued in collision avoidance studies

from datetime import datetime, timedelta

def spock_swpc_f107_ap_file(f107_ap_pred_filename, f107, ap, start_date): # start_date has to be "%y-%m-%d" - ex: '17-04-26'
    
# f107_ap_pred_filename = 'test.txt'
# f107 = 100
# ap = 15
# start_date = '20-12-01'
    start_date_datetime = datetime.strptime(start_date, "%y-%m-%d")
    nb_day = 45 # !!!!! can't change that (SpOCK expects a certain format for this file)
    date = []
    for iday in range(nb_day):
        date_pred = start_date_datetime + timedelta(days = iday)
        date_pred_str = datetime.strftime(date_pred, "%y-%m-%d")

        year = date_pred_str.split('-')[0]
        day = date_pred_str.split('-')[2]
        month_nb = (int)(date_pred_str.split('-')[1])
        # convert month number to str
        if (month_nb == 1):
            month = 'Jan'
        elif (month_nb == 2):
            month = 'Feb'
        elif (month_nb == 3):
            month = 'Mar'
        elif (month_nb == 4):
            month = 'Apr'
        elif (month_nb == 5):
            month = 'May'
        elif (month_nb == 6):
            month = 'Jun'
        elif (month_nb == 7):
            month = 'Jul'
        elif (month_nb == 8):
            month = 'Aug'
        elif (month_nb == 9):
            month = 'Sep'
        elif (month_nb == 10):
            month = 'Oct'
        elif (month_nb == 11):
            month = 'Nov'
        elif (month_nb == 12):
            month = 'Dec'
        else:
            print "***! Problem with the date. The program will stop.!***";raise Exception

        date.append( day + month + year )


    f107_ap_pred_file = open(f107_ap_pred_filename, "w")

    #header
    print >> f107_ap_pred_file, ':Product: 45 Day AP Forecast  45DF.txt\n:Issued: 2017 Nov 30 2120 UTC\n# Prepared by the U.S. Air Force.\n# Retransmitted by the Dept. of Commerce, NOAA, Space Weather Prediction Center\n# Please send comments and suggestions to SWPC.Webmaster@noaa.gov\n#\n#\n#          45-Day AP and F10.7cm Flux Forecast\n#-------------------------------------------------------------'
    nb_line = 9 # can't cahnge that
    nb_col = 5 # can't cahnge that

    # ap prediction
    print >> f107_ap_pred_file, '45-DAY AP FORECAST'
    iday = 0
    for iline in range(nb_line):
        print >> f107_ap_pred_file, date[iday] + ' ' + str(ap).zfill(3) + ' ' + date[iday+1] + ' ' + str(ap).zfill(3) + ' ' + date[iday+2] + ' ' + str(ap).zfill(3) + ' ' +date[iday+3] + ' ' + str(ap).zfill(3) + ' '  +date[iday+4] + ' ' + str(ap).zfill(3)
        iday = iday + 5
    # f107 prediction
    print >> f107_ap_pred_file, '45-DAY F10.7 CM FLUX FORECAST'
    iday = 0
    for iline in range(nb_line):
        print >> f107_ap_pred_file, date[iday] + ' ' + str(f107).zfill(3) + ' ' + date[iday+1] + ' ' + str(f107).zfill(3) + ' ' + date[iday+2] + ' ' + str(f107).zfill(3) + ' ' +date[iday+3] + ' ' + str(f107).zfill(3) + ' '  +date[iday+4]+ ' ' + str(f107).zfill(3)
        iday = iday + 5

    print >> f107_ap_pred_file, 'FORECASTER:  MAXWELL / AGUILERA\n99999'


    f107_ap_pred_file.close()


    return 0
