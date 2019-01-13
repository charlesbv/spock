# this script converts the date from a TLE into a datetime object

from datetime import datetime, timedelta
import numpy as np

def convert_tle_date_to_date(date_tle):
    #date_tle = '17022.87851213'
    date = datetime.strptime(date_tle.split('.')[0] + "T00:00", "%y%jT%H:%M")
    date = date + timedelta(hours= np.float( '0.' + date_tle.split('.')[1] ) * 24) 

    return date
