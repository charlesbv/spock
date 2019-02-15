from datetime import datetime
from dateutil import tz
from pytz import timezone

# This script converts a time from UTC to Eastern time (EST), taking into account daylight saving time
# input: utc time as a string '%Y-%m-%dT%H:%M:%S'
def utc2est(utc_str):

    # METHOD 1: Hardcode zones:
    from_zone = tz.gettz('UTC')
    utc = datetime.strptime(utc_str, '%Y-%m-%dT%H:%M:%S')
    utc_ok = utc.replace(tzinfo=from_zone)
    est = utc_ok.astimezone(timezone('US/Eastern'))
    est_str = datetime.strftime(est, '%Y-%m-%dT%H:%M:%S')   
    return est_str
