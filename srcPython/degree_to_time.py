# THIS SCRIPT 
# Set twelve_hour_offset to 0 if you want to start at 12 am: 0 degree is midnight, 90 degrees is 6 am, 180 degrees is 12 pm, 270 degrees is 6 pm  
# Set twelve_hour_offset to 1 if you want to start at 12 pm and not 12 am: 0 degree is 12 pm, 90 degrees is 6 pm, 180 degrees is midnight, 270 degrees is 6 am

from datetime import datetime, timedelta

def degree_to_time(var_in_deg, twelve_hour_offset):
    # var_in_deg = [0, 90, 180, 270, 225]
    # twelve_hour_offset = 1
    n = len(var_in_deg)
    var_in_time = []
    if twelve_hour_offset == 0: # 0 degree is midnight, 90 degrees is 6 am, 180 degrees is 12 pm, 270 degrees is 6 pm
        time_start = datetime.strptime("00:00:00", "%H:%M:%S")
    if twelve_hour_offset == 1: # 0 degree is 12 pm, 90 degrees is 6 pm, 180 degrees is midnight, 270 degrees is 6 am
        time_start = datetime.strptime("12:00:00", "%H:%M:%S")
    for istep in range(n):
        var_in_hour = var_in_deg[istep] / 360. * 24
        var_in_hour_date_format = time_start + timedelta(hours=var_in_hour) 
        var_in_time.append( datetime.strftime(var_in_hour_date_format, "%H:%M:%S") )

    return var_in_time
