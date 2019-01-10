import numpy as np

def orbit_average_with_value_at_equator(var_to_average, latitude, time_before_average):
    # isat = 0
    # var_to_average = power[isat, :, 0]
    # time_before_average = date
    #latitude = latitude[isat, :] 
    n = len(var_to_average)
    var_orbit_averaged = []
    var_at_equator = []
    time_averaged = []
    index_time_averaged = []
    istep = 0
    passed_first_ascending_node = 0
    while istep < n:
        var_orbit_averaged_temp = []
        time_averaged_temp = []
        index_time_averaged_temp = []
        while ( ( latitude[istep ] > 0 ) | ( latitude[istep  + 1] <= 0 ) ):
            # start average at the first ascending node (ignore the steps before crossing the ascending node for the first time)
            if ( passed_first_ascending_node == 1 ):
                var_orbit_averaged_temp.append( var_to_average[ istep ] )
                if len(var_orbit_averaged_temp) == 1:
                    time_averaged_temp.append( time_before_average[istep] )
                    index_time_averaged_temp.append( istep )
                    # linear interpolate to get the value of var at the equator
                    x0 = latitude[istep-1]
                    x1 = latitude[istep]
                    y0 = var_to_average[ istep -1]
                    y1 = var_to_average[ istep ]
                    slope = (y1 - y0) / (x1 - x0)
                    b = y0 - slope * x0
                    y = slope * 0 + b # latitue = 0 at the equator
                    #print x0,y0,y,x1,y1
                    #var_at_equator.append(y)
                    var_at_equator.append(var_to_average[ istep ])
                    istep_beginning = istep
            istep = istep + 1
            if istep + 1 >= n:
                break
        if ( passed_first_ascending_node == 1 ):
            var_orbit_averaged_temp.append( var_to_average[ istep ] )
            time_averaged_temp.append( time_before_average[(istep_beginning + istep) / 2] )
            time_averaged_temp.append( time_before_average[istep] )
            index_time_averaged_temp.append( (istep_beginning + istep) / 2 )
            index_time_averaged_temp.append( istep )
            var_orbit_averaged.append( np.mean( var_orbit_averaged_temp ) )
            time_averaged.append( time_averaged_temp )
            index_time_averaged.append( index_time_averaged_temp )
        passed_first_ascending_node = 1
        istep = istep + 1
        if istep + 1 >= n:
            break

    # Delete the last element because it likely corresponds to an average over only a part of the orbit
    del time_averaged[-1]
    del var_orbit_averaged[-1]
    del index_time_averaged[-1]
    del var_at_equator[-1]

    return var_orbit_averaged, time_averaged, index_time_averaged, var_at_equator  # time_averaged[0] start of the orbit over which the average was made
                                                                 # time_averaged[1] middle of the orbit over which the average was made
                                                                 # time_averaged[2] end of the orbit over which the average was made
                                                                 # index_time_averaged index in time_before_average of the orbit over which the average was made
