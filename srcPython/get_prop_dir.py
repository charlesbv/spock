# This function is pretty useless but just made if other users use this propagator and the python scripts I wrote

import os 

def get_prop_dir(dir_offset): # dir_offset is how deep the current directory is from the propagator directory
    propagator_directory_list = os.getcwd().split('/')[0:-dir_offset]
    propagator_directory = ""
    for i in range(len(propagator_directory_list)):
        if (i > 0):
            propagator_directory = propagator_directory + "/" + propagator_directory_list[i]
            
    propagator_directory = propagator_directory + '/'
    return propagator_directory
