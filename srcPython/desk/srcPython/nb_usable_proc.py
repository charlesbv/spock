# This script finds the number of processors you can use when runninng mpirun. Other solutions (such as multiprocessing.cpu_count()) return more processors than actually available from mpirun.

import os
def nb_usable_proc(path_mpirun):
    os.system(path_mpirun + " hostname > nbcore_temp")
    file_nb_core = open("nbcore_temp")
    read_file_nb_core = file_nb_core.readlines()
    nb_core = len(read_file_nb_core)
    os.system("rm -f nbcore_temp")
    return nb_core

