import numpy as np

def compute_power(filename):
    #filename = '/home/cbv/PropSim/output/run_aerie_power_only_tail_and_normalized/aerie_power_only_tail_and_normalized/power_aerie_power_only_tail_and_normalized.txt'
    file = open(filename, "r")
    read_file = file.readlines()
    # SKIP HEADER
    n_header = 0
    while (read_file[n_header].split()[0] != '#START'):
        n_header = n_header + 1
    n_header = n_header + 1
    n = len(read_file) - n_header
    nb_surfaces = len(read_file[n_header].split()) - 3 # -3 to skip the time and the Sun elevation angle
    power = np.zeros([n, nb_surfaces])
    sun_elevation = np.zeros([n])
    for istep in range(1,n):
        for isurf in range(nb_surfaces):
            power[istep, isurf] = read_file[istep + n_header].split()[isurf + 2]
        sun_elevation[istep] = read_file[istep + n_header].split()[nb_surfaces+1 + 2]
    return power, sun_elevation
