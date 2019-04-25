# This script creates the structure of the package to be given to the team following the specification "Beacon Signal Package Spec" written by Andrew on 04/16/2019
# Inputs:
# - start_time (%Y-%m-%dT%H:%M:%S) 
# - end_time

import sys
import os
sys.path.append('/Users/cbv/work/spock/srcPython')



def cygnss_beacon_package(start_time, end_time):
    # main directory
    main_package_dir = 'package/' + start_time[:10].replace('-', '') + '/'
    if (os.path.isdir(main_package_dir) == False):
        os.system('mkdir ' + main_package_dir)

    # main_package_dir/other directory
    other_package_dir = main_package_dir + 'other/'
    if (os.path.isdir(other_package_dir) == False):
        os.system('mkdir ' + other_package_dir)

    # main_package_dir/other/sat-bop directory
    satbop_package_dir = other_package_dir + 'sat-bop/'
    if (os.path.isdir(satbop_package_dir) == False):
        os.system('mkdir ' + satbop_package_dir)

    # main_package_dir/other/waveform_gen directory
    waveform_gen_package_dir = other_package_dir + 'waveform_gen/'
    if (os.path.isdir(waveform_gen_package_dir) == False):
        os.system('mkdir ' + waveform_gen_package_dir)

    # main_package_dir/other/waveform_combiner directory
    waveform_combiner_package_dir = other_package_dir + 'waveform_combiner/'
    if (os.path.isdir(waveform_combiner_package_dir) == False):
        os.system('mkdir ' + waveform_combiner_package_dir)

    # main_package_dir/other/waveform_verification directory
    waveform_verification_package_dir = other_package_dir + 'waveform_verification/'
    if (os.path.isdir(waveform_verification_package_dir) == False):
        os.system('mkdir ' + waveform_verification_package_dir)

    # campaign directory
    campaign_package_dir = 'package/campaign/'
    if (os.path.isdir(campaign_package_dir) == False):
        os.system('mkdir ' + campaign_package_dir)

    # campaign_package_dir/sat-bop directory
    campaign_satbop_package_dir = campaign_package_dir + 'sat-bop/'
    if (os.path.isdir(campaign_satbop_package_dir) == False):
        os.system('mkdir ' + campaign_satbop_package_dir)

    # campaign_package_dir/waveform_gen directory
    campaign_waveform_gen_package_dir = campaign_package_dir + 'waveform_gen/'
    if (os.path.isdir(campaign_waveform_gen_package_dir) == False):
        os.system('mkdir ' + campaign_waveform_gen_package_dir)

    # campaign_package_dir/waveform_combine directory
    campaign_waveform_combine_package_dir = campaign_package_dir + 'waveform_combine/'
    if (os.path.isdir(campaign_waveform_combine_package_dir) == False):
        os.system('mkdir ' + campaign_waveform_combine_package_dir)

    ###############
    # SAME STURCUTRE FOR THE RUN DIRECTORY
    # main directory
    main_simu_dir = 'simu/' + start_time[:10].replace('-', '') + '/'
    if (os.path.isdir(main_simu_dir) == False):
        os.system('mkdir ' + main_simu_dir)

    # main_simu_dir/other directory
    other_simu_dir = main_simu_dir + 'other/'
    if (os.path.isdir(other_simu_dir) == False):
        os.system('mkdir ' + other_simu_dir)

    # main_simu_dir/other/sat-bop directory
    satbop_simu_dir = other_simu_dir + 'sat-bop/'
    if (os.path.isdir(satbop_simu_dir) == False):
        os.system('mkdir ' + satbop_simu_dir)
    fm_dir = []
    for cygfm in range(1,9):
        fm_dir.append( satbop_simu_dir + 'fm0' + str(cygfm) + '/' )
        if (os.path.isdir(fm_dir[-1]) == False):
            os.system('mkdir ' + fm_dir[-1])
        
    # main_simu_dir/other/waveform_gen directory
    waveform_gen_simu_dir = other_simu_dir + 'waveform_gen/'
    if (os.path.isdir(waveform_gen_simu_dir) == False):
        os.system('mkdir ' + waveform_gen_simu_dir)

    # main_simu_dir/other/waveform_combiner directory
    waveform_combiner_simu_dir = other_simu_dir + 'waveform_combiner/'
    if (os.path.isdir(waveform_combiner_simu_dir) == False):
        os.system('mkdir ' + waveform_combiner_simu_dir)

    # main_simu_dir/other/waveform_verification directory
    waveform_verification_simu_dir = other_simu_dir + 'waveform_verification/'
    if (os.path.isdir(waveform_verification_simu_dir) == False):
        os.system('mkdir ' + waveform_verification_simu_dir)
        
    return satbop_package_dir, satbop_simu_dir, fm_dir
