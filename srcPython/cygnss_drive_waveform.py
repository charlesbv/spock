# This script:
# - sets up all 8 FMs waveform generation and combination waveforms
# - runs them
# Notes:
# - For each FM, the script:
#    1- writes the two PRN waveform generation input files
#    2- successively runs beacon-sig-gen to generate each of the 2 PRN waveform binary files
#    3- runs beacon-sig-combine to combine the 2 PRN waveform binary files into one binary file
#   Each step is performed in the proper directory:
#    - steps 1- and 2- are performed in simu/YYYYMMDD/other/waveform_gen/fm0X/prn_Y, where YYYYMMDD is taken from start_time, X is 1, 2, ..., or 8, and Y is either 'i' (first PRN) or 'ii' (second PRN)
#    - step 3- is performed in simu/YYYYMMDD/other/waveform_combiner
# - In step 2-, the names of each of the 2 PRN output file follows the format cygnss_YYYYMMDD_passM_fmX_prn_ZZ.bin, where YYYYMMDD is taken from start_time, M is a number between 1 and 8, X is 1, 2, ..., or 8, and ZZ is 01, 02, ..., or 32
# In step 3-, the name of the combined output file follows the format cygnss_YYYYMMDD_passM_fmX_prns_ZZ_UU.bin, where YYYYMMDD is taken from start_time, M is a number between 1 and 8, X is 1, 2, ..., or 8, and ZZ and UU are 01, 02, ..., or 32

# Inputs:
# - start_time (%Y-%m-%dT%H:%M:%S) 
# - csv_filename_list (1st component is fm, 2nd component is prn csv filename -> 8 * 2 components)
# - prn_list (1st component is fm, 2nd component is prn -> 8 * 2 components). Important: format must be ZZ (e.g.: 03, not 3)
# - sp3_filename
# - fm_pass (1 to 8. For example, if FM07 is the third satellite to overpass the ebacon station then fm_pass[6] = 3)
import sys
import os
sys.path.append('/Users/cbv/work/spock/srcPython')
from cygnss_beacon_write_input_satbop import *
from cygnss_beacon_package import *
from pathlib import Path
from shutil import copyfile, move
from distutils.dir_util import copy_tree
import ipdb
from cygnss_beacon_write_input_waveform import *

# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT
start_time = '2018-10-31T00:00:00'
csv_filename_list = [['PRN_20_truncated.csv', 'PRN_21_truncated.csv'], # FM01
                     ['', ''], # FM02
                     ['', ''], # FM03
                     ['', ''], # FM04
                     ['', ''], # FM05
                     ['', ''], # FM06
                     ['', ''], # FM07
                     ['', '']] # FM08
prn_list = [['20', '21'],  # FM01
            ['', ''],  # FM02
            ['', ''],  # FM03
            ['', ''],  # FM04
            ['', ''],  # FM05
            ['', ''],  # FM06
            ['', ''],  # FM07
            ['', '']]  # FM08
sp3_filename = 'igs20254.sp3'
fm_pass = ['4', # FM01
           '', # FM02
           '', # FM03
           '', # FM04
           '', # FM05
           '', # FM06
           '', # FM07
           ''] # FM08

# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT


current_dir = str(Path().absolute()) + '/'
#os.chdir(current_dir) if script crases afet chanding dir, can go back to initial dir
start_time_dir = start_time[0:10].replace('-', "") + '/'

binary_output_filename = []
combined_binary = []
# For each FM
for cygfm in range(1,2):#!!!!! should be range(1,9):
    if fm_pass[cygfm-1] != '':
        cygfm_str = str(cygfm)
        ## For each of the 2 PRNs
        binary_output_filename_prn = []
        for iprn in range(2):
            ### Write the two PRN waveform generation input files
            if iprn == 0:
                prn_dir = 'prn_i/'
            else:
                prn_dir = 'prn_ii/'
            waveform_gen_dir = 'simu/' + start_time_dir + 'other/waveform_gen/fm0' + cygfm_str + '/' + prn_dir
            os.chdir(waveform_gen_dir)
            copy_tree(current_dir + 'simu/waveform_gen_original_adapted', '.') # sat-bop exe and files for the exe
            csv_filename = csv_filename_list[cygfm-1][iprn]
            waveform_gen_output_dir = cygnss_beacon_write_input_waveform(csv_filename, sp3_filename)
            ### Run beacon-sig-gen to generate each of the 2 PRN waveform binary files
            #os.system('./beacon-sig-gen.exe')
            ### Rename the output binary file in the format cygnss_YYYYMMDD_passM_fmX_prn_ZZ.bin
            for file in os.listdir(waveform_gen_output_dir + '/'): # only one file should have the .bin extension. That's the ouptut file that just got created by beacon-sig-gen, and that we want to rename
                if file.endswith(".bin"):
                    old_binary_output_filename = file
            new_binary_output_filename = 'cygnss_' + start_time[0:10].replace('-','') + '_pass' + fm_pass[cygfm-1] + '_fm' + cygfm_str + '_prn_' + prn_list[cygfm-1][iprn] + '.bin'
            move(waveform_gen_output_dir + '/' + old_binary_output_filename, waveform_gen_output_dir + '/' + new_binary_output_filename)
            binary_output_filename_prn.append(waveform_gen_output_dir + '/' + new_binary_output_filename)
            os.chdir(current_dir)
        binary_output_filename.append(binary_output_filename_prn)
        ## Run beacon-sig-combine to combine the 2 PRN waveform binary files into one binary file
        waveform_combiner_dir = 'simu/' + start_time_dir + 'other/waveform_combiner/'
        os.chdir(waveform_combiner_dir)
        ### Andrew wants a shell script that combines the two PRN waveforms so we make this script and run it
        combine_script_filename = 'combine_fm0' + cygfm_str + '.sh'
        combine_script_file = open(combine_script_filename, 'w')
        combined_binary_here = 'cygnss_' + start_time[0:10].replace('-','') + '_pass' + fm_pass[cygfm-1] + '_fm' + cygfm_str + '_prns_' + prn_list[cygfm-1][0] + '_' + prn_list[cygfm-1][1] +'.bin' #cygnss_YYYYMMDD_passM_fmX_prns_ZZ_UU.bin
        combined_binary.append( combined_binary_here )
        print >> combine_script_file, 'beacon-sig-combine.exe ' + combined_binary_here + ' ' + '../waveform_gen/fm0' + cygfm_str + '/prn_i/' + binary_output_filename_prn[0] +\
             ' ../waveform_gen/fm0' + cygfm_str + '/prn_ii/' + binary_output_filename_prn[1] 
        combine_script_file.close()
        os.chmod(combine_script_filename, 0o777)
        raise Exception
        os.system('./' + combined_binary_here)
