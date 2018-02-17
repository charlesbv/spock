from get_prop_dir import *
from os import listdir
from os.path import isfile, join

dir_to_open = get_prop_dir(1) + "run_aerie/input/density/density_NRLMSIS00e/"
file_to_modify_arr = [f for f in listdir(dir_to_open) if isfile(join(dir_to_open, f))]

for ifile in range(len(file_to_modify_arr)):
    file_to_modify = open(dir_to_open + file_to_modify_arr[ifile])
    file_out = open(dir_to_open + 'modif_' + file_to_modify_arr[ifile], "w+")
    read_file_to_modify = file_to_modify.readlines()
    n_header = 3
    n = len(read_file_to_modify) - n_header - 1 #  -1 to remove last line '#ENDOFFILE'
    print >> file_out, "#BEGINNINGOFHEADER" 
    print >> file_out, "#ENDOFHEADER" 
    print >> file_out, "YEAR DOY HR    1" 
    for istep in range(n):
        print >> file_out, read_file_to_modify[istep + n_header].split()[0] + " " + read_file_to_modify[istep + n_header].split()[1] + " 0 " + read_file_to_modify[istep + n_header].split()[2]
    print >> file_out, "#ENDOFFILE" 
    file_out.close()
    file_to_modify.close()
