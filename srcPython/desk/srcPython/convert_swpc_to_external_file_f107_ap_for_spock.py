# This script converts a file from SWPC into the format for SpOCK for external files (Ap and F10.7). The script figures out itself if the input file is F10.7 or Ap (as long as it's a file downloaded from SWPC)
# ASSUMPTIONS

filename_to_convert = "2017Q3_DGD.txt"
file_to_convert = open(filename_to_convert)
read_file_to_convert = file_to_convert.readlines()
file_to_convert.close()

# Figure out if the file is Ap or F10.7: 1st line has 'Geomagnetic' if Ap, 'Solar' if F10.7
if 'Geomagnetic' in read_file_to_convert[0]:
    param = 'ap'
else:
    param = 'f107'

# skip header (first line not header starts with a '2' (for the year))
nb_header = 0
while (read_file_to_convert[nb_header][0] != '2'):
    nb_header = nb_header + 1
nb_header = nb_header + 1

