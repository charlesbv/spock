import os
import sys
filename_spock = sys.argv[2]
raise Exception
os.system("python concatenate_proc.py " + sys.argv[1] + " " + filename_spock)
os.system("python distance_ensemble_to_main_sc.py " + sys.argv[1] + " " + filename_spock + " 1")
os.system("python collision.py " + sys.argv[1] + " " + filename_spock)
