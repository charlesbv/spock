import os
os.system("qstat -u cbussy > log.qsat")
file = open('log.qsat')
read_file = file.readlines()
n = len(read_file) - 3

for ifile in range(n):
    job_id = read_file[3 + ifile].split('.')[0]
    os.system("qdel " + job_id)
