# This script creates a job file that runs SpOCK with the main input file chosen as argument of this script. The name of the job file is "job." + spock_main_input_filename (the job file name is without the extension of SpOCK main input file). The second argument of the script is the queue to use.


import os
def spock_job_file(spock_main_input_file, queue, walltime_in, select_in, ncpus_in):

    # parameters of the run
    walltime=walltime_in
    select=select_in
    ncpus=ncpus_in
    # end of parameters of the run

    if ( (queue != "low") & (queue != "devel") & (queue != "debug")): # in case the second argument is messed up in the call of this script
        queue = "low"

    os.chdir("../run_collision")
    filename_job = "job." + spock_main_input_file.replace(".txt", "")
    file_job = open(filename_job, "w+")
    
    print >> file_job, "#!/bin/csh\n#PBS -W group_list=g26192\n#PBS -S /bin/csh\n#PBS -N SpOCK." + spock_main_input_file.replace(".txt", "") + "\n#PBS -l walltime=" + walltime + "\n#PBS -l select=" + str(select) + ":ncpus=" + str(ncpus) + ":model=bro\n#PBS -j oe\n#PBS -m e\n#PBS -q " + queue + "\n\ncd $PBS_O_WORKDIR"
    print >> file_job, "setenv MPI_UNBUFFERED_STDIO 1"
    print >> file_job, 'set run_spock="' + spock_main_input_file.replace(".txt", "") + '"'
    print >> file_job, "/usr/bin/time -p mpiexec spock $run_spock.txt > log.$run_spock"
    print >> file_job, 'cd output\ntar -zcvf ./"$run_spock"_out.tar.gz ./$run_spock\ncp ../input/main_input/$run_spock.txt ./"$run_spock"_in.txt\nrm -Rf $run_spock\nchmod 777 ./"$run_spock"_out.tar.gz\nchmod 777 ./"$run_spock"_in.txt\ncd ..'

    file_job.close()
    os.chdir("../srcPython")
    return filename_job
