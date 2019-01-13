# This script creates a job file that runs SpOCK with the main input file chosen as argument of this script. The name of the job file is "job." + spock_main_input_filename (the job file name is without the extension of SpOCK main input file). The second argument of the script is the queue to use.
# i couldnt find spock_job_file_new.py which i used for CA2 so I copied spock_job_file.py and called it spock_job_file_new_july18.py (in jul 2018)


import os
def spock_job_file_new_july18(spock_main_input_file, queue, walltime_in, select_in, ncpus_in, folder_tar):
    if folder_tar[-1] != '/':
        folder_tar = folder_tar + '/'

    # parameters of the run
    walltime=walltime_in
    select=select_in
    ncpus=ncpus_in
    # end of parameters of the run

    if ( (queue != "low") & (queue != "devel") & (queue != "debug")): # in case the second argument is messed up in the call of this script
        queue = "low"


    filename_job = "job." + spock_main_input_file.replace(".txt", "")
    file_job = open(filename_job, "w+")
    
    print >> file_job, "#!/bin/csh\n#PBS -W group_list=s1815\n#PBS -S /bin/csh\n#PBS -N SpOCK." + spock_main_input_file.replace(".txt", "") + "\n#PBS -l walltime=" + walltime + "\n#PBS -l select=" + str(select) + ":ncpus=" + str(ncpus) + ":model=bro\n#PBS -j oe\n#PBS -m e\n#PBS -q " + queue + "\n\nmodule load mpi-sgi/mpt\ncd $PBS_O_WORKDIR"
    print >> file_job, "setenv MPI_UNBUFFERED_STDIO 1"
    print >> file_job, 'set run_spock="' + spock_main_input_file.replace(".txt", "") + '"'
    print >> file_job, "mpiexec /home1/cbussy/spock_eq6 $run_spock.txt > log.$run_spock"
    print >> file_job, 'tar -zcvf ' + folder_tar + '"$run_spock"_out.tar.gz $run_spock\ncp $run_spock.txt ' + folder_tar + '"$run_spock"_in.txt\n# rm -Rf $run_spock\nchmod 777 ' + folder_tar + '"$run_spock"_out.tar.gz\nchmod 777 ' + folder_tar + '"$run_spock"_in.txt'
    
    print >> file_job, 'scp -p ' + folder_tar + '"$run_spock"_out.tar.gz '  + folder_tar + '"$run_spock"_in.txt lfe:./ca3/' + folder_tar  #!! eneeed to make sure that folder_tar exissts on lfe in ./ca3
    file_job.close()
    return filename_job
