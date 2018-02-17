import os
import sys

# list_run = ['devel_1126_ok_quartile_f107_1_quartile_ap_1.txt',
#              '1126_ok_quartile_f107_2_quartile_ap_2.txt',
#              '1126_ok_quartile_f107_3_quartile_ap_3.txt',
#              '1126_ok_quartile_f107_4_quartile_ap_4.txt',
#              '1126_ok_quartile_f107_5_quartile_ap_5.txt',
#              '1126_ok_quartile_f107_6_quartile_ap_6.txt',
#              '1126_ok_quartile_f107_7_quartile_ap_7.txt',
#              '1126_ok_quartile_f107_8_quartile_ap_8.txt',
#              '1126_ok_quartile_f107_9_quartile_ap_9.txt']

list_run = ['with_tca_storm_dynamic_32600ens_min_dist_coll_1_2m.txt']#with_tca_storm_static_32600ens_min_dist_coll_1_2m.txt']
n = len(list_run)
for i in range(n):
    print i, n-1
    # os.chdir("/raid4/cbv/installation_spock/spock/run_collision/output")
    # os.system("tar -zxvf " + list_run[i].replace(".txt","_out.tar.gz"))

    os.system("/usr/local/bin/mpirun -np 12 python mpi_concatenate_proc.py run_collision " + list_run[i])
    # also mpi_distance_ensemble_to_main_sc
#    os.system("/usr/local/bin/mpirun -np 12 python mpi_distance_ensemble_to_main_sc.py run_collision " + list_run[i] + " x_eci y_eci z_eci f107 ap rho latitude true_anomaly argument_perigee")

    os.system("/usr/local/bin/mpirun -np 12 python mpi_distance_ensemble_to_main_sc.py run_collision " + list_run[i] + " tca")
