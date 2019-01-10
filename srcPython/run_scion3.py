import os 

#run_names = ['SCION_2_5_deg_s_f107_120_ap_15.txt','SCION_4_5_deg_s_f107_120_ap_15.txt','SCION_3_5_deg_s_f107_120_ap_15.txt','SCION_0_5_deg_s_f107_120_ap_15.txt','SCION_1_5_deg_s_f107_120_ap_15.txt']

run_names = ['SCION_2_deg_s_f107_90_ap_7_90_days_10_s_500_ensembles.txt', 'SCION_2_deg_s_f107_120_ap_15_90_days_10_s_500_ensembles.txt', 'SCION_2_deg_s_f107_200_ap_80_90_days_10_s_500_ensembles.txt']

n = len(run_names)
for i in range(n):
    print run_names[0:i+1]
    os.chdir("../run_scion/")
    os.system("time /usr/local/bin/mpirun -np 10 run_moat " + run_names[i])
    os.chdir("../srcPython/")
    os.system("python concatenate_proc.py run_scion " + run_names[i])
    os.system("python distance_ensemble_to_main_sc.py run_scion " + run_names[i] + " 1")








