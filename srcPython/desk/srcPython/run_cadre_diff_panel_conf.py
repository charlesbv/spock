import os

#panel_conf = ['0panel', '1panel', '2panels', '3panels', '4panels']
panel_conf = ['1panel']


for ipanel in range(len(panel_conf)):
    os.chdir("../../")
    os.system("/usr/local/bin/mpirun -np 10 run_moat " + 'cadre_tumbling_' + panel_conf[ipanel] + '.txt')
    os.chdir("./code/python_propagator/")
    os.system("python concatenate_proc.py " + "cadre_tumbling_" + panel_conf[ipanel] + ".txt")
    os.system("python distance_ensemble_to_main_sc.py " + "cadre_tumbling_" + panel_conf[ipanel] + ".txt 1")
