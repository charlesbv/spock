import os 

run_names = [filename for filename in os.listdir('../../input/main_input') if ('SCION' in filename and filename.endswith('.txt'))]
n = len(run_names)
for i in range(n):
    print run_names[0:i+1]
    os.chdir("../../")
    os.system("time /usr/local/bin/mpirun -np 10 run_moat " + run_names[i])
    os.chdir("code/python_propagator")
    os.system("python concatenate_proc.py " + run_names[i])
    os.system("python distance_ensemble_to_main_sc.py " + run_names[i] + " 1")








