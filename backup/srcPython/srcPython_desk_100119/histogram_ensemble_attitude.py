import fileinput
import numpy as np
from matplotlib import pyplot as P
#import pylab as P
import os


# file_mean_satellite = open('../output/power_CYGNSS_1.txt')
# aa = file_mean_satellite.readlines()
# total_power = 0L
# len_file_mean_satellite = len(aa) - 13
# power_mean_satellite = np.zeros([len_file_mean_satellite,5])
# total_power_mean_satellite = np.zeros(len_file_mean_satellite)
# for i in range(len_file_mean_satellite):
#     power_mean_satellite[i,0] = float(aa[13+i].split()[2]) # "/ 3.0" because STK outputs the power for each of the 3 top solar panels. The propagator outputs the SUM of the 3 top solar panels' power.
#     power_mean_satellite[i,1] = float(aa[13+i].split()[2+4])
#     power_mean_satellite[i,2] = float(aa[13+i].split()[2+14])
#     total_power_mean_satellite[i] = power_mean_satellite[i,0] * 3.0 + power_mean_satellite[i,1] + power_mean_satellite[i,2] # "* 3.0" because STK outputs the power for each of the 3 top solar panels. The propagator outputs the SUM of the 3 top solar panels' power.

# P.plot(power_mean_satellite[:,0], color = 'm', linewidth = 2)
# P.plot(power_mean_satellite[:,1], color = 'b', linewidth = 2)
# P.plot(power_mean_satellite[:,2], color = 'g', linewidth = 2)
# P.plot(total_power_mean_satellite, color = 'r', linewidth = 2)
# P.ylabel('Power (W)')
# P.xlabel('Time (minutes)')
# P.title('Nadir pointing - fully deployed - S/A eff = 28% - our propagator')
# P.yticks(np.arange(min(total_power_mean_satellite), max(total_power_mean_satellite), 5))
# P.show()
# mean_power_mean_satellite = total_power / len_file_mean_satellite
# # #ENDOFHEADER

# raise Exception

ensemble_to_plot = ['x_eci', 'y_eci', 'z_eci']
nb_histo =  1
time_histo_array = [645]

# read input file
filename = '../input.d'
file = open(filename,'r')
a = file.readlines()

nb_hours = 5#int(a[1][3:8])
nb_seconds = nb_hours * 3600L
dt = int(float(a[3][0:10]))
nb_time_steps = int(np.floor(nb_seconds / dt)) + 1
nb_spacecraft = int(a[9][0:6])
satellite_to_plot= []
for i in range(nb_spacecraft):
    satellite_to_plot.append( a[11][0:30*nb_spacecraft].split(',')[i].strip() )

nb_ensembles_coe = int(a[28][0:6])
nb_ensembles_attitude = int(a[42][0:6])
# calculate nb_ensembles_cd
name_file_surface = a[14].split()[0]
file_surface = open("../moat/input/"+name_file_surface)
a_file_surface = file_surface.readlines()
read_next_line = 0
for line_geo in fileinput.input("../moat/input/"+name_file_surface, inplace = False):
    if (read_next_line == 1):
        nb_ensembles_cd = (int)(line_geo[0:10])
        break
    if ( line_geo[0:14] == "# NB_ENSEMBLES" ):
        read_next_line = 1
nb_ensembles_array = [nb_ensembles_coe, nb_ensembles_attitude, nb_ensembles_cd]
nb_ensembles = np.max(nb_ensembles_array)
for i in range(len(nb_ensembles_array)):
    if ( (nb_ensembles_array[i] > 0 ) &  (nb_ensembles_array[i] < nb_ensembles) ) :
        nb_ensembles = nb_ensembles_array[i]

## number of processors
nb_ensemble_to_plot = len(ensemble_to_plot)
nProc = 0
for file in os.listdir("../output/ensemble/proc_files/"):
    nProc = nProc + 1
    
nProc = nProc / ( nb_ensemble_to_plot * nb_spacecraft)
if ( nProc == 0 ):
    print "No iproc files found."
    quit()
nb_ensemble_per_proc = int(np.floor(nb_ensembles / nProc))
nb_ensembles_ok = nb_ensemble_per_proc * nProc

# read ensemble files
for sss in range(nb_spacecraft):
    P.figure(figsize=(23, 12))
    for eee in range(nb_ensemble_to_plot):
            
        location_subplot = eee + 1
#        P.subplot(2,3,location_subplot)

        filename = 'ensemble_'+ensemble_to_plot[eee]+'_'+satellite_to_plot[sss]
        file = open('../output/ensemble/'+filename,'r')
        a = file.readlines()
        line = ["" for x in range(nb_time_steps)]
        time = ["" for x in range(nb_time_steps)]
        x_alltime_ensemble = np.zeros([nb_time_steps, nb_ensembles_ok-1])

        for line_count in range(nb_time_steps):
            line[line_count] = (a[13+line_count]).split(' ')
            for ensemble in range(nb_ensembles_ok-1):
                x_alltime_ensemble[line_count][ensemble] = float( line[line_count][ensemble+2]  )
        file.close()

    x_mean_ensemble_over_time = np.zeros([nb_ensembles_ok-1])
    for mmm in range(nb_ensembles_ok-1):
        x_mean_ensemble_over_time[mmm] = np.mean(x_alltime_ensemble[:, mmm])
    
    name_satellite_to_plot = satellite_to_plot[sss][:-4]
  #  P.title(name_satellite_to_plot+' - '+ensemble_to_plot[eee])


# mmm= 458
# P.plot(x_alltime_ensemble[:,mmm],linewidth=2, color = 'b')
# P.plot([0,nb_time_steps-1],[x_mean_ensemble_over_time[mmm],x_mean_ensemble_over_time[mmm]], linewidth=2, color = 'r')
# P.ylim([0,max(x_alltime_ensemble[:,mmm])])
# P.yticks(np.arange(0, 70, 5))
# P.ylabel('Power (W)')
# P.xlabel('Time (seconds)')
# P.title('Our propagator')
# print x_mean_ensemble_over_time[mmm]
# P.show()

# raise Exception


# ################ PLOT AVERAGE OVER TIME
    x_histo = x_mean_ensemble_over_time

    n, bins, patches = P.hist(x_histo, 100,  histtype='stepfilled', alpha = 0.7) # 
    P.xlim([min(bins), max(bins)])
    P.xlabel("Power (W) - bins", size = 20)
    P.ylabel("Distribution: # per bin", size = 20)
    P.title("Average power distribution - S/A eff = 28% - sigma_omega = 4.8 deg/s", size = 24)
    P.xticks(np.arange(10, 28, 2))
    P.show()
# ################ PLOT AT A GIVEN TIME
# #         for time_histo_temp in range(nb_histo):
# # #            if (nb_histo == 1):
# #             time_histo = time_histo_array[time_histo_temp]#np.floor( ( ( nb_seconds ) * time_histo_temp  ) / dt )
# #  #           else:
# #   #              time_histo = time_histo_array[time_histo_temp]#np.floor( ( ( nb_seconds / ( nb_histo - 1) ) * time_histo_temp  ) / dt )
# #             # if ( time_histo_temp == nb_histo - 1 ):
# #             #     time_histo = time_histo - 1
# #             x_histo = np.zeros([nb_ensembles_ok])
# #             x_histo = x_alltime_ensemble[time_histo,:]

# #             if (nb_histo == 1 ):
# #                 n, bins, patches = P.hist(x_histo, 100, normed=1, histtype='stepfilled', label=[str((int)(time_histo * dt / 60)) + ' minutes'], alpha = 0.7)
# #             else:
# #                 n, bins, patches = P.hist(x_histo, 100, normed=1, histtype='stepfilled',label=[str((int)(time_histo * dt / 60)) + ' minutes'], alpha = 0.7)
# # #P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

# #             P.legend(loc = 'upper right')

#     P.tight_layout()
# #P.show()



# file_mean_satellite = open('../output/power_CYGNSS_1.txt')
# aa = file_mean_satellite.readlines()
# total_power = 0L
# len_file_mean_satellite = len(aa) - 13
# power_mean_satellite = np.zeros(len_file_mean_satellite)
# for i in range(len_file_mean_satellite):
#     for j in range(13):
#         total_power = total_power + float(aa[13+i].split()[2 + j])
#         power_mean_satellite[i] = power_mean_satellite[i] + float(aa[13+i].split()[2 + j])
        
# mean_power_mean_satellite = total_power / len_file_mean_satellite
