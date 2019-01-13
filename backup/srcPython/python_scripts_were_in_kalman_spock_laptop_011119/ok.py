
r_spock_ref = var_out[find_in_read_input_order_variables(var_out_order, 'position')]
v_spock_ref = var_out[find_in_read_input_order_variables(var_out_order, 'velocity')]


# Compare SpOCK and data
# Assumption: SpOCK was run with a 1s time step to avoid having to do interpolation here: the steps in SpOCK falls at the same time as the steps in data 
## Select the time where date_spock = date_obs 


print 'Computing distances between ensembles and observations'
if it == 0:
    index_spock_same_date_as_obs = []
    index_obs_kept = []
    iobs = 0
    while iobs < nb_obs:
        if date_obs[iobs] >= date_datetime_round_sec_spock[-1]:
            break
        else:
            index_spock_same_date_as_obs.append(np.where(date_datetime_round_sec_spock == date_obs[iobs])[0][0])
            index_obs_kept.append(iobs)
            nb_seconds_since_start.append( ( date_obs[iobs] - date_obs[0] ).total_seconds() )
        iobs = iobs + 60 # jump by 60 observation time steps (uusualyy one time step is one second)    #!!!!!!! used to be iobs + 1

date_obs_arr = np.array(date_obs)
date_obs_ok = date_obs_arr[index_obs_kept]
r_obs_ok = r_obs[index_obs_kept]
v_obs_ok = v_obs[index_obs_kept]

n = len(index_obs_kept)# !!!!! used to be iobs #!!!!!!!!!! j-index_interval[iinter]

date_spock_ok = date_spock[index_spock_same_date_as_obs]
nb_seconds_since_start_spock_ok = nb_seconds_since_start_spock[index_spock_same_date_as_obs]

r_spock_ref_ok = np.zeros([n, 3])
r_spock_ref_ok[:, 0] = r_spock_ref[index_spock_same_date_as_obs, 0]
r_spock_ref_ok[:, 1] = r_spock_ref[index_spock_same_date_as_obs, 1]
r_spock_ref_ok[:, 2] = r_spock_ref[index_spock_same_date_as_obs, 2]
v_spock_ref_ok = np.zeros([n, 3])
v_spock_ref_ok[:, 0] = v_spock_ref[index_spock_same_date_as_obs, 0]
v_spock_ref_ok[:, 1] = v_spock_ref[index_spock_same_date_as_obs, 1]
v_spock_ref_ok[:, 2] = v_spock_ref[index_spock_same_date_as_obs, 2]

