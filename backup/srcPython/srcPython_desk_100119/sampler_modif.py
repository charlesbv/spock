import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
from scipy.stats import norm
sns.set_style('white')
sns.set_context('talk')

def sampler_modif(data, total_sample_target = 10000, mu_init=.5, proposal_width=.5, plot=False, accept_criteria  = 500):
    mu_current = mu_init
    posterior = [mu_current]
    nb_accept = 0
    nb_tries = 0
    while ( len(posterior) < total_sample_target ):
        # suggest new position
        mu_proposal = norm(mu_current, proposal_width).rvs()

        if ( ( accept_criteria- mu_proposal >=0) & (mu_proposal >=0) ): 
            accept = 1
        else:
            accept = 0
#        print len(posterior),mu_proposal
        if accept:
            # Update position
            mu_current = mu_proposal
            nb_accept = nb_accept + 1
            posterior.append(mu_current)
        nb_tries = nb_tries + 1
    print 'acceptance rate is:', nb_accept / np.double(nb_tries) * 100
        
    return posterior


