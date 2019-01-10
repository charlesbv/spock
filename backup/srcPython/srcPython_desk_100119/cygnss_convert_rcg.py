# This script converts the CYGNSS RCG to/from the 0-15 scale from/to the 0/220 scale, based on the main.cpp script (function writeDdmiFile) that Tim Butler sent me by email on 09/21/2018)
# Inputs:
# - gain_in: numpy array of gain (n*1) from orginal scale (see second input)
# - scale: name of orignial scale. 'small' if orginail scale is 0-15, 'large' is orignial scale is 0-220
# Output:
# - gain_out: array of gain from final scale (see input "scale")

import numpy as np
import ipdb
def cygnss_convert_rcg(gain_in, scale):
    if  ( ( scale != "small" ) & ( scale != "large" ) ):
        print "***! The scale needs to be either 'small' or 'large'. The program will stop.!***"; raise Exception
    
    scale_factor = 16./np.log(220) - 0.001; #     double scale_factor = (16 / log(maxRcg_)) - 0.001;
    n = len(gain_in)
    gain_out = np.zeros([n])
    if scale == "small":
        where_0 = np.where(gain_in == 0)[0]
        where_gt_0 = np.where(gain_in > 0)[0]

        gain_out[where_0] = 0.5 # arbitraty, could be any value between 0 and 1. 
                                # this is because if the high scale (0-220) gain is lower than 1, it's scled to 0
        gain_out[where_gt_0] = np.exp(gain_in[where_gt_0] / scale_factor)

    if scale == "large":
        ipdb.set_trace()
        where_ge_1 = np.where(gain_in >= 1)[0]
        where_lt_1 = np.where(gain_in < 1)[0]
        gain_out[where_lt_1] = 0
        gain_out[where_ge_1] = np.floor(np.log(gain_in[where_ge_1]) * scale_factor);        
#             int scaled =
#                 (orig >= 1) ? floor(log(orig) * scale_factor) : 0;


    return gain_out
