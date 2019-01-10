from datetime import datetime, timedelta
import numpy as np

def seconds_since_epoch_start_to_utc(epoch_start_in_utc, seconds_since_epoch_start): # seconds_since_epoch_start can be a decimal number (precision: microseconds). IMPORTANT: epoch_start_in_utc has to be datetime object 

    utc_date = epoch_start_in_utc + timedelta(seconds = seconds_since_epoch_start)

    return utc_date
