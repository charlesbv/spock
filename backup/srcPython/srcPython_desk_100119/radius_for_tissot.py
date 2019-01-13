import numpy as np
def radius_for_tissot(dist_km):
    return np.rad2deg(dist_km/6378.) # this is calculating using the Haversine formula
