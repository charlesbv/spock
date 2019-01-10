import numpy as np
import numpy.linalg as li

# same as eci_to_lvlh but replaced eci with ecef (i could directly input ecef into eci_to_lvlh but to make it consistent with the name i copy pasted eci_to_lvlh here and chaned eci to ecef)
def ecef_to_lvlh(position_ecef, velocity_ecef, vector_ecef):
    "This function transforms a vector from the interial ECEF coordinates to the LVLH coordinates "
    # Radial (1/2)
    radial_minus = -position_ecef / li.norm(position_ecef)

    # Cross-track
    cross_track_temp = np.cross(radial_minus, velocity_ecef)
    cross_track = -cross_track_temp / li.norm(cross_track_temp)

    # Along-track
    along_track_temp = np.cross(cross_track, radial_minus)
    along_track = - along_track_temp / li.norm(along_track_temp)

    # Radial (2/2)
    radial = -radial_minus

    # Matrix conversion
    T_ecef_2_lvlh = np.zeros([3, 3])
    T_ecef_2_lvlh[0, :] = along_track
    T_ecef_2_lvlh[1, :] = cross_track 
    T_ecef_2_lvlh[2, :] = radial

    # Matrix * vector_ecef
    vector_lvlh = np.dot(T_ecef_2_lvlh , vector_ecef)

    return vector_lvlh;
