# This script converts quaternions from the engineering convention to the SPICE convention.
# Engineering convnetion:
# qeng = ( sin(theta/2) ex, sin(theta/2) ey, sin(theta/2) ez, cos(theta/2))
# SPICE convention:
# qspice = ( cos(theta/2), sin(theta/2) ex, sin(theta/2) ey, sin(theta/2) ez)
# Careful: according to the SPICE script q2m_c.c, the engineering convention is  qeng = ( -sin(theta/2) ex, -sin(theta/2) ey, -sin(theta/2) ez, cos(theta/2)) but it's not the convention assumed here.

def convert_quaternion_eng_to_spice(qeng):
    qspice = [qeng[3], qeng[0], qeng[1], qeng[2]]
    return qspice
