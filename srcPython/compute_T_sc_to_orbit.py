import numpy as np
import sys
def compute_T_sc_to_orbit(v_angle, order_rotation): # v_angle is [pitch, roll, yaw] (in deg), order_rotation order 1 means firt you do this rotation. for example: [2,1,3] -> first you roll then pitch then yaw
# # !!!!!!!!!! REMOVE BLOCK
# v_angle = [0,0,-90] # pitch roll yaw
# order_rotation = [3,2,1] # order 1 means firt you do this rotation. for example: [2,1,3] -> first you roll then pitch then yaw
# # end of REMOVE BLOCK


    T_sc_to_orbit = np.zeros([3,3])


    dtor = np.pi / 180
    v_angle_rad = np.zeros([3])
    for i in range(3):
        v_angle_rad[i] = v_angle[i] *  dtor;


    theta = v_angle_rad[0]; # pitch
    phi = v_angle_rad[1]; # roll
    psi = v_angle_rad[2]; # yaw

    pitch_mat = np.zeros([3,3]); roll_mat =  np.zeros([3,3]); yaw_mat = np.zeros([3,3])

    pitch_mat[0, 0] = np.cos(theta);   pitch_mat[0, 1] = 0;   pitch_mat[0, 2] = np.sin(theta);
    pitch_mat[1, 0] = 0;   pitch_mat[1, 1] = 1;   pitch_mat[1, 2] = 0;
    pitch_mat[2, 0] = -np.sin(theta);   pitch_mat[2, 1] = 0;   pitch_mat[2, 2] = np.cos(theta);

    # Roll rotation
    roll_mat[0, 0] = 1;   roll_mat[0, 1] = 0;   roll_mat[0, 2] = 0;
    roll_mat[1, 0] = 0;   roll_mat[1, 1] = np.cos(phi);   roll_mat[1, 2] = -np.sin(phi);
    roll_mat[2, 0] = 0;   roll_mat[2, 1] = np.sin(phi);   roll_mat[2, 2] = np.cos(phi);

    # Yaw rotation
    yaw_mat[0, 0] = np.cos(psi);   yaw_mat[0, 1] = -np.sin(psi);   yaw_mat[0, 2] = 0;
    yaw_mat[1, 0] = np.sin(psi);   yaw_mat[1, 1] = np.cos(psi);   yaw_mat[1, 2] = 0;
    yaw_mat[2, 0] = 0;   yaw_mat[2, 1] = 0;   yaw_mat[2, 2] = 1;

    # Order of rotation
    # # First matrix
    if (order_rotation[0] == 1):
        first_mat = pitch_mat
    elif (order_rotation[1] == 1):
        first_mat = roll_mat
    elif (order_rotation[2] == 1):
        first_mat = yaw_mat
    else:
        print "***! You did not correctly choose the first rotation. The program will stop. !***\n"; sys.exit()


    # # Second matrix
    if (order_rotation[0] == 2):
        second_mat = pitch_mat
    elif (order_rotation[1] == 2):
        second_mat = roll_mat
    elif (order_rotation[2] == 2):
        second_mat = yaw_mat
    else:
        print "***! You did not correctly choose the second rotation. The program will stop. !***\n"; sys.exit()

    # # Third matrix
    if (order_rotation[0] == 3):
        third_mat = pitch_mat
    elif (order_rotation[1] == 3):
        third_mat = roll_mat
    elif  (order_rotation[2] == 3):
        third_mat = yaw_mat
    else:
        print "***! You did not correctly choose the third rotation. The program will stop. !***\n"; sys.exit()

    # # Final matrix: T_sc_to_orbit
    T_sc_to_orbit_temp = np.matmul(second_mat, first_mat)
    T_sc_to_orbit = np.matmul( third_mat, T_sc_to_orbit_temp)

    return T_sc_to_orbit
