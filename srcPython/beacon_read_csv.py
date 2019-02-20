# This script reads the csv files output by tds-bop.exe
import sys
sys.path.append('/Users/cbv/work/spock/srcPython')
import csv
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
from ecef_to_lvlh import *
def beacon_read_csv(filename):
    #filename = 'outputCygnssMar/pass_9_PRN_31.csv'
    deg2rad = np.pi/180
    irow = -1
    date = []
    with open(filename) as csvfile:
        readCSV = csv.reader(csvfile, delimiter=',')
        yy = []; mm = []; dd = []; hr = []; minute = []; sec = []; elapsed_sec = [];
        prn = []; target_lat = []; target_lon = []; target_alt = [];
        target_ecef_x = []; target_ecef_y = []; target_ecef_z = [];
        target_rx_sat_look_angle_az = []; target_rx_sat_look_angle_el = [];
        target_rx_sat_range = []
        sp_lat = []; sp_lon = []; sp_ecef_pos_x = []; sp_ecef_pos_y = [];
        sp_ecef_pos_z = []; sp_gain = []; rx_sub_sat_lat = []; rx_sub_sat_lon = [];
        rx_sat_ecef_pos_x = []; rx_sat_ecef_pos_y = []; rx_sat_ecef_pos_z = [];
        rx_sat_ecef_vel_x = []; rx_sat_ecef_vel_y = []; rx_sat_ecef_vel_z = [];
        tx_sat_ecef_pos_x = []; tx_sat_ecef_pos_y = []; tx_sat_ecef_pos_z = [];
        rx_power= [];
        for row in readCSV:
            irow = irow + 1
            if irow != 0:
                date.append(row[0] + '-' + row[1].zfill(2) + '-' + row[2].zfill(2) + \
                            'T' + row[3].zfill(2) + ':' + row[4].zfill(2) + ':' + row[5].zfill(2))
                yy.append( np.float(row[0]) )
                mm.append( np.float(row[1]) )
                dd.append( np.float(row[2]) )
                hr.append( np.float(row[3]) )
                minute.append( np.float(row[4]) )
                sec.append( np.float(row[5]) )
                elapsed_sec.append( np.float(row[6]) )
                prn.append( np.float(row[7]) )
                target_lat.append( np.float(row[8]) )
                target_lon.append( np.float(row[9]) )
                target_alt.append( np.float(row[10]) )
                target_ecef_x.append( np.float(row[11]) )
                target_ecef_y.append( np.float(row[12]) )
                target_ecef_z.append( np.float(row[13]) )
                target_rx_sat_look_angle_az.append( np.float(row[14]) )
                target_rx_sat_look_angle_el.append( np.float(row[15]) )
                target_rx_sat_range.append( np.float(row[16]) )
                sp_lat.append( np.float(row[19]) )
                sp_lon.append( np.float(row[20]) )
                sp_ecef_pos_x.append( np.float(row[21]) )
                sp_ecef_pos_y.append( np.float(row[22]) )
                sp_ecef_pos_z.append( np.float(row[23]) )
                sp_gain.append( np.float(row[27]) )
                rx_sub_sat_lat.append( np.float(row[28]) )
                rx_sub_sat_lon.append( np.float(row[29]) )
                rx_sat_ecef_pos_x.append( np.float(row[30]) )
                rx_sat_ecef_pos_y.append( np.float(row[31]) )
                rx_sat_ecef_pos_z.append( np.float(row[32]) )
                rx_sat_ecef_vel_x.append( np.float(row[33]) )
                rx_sat_ecef_vel_y.append( np.float(row[34]) )
                rx_sat_ecef_vel_z.append( np.float(row[35]) )
                tx_sat_ecef_pos_x.append( np.float(row[36]) )
                tx_sat_ecef_pos_y.append( np.float(row[37]) )
                tx_sat_ecef_pos_z.append( np.float(row[38]) )
                rx_power.append( np.float(row[42]) )

    yy = np.array(yy); mm = np.array(mm); dd = np.array(dd); hr = np.array(hr);
    minute = np.array(minute); sec = np.array(sec); elapsed_sec = np.array(elapsed_sec);
    prn = np.array(prn); target_lat = np.array(target_lat);
    target_lon = np.array(target_lon); target_alt = np.array(target_alt);
    target_ecef_x = np.array(target_ecef_x); target_ecef_y = np.array(target_ecef_y);
    target_ecef_z = np.array(target_ecef_z);
    target_rx_sat_look_angle_az = np.array(target_rx_sat_look_angle_az);
    target_rx_sat_look_angle_el = np.array(target_rx_sat_look_angle_el);
    target_rx_sat_range = np.array(target_rx_sat_range);
    sp_lat = np.array(sp_lat); sp_lon = np.array(sp_lon);
    sp_ecef_pos_x = np.array(sp_ecef_pos_x); sp_ecef_pos_y = np.array(sp_ecef_pos_y);
    sp_ecef_pos_z = np.array(sp_ecef_pos_z); sp_gain = np.array(sp_gain);
    rx_sub_sat_lat = np.array(rx_sub_sat_lat); rx_sub_sat_lon = np.array(rx_sub_sat_lon);
    rx_sat_ecef_pos_x = np.array(rx_sat_ecef_pos_x);
    rx_sat_ecef_pos_y = np.array(rx_sat_ecef_pos_y); rx_sat_ecef_pos_z = np.array(rx_sat_ecef_pos_z);
    rx_sat_ecef_vel_x = np.array(rx_sat_ecef_vel_x); rx_sat_ecef_vel_y = np.array(rx_sat_ecef_vel_y);
    rx_sat_ecef_vel_z = np.array(rx_sat_ecef_vel_z);
    tx_sat_ecef_pos_x = np.array(tx_sat_ecef_pos_x); tx_sat_ecef_pos_y = np.array(tx_sat_ecef_pos_y);
    tx_sat_ecef_pos_z = np.array(tx_sat_ecef_pos_z); rx_power= np.array(rx_power);

    return date, prn, target_lat, target_lon, target_alt, target_ecef_x, target_ecef_y, target_ecef_z, target_rx_sat_look_angle_az, target_rx_sat_look_angle_el, target_rx_sat_range, sp_lat, sp_lon, sp_ecef_pos_x, sp_ecef_pos_y, sp_ecef_pos_z, sp_gain, rx_sub_sat_lat, rx_sub_sat_lon, rx_sat_ecef_pos_x, rx_sat_ecef_pos_y, rx_sat_ecef_pos_z, rx_sat_ecef_vel_x, rx_sat_ecef_vel_y, rx_sat_ecef_vel_z, tx_sat_ecef_pos_x, tx_sat_ecef_pos_y, tx_sat_ecef_pos_z, rx_power

