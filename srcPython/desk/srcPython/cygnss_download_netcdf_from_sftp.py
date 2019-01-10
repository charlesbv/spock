# This scriptdownloads L1 netcdf files from cygnss-sftp-2.engin.umich.edu:/data/from_soc/l1
# inputs: date_start, date_stop, cygfm, save_dir


import os
from datetime import datetime, timedelta

date_start = '2018-09-20'
date_stop = '2018-09-30'
for cygfm in range(1, 9):
    #cygfm = 1 # 1 to 8
    save_dir = '/Users/cbv/cygnss/netcdf_sand39/' # with slash at the end. the year will be added as a sub-directory


    date_start_date = datetime.strptime(date_start, "%Y-%m-%d")
    date_stop_date = datetime.strptime(date_stop, "%Y-%m-%d")


    date_here = date_start_date
    while date_here <= date_stop_date:
        yy = date_here.strftime('%Y')
        doy = date_here.strftime('%j')

        if (os.path.isdir(save_dir + yy + '/' + str(doy).zfill(3)) == False):
            os.system("mkdir " + save_dir + yy + '/' + str(doy).zfill(3))
        # print "scp -p cygnss-sftp-2.engin.umich.edu:/data/from_soc/l1/" + yy + '/'  + str(doy).zfill(3) + "/cyg0" + str(cygfm) + "*37.nc " + save_dir + yy + '/'  + str(doy).zfill(3)
        # os.system("scp -p cygnss-sftp-2.engin.umich.edu:/data/from_soc/l1/" + yy + '/'  + str(doy).zfill(3) + "/cyg0" + str(cygfm) + "*37.nc " + save_dir + yy + '/'  + str(doy).zfill(3))
        print "scp -p cygnss-sftp-1.engin.umich.edu:/data/temp/v2_1_preview/l1/" + yy + '/'  + str(doy).zfill(3) + "/cyg0" + str(cygfm) + "*39.nc " + save_dir + yy + '/'  + str(doy).zfill(3)
        os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/temp/v2_1_preview/l1/" + yy + '/'  + str(doy).zfill(3) + "/cyg0" + str(cygfm) + "*39.nc " + save_dir + yy + '/'  + str(doy).zfill(3))

        highest_sand = -1
        for file in os.listdir(save_dir + yy + '/'  + str(doy).zfill(3)):
            sand = (int)(file.split('.sand')[1].split('.')[0])
            if sand >= highest_sand:
                highest_sand = sand
        for file in os.listdir(save_dir + yy + '/'  + str(doy).zfill(3)):
            sand = (int)(file.split('.sand')[1].split('.')[0])
            if sand != highest_sand:
                #print "rm " + save_dir + yy + '/'  + str(doy).zfill(3) + '/' + file
                os.system("rm " + save_dir + yy + '/'  + str(doy).zfill(3) + '/' + file)


        date_here = date_here + timedelta(days = 1)



    #         for iday in range(nb_day):
    #             if (os.path.isdir("/Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3)) == False):
    #                 os.system("mkdir /Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3))
    #             if download_netcdf == 1:
    #                 os.system("scp -p cygnss-sftp-1.engin.umich.edu:/data/cygnss/products/l1/" + datetime.strptime( date_stop_val, "%Y-%m-%dT%H:%M:%S" ).strftime('%Y') + " /" + str(doy_array[iday]).zfill(3) + "/cyg0" + str(cygfm) + "* /Users/cbv/cygnss/netcdf/" + yy  + str(doy_array[iday]).zfill(3))

