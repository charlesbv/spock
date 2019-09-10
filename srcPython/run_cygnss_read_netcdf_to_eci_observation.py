import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cygnss_read_netcdf_to_eci_observation import *

filename_list = ['326/cyg03.ddmi.s20181122-000000-e20181122-235959.l1.power-brcs.a21.d21.nc',
                 '327/cyg03.ddmi.s20181123-000000-e20181123-235959.l1.power-brcs.a21.d21.nc',
                 '328/cyg03.ddmi.s20181124-000000-e20181124-235959.l1.power-brcs.a21.d21.nc',
                 '329/cyg03.ddmi.s20181125-000000-e20181125-235959.l1.power-brcs.a21.d21.nc',
                 '330/cyg03.ddmi.s20181126-000000-e20181126-235959.l1.power-brcs.a21.d21.nc',
                 '331/cyg03.ddmi.s20181127-000000-e20181127-235959.l1.power-brcs.a21.d21.nc',
                 '332/cyg03.ddmi.s20181128-000000-e20181128-235959.l1.power-brcs.a21.d21.nc',
                 '333/cyg03.ddmi.s20181129-000000-e20181129-235959.l1.power-brcs.a21.d21.nc',
                 '334/cyg03.ddmi.s20181130-000000-e20181130-235959.l1.power-brcs.a21.d21.nc',
]

cygnss_read_netcdf_to_eci_observation(filename_list)
