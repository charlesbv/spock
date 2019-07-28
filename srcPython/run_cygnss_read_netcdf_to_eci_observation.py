import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cygnss_read_netcdf_to_eci_observation import *
import ipdb
filename_list = ['244/cyg08.ddmi.s20170901-000000-e20170901-235959.l1.power-brcs.a21.d21.nc']


# ['229/cyg01.ddmi.s20170817-000000-e20170817-235959.l1.power-brcs.a21.d21.nc',
#                  '230/cyg01.ddmi.s20170818-000000-e20170818-235959.l1.power-brcs.a21.d21.nc',
#                  '231/cyg01.ddmi.s20170819-000000-e20170819-235959.l1.power-brcs.a21.d21.nc',
#                  '232/cyg01.ddmi.s20170820-000000-e20170820-235959.l1.power-brcs.a21.d21.nc',
#                  '233/cyg01.ddmi.s20170821-000000-e20170821-235959.l1.power-brcs.a21.d21.nc',
#                  '234/cyg01.ddmi.s20170822-000000-e20170822-235959.l1.power-brcs.a21.d21.nc',
#                  '235/cyg01.ddmi.s20170823-000000-e20170823-235959.l1.power-brcs.a21.d21.nc',
#                  '236/cyg01.ddmi.s20170824-000000-e20170824-235959.l1.power-brcs.a21.d21.nc',
#                  '237/cyg01.ddmi.s20170825-000000-e20170825-235959.l1.power-brcs.a21.d21.nc',
#                  '238/cyg01.ddmi.s20170826-000000-e20170826-235959.l1.power-brcs.a21.d21.nc',
# ]

cygnss_read_netcdf_to_eci_observation(filename_list)
