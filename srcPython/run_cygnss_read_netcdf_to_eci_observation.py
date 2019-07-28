import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cygnss_read_netcdf_to_eci_observation import *
<<<<<<< HEAD
filename_list = ['244/cyg08.ddmi.s20170901-000000-e20170901-235959.l1.power-brcs.a21.d21.nc',
                 '245/cyg08.ddmi.s20170902-000000-e20170902-235959.l1.power-brcs.a21.d21.nc',
                 '246/cyg08.ddmi.s20170903-000000-e20170903-235959.l1.power-brcs.a21.d21.nc',
                 '247/cyg08.ddmi.s20170904-000000-e20170904-235959.l1.power-brcs.a21.d21.nc',
                 '248/cyg08.ddmi.s20170905-000000-e20170905-235959.l1.power-brcs.a21.d21.nc',
                 '249/cyg08.ddmi.s20170906-000000-e20170906-235959.l1.power-brcs.a21.d21.nc',
                 '250/cyg08.ddmi.s20170907-000000-e20170907-235959.l1.power-brcs.a21.d21.nc',
                 '251/cyg08.ddmi.s20170908-000000-e20170908-235959.l1.power-brcs.a21.d21.nc',
                 '252/cyg08.ddmi.s20170909-000000-e20170909-235959.l1.power-brcs.a21.d21.nc',
                 '253/cyg08.ddmi.s20170910-000000-e20170910-235959.l1.power-brcs.a21.d21.nc',
                 '254/cyg08.ddmi.s20170911-000000-e20170911-235959.l1.power-brcs.a21.d21.nc'
]
=======
import ipdb
filename_list = ['244/cyg08.ddmi.s20170901-000000-e20170901-235959.l1.power-brcs.a21.d21.nc']
>>>>>>> b85f21d3ac71e35cdfaf355470dee7eaea97494e


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
