import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cygnss_read_netcdf_to_eci_observation import *

filename_list = ['244/cyg03.ddmi.s20180901-000000-e20180901-235959.l1.power-brcs.a21.d21.nc',
                 '245/cyg03.ddmi.s20180902-000000-e20180902-235959.l1.power-brcs.a21.d21.nc',
                 '246/cyg03.ddmi.s20180903-000000-e20180903-235959.l1.power-brcs.a21.d21.nc',
                 '247/cyg03.ddmi.s20180904-000000-e20180904-235959.l1.power-brcs.a21.d21.nc',
                 '248/cyg03.ddmi.s20180905-000000-e20180905-235959.l1.power-brcs.a21.d21.nc',
                 '249/cyg03.ddmi.s20180906-000000-e20180906-235959.l1.power-brcs.a21.d21.nc',
                 '250/cyg03.ddmi.s20180907-000000-e20180907-235959.l1.power-brcs.a21.d21.nc',
                 '251/cyg03.ddmi.s20180908-000000-e20180908-235959.l1.power-brcs.a21.d21.nc',
                 '252/cyg03.ddmi.s20180909-000000-e20180909-235959.l1.power-brcs.a21.d21.nc',
]

cygnss_read_netcdf_to_eci_observation(filename_list)
