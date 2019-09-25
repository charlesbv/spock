import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cygnss_read_netcdf_to_eci_observation import *

filename_list = ['../2019/099/cyg03.ddmi.s20190409-000000-e20190409-235959.l1.power-brcs.a21.d21.nc',
                 '../2019/100/cyg03.ddmi.s20190410-000000-e20190410-235959.l1.power-brcs.a21.d21.nc',
                 '../2019/101/cyg03.ddmi.s20190411-000000-e20190411-235959.l1.power-brcs.a21.d21.nc',
                 '../2019/102/cyg03.ddmi.s20190412-000000-e20190412-235959.l1.power-brcs.a21.d21.nc',
                 '../2019/103/cyg03.ddmi.s20190413-000000-e20190413-235959.l1.power-brcs.a21.d21.nc',
                 '../2019/104/cyg03.ddmi.s20190414-000000-e20190414-235959.l1.power-brcs.a21.d21.nc',
]

cygnss_read_netcdf_to_eci_observation(filename_list)
