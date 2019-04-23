import sys
sys.path.append("/Users/cbv/work/spock/srcPython")
from cygnss_read_netcdf_to_eci_observation import *
filename_list = ['013/cyg03.ddmi.s20180113-000000-e20180113-235959.l1.power-brcs.a21.d21.nc',
                 '014/cyg03.ddmi.s20180114-000000-e20180114-235959.l1.power-brcs.a21.d21.nc',
                 '015/cyg03.ddmi.s20180115-000000-e20180115-235959.l1.power-brcs.a21.d21.nc',
                 '016/cyg03.ddmi.s20180116-000000-e20180116-235959.l1.power-brcs.a21.d21.nc',
                 '017/cyg03.ddmi.s20180117-000000-e20180117-235959.l1.power-brcs.a21.d21.nc',
                 '018/cyg03.ddmi.s20180118-000000-e20180118-235959.l1.power-brcs.a21.d21.nc',
                 '019/cyg03.ddmi.s20180119-000000-e20180119-235959.l1.power-brcs.a21.d21.nc',
                 '020/cyg03.ddmi.s20180120-000000-e20180120-235959.l1.power-brcs.a21.d21.nc',
                 '021/cyg03.ddmi.s20180121-000000-e20180121-235959.l1.power-brcs.a21.d21.nc',
                 '022/cyg03.ddmi.s20180122-000000-e20180122-235959.l1.power-brcs.a21.d21.nc',
]

cygnss_read_netcdf_to_eci_observation(filename_list)
