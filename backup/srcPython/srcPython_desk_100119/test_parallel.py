
import socket
import os, sys, re

import numpy as np

from IPython import parallel


mpi_profile = 'mpi'
rc = parallel.Client(profile=mpi_profile)
eall = rc[:]
root = rc[-1]

# #  from IPython.parallel import Client
# # rc = Client()
# # rc.ids
# # dview = rc[:]
# # @dview.remote(block=True)
# # def getpid():
# #     import os
# #     return os.getpid()

# from IPython.parallel import Client

# c = Client(profile='mpi')

# view = c[:]

# view.activate() # enable magics

# # run the contents of the file on each engine:
# view.run('psum.py')

# view.scatter('a',np.arange(16,dtype='float'))

# view['a']

# %px totalsum = psum(a)

# view['totalsum']

