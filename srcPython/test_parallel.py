# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at

#   http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.


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

