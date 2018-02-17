import os
import socket
from mpi4py import MPI

MPI = MPI.COMM_WORLD

print("{host}[{pid}]: {rank}/{size}".format(
    host=socket.gethostname(),
    pid=os.getpid(),
    rank=MPI.rank,
    size=MPI.size,
))
