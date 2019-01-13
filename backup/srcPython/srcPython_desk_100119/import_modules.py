import sys
from out import *
from matplotlib import pyplot as plt

if 'out' in sys.argv:
    import out
    reload(out)
    from out import *
