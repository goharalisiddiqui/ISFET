import numpy as np
import math


def loggrid(L,min_spacing):
    dm = min_spacing
    ff=1e-9
    lndm = math.log((dm/ff)+1.0)
    numer = math.log(L+ff)-math.log(ff)
    NP = int(numer/lndm)+2
    x = np.logspace(math.log(ff),math.log(L+ff),NP,base=math.e)-ff
    return x
