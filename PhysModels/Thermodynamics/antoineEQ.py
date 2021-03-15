
import numpy as np

def antoine(para,T):
    # antoine eq, T in Kelvin
    p = para[0]-para[1]/(T+para[2])
    p = 10**p
    return p