
import numpy as np
import pandas as pd
import scipy.optimize as opt

def clapeyronFit( T,p):
    # fit function for clapeyron
    import numpy as np
    from scipy import stats
    x = 1/(T)
    y =np.log(p)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return [slope, intercept]
    
def clapeyron(para,T):
    #clapeyron function
    import numpy as np
    p = (para[0]/(T)+para[1])
    # finally convert to mbar
    p = np.exp(p)
    return p