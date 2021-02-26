# -*- coding: utf-8 -*-
"""
Created on Tue May 26 07:49:48 2020
this is only for use on HFUS
@author: X202722
"""


def make_parameter_HFUS_fit (T_sim, T_exp, method, na, M):
    import numpy as np
    # fit if possible BPT guess to 
    OffsetP = 0
    OffsetP = T_exp-T_sim
    return OffsetP

def useOffset_HFUS (OffsetP, T_sim, method, na, M):
    import numpy as np
    # use fit to create new value
    #nannoolal
    Topt = 0
    
    Topt = T_sim + OffsetP
    
        
    return Topt

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        