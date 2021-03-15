# -*- coding: utf-8 -*-
"""
Created on Tue May 26 07:49:48 2020
this is only for use on MPT
@author: X202722
"""


def make_parameter_MPT_fit (T_sim, T_exp, method, na, M):
    import numpy as np
    # fit if possible BPT guess to 
    OffsetP = 0
    if method == 0:
        #joback
        OffsetP = T_exp-T_sim
    elif method == 1:
        #GVS
        OffsetP = T_exp-T_sim
    elif method == 2:
        #gani
        k = 102.425
        OffsetP = np.exp(T_exp/k) -np.exp(T_sim/k)
        # OffsetP = T_exp-T_sim
   
        
    return OffsetP

def useOffset_MPT (OffsetP, T_sim, method, na, M):
    import numpy as np
    # use fit to create new value
    #nannoolal
    Topt = 0
    
    if method == 0:
        #joback
        Topt = T_sim + OffsetP
    elif method == 1:
        #GVS
        Topt = T_sim + OffsetP
    elif method == 2:
        #gani
        k = 102.425
        value = np.exp(T_sim/k)
        # print(value + OffsetP, 'Gani shouts help')
        # if value + OffsetP<0:
        #     print(value, 'Gani shouts help')
        #     print(OffsetP)
            
        Topt = k * np.log(value + OffsetP)
        # Topt = T_sim + OffsetP
        # Topt = T_sim + OffsetP
   
        
    return Topt

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        