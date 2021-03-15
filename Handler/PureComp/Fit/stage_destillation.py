# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:48:06 2021
Stage construction distillation
Test missing and check fenkse underwood from gorak - maybe implement from tom
@author: X202722
"""

import numpy as np

def destination_Stage_C1(x_Vu,x_Prod,S,Z, alpha):
    #this function calculates stage number according to fenske underwood. Case 1 in doc
    # S is loss of product and Z defined purity by reduction of VU -> x1/S=x2
    try:
        x_ProdS = 1-x_Vu/100
        x_ProdX = (1-x_Vu/100*Z)
        part0 = (x_ProdS*S)/(1-(1-S)*x_ProdS/x_ProdX)
        part1 = (x_ProdX)/(1-x_ProdX)*(1/part0-1)        
    except ZeroDivisionError:
        part1 = 1
        
    return np.log(part1)/np.log(alpha)


def destination_Stage_C2(x_Vu,x_Prod,S,Z, alpha):
    #this function calculates stage number according to fenske underwood.Case 2 in doc
    # S is loss of product and Z defined purity by reduction of VU -> x1/S=x2
    try:
        x_ProdS = 1-x_Vu/100
        x_ProdX = (1-x_Vu/100*Z)
        part0 = S/(1/x_ProdS-(1-S)/x_ProdX)
        part1 = (x_ProdX)/(1-x_ProdX)*(1/part0-1)  
    except ZeroDivisionError:
        part1 = 1
        
    return np.log(part1)/np.log(1/alpha)