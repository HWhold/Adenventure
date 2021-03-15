# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 23:39:32 2021
This is the function for the gmehling algorith as described in the paper doi/abs/10.1021/
@author: X202722
Vorgehen:
    - Ghemling function schreiben und validieren
    - Gmehling function in Pandas abbilden
"""
import itertools
from itertools import repeat
from Adenventure.PhysModels.Thermodynamics.BinaryActivityCoefficients import CalcBinary
import cCOSMO
import numpy as np
import pandas as pd
# Hier Annahme Produkt geht immer ins Raffinat. Der Rest soll im Extrakt bleiben. Test werden über Iter abgehandelt, sodass auch Mehrkomponentenmischungen abgebildet werden können

def calc_KPI_gmehling(K_E_i, K_R_j, x_SR, alpha_dis, vec ):
    #Four parts first capacaity
    # solvent loss
    # this is the validation with data from gmeling

    ws = vec[0]
    wc = vec[1]
    wsl = vec[2]
    wsd = vec[3]
    
    #number of components
    num_components = len(K_E_i)+1
    #capacity = first part of formula
    capacity =  sum(K_E_i/K_R_j)/(num_components-1)
    
    # selectivity
    selectivity = sum(K_E_i)/len(K_E_i)
    
    # solvent loss of raffinat solvet
    solvent_loss = 1-x_SR
    
    #  distillation difficulty
    distillation = np.log10(alpha_dis)
    
    return capacity**ws * selectivity**wc * solvent_loss**wsl * distillation**wsd




if __name__ == '__main__':
    #Four parts first capacaity
    # solvent loss
    # this is the validation with data from gmeling
    x_SR = 0.01
    alpha_dis = 1/0.46
    ws = 1
    wc = 1
    wsl = 3
    wsd = 0.75
    # distribution coefficient of components in extract =gamma_LM1/gamma_LM2
    K_E_i = np.array([0.4843,	0.3625])
    # distribution coefficient of components in raffinat - here product
    K_R_j = 0.0106
    
    #number of components
    num_components = len(K_E_i)+1
    #capacity = first part of formula
    capacity =  sum(K_E_i/K_R_j)/(num_components-1)
    
    # selectivity
    selectivity = sum(K_E_i)/len(K_E_i)
    
    # solvent loss of raffinat solvet
    solvent_loss = 1-x_SR
    
    #  distillation difficulty
    distillation = np.log10(alpha_dis)
    
    res = capacity**ws * selectivity**wc * solvent_loss**wsl * distillation**wsd
    vec = [ws, wc, wsl, wsd]
    compare = res - calc_KPI_gmehling(K_E_i, K_R_j, x_SR, alpha_dis, vec )
    
    print(compare)
