# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 09:00:50 2020

@author: M06003
This function should caculate solubility SLE for a given binary pair. Activity coefficients
are calculated by cosmo sac
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import math
import cCOSMO

def calc_solubility (T,Hm, Tm):
    # calculate right hand side of solubilty accordingly to gmehling (ideal eutektikum) p 376 input is Temperature, Heat of fusion and Melting point
    #Gaskonstante, Gamma, e
    #  in J/mol/K
    R = 8.314446262
    df_a1 = -(Hm/(R*T))*(1-T/Tm)
    df_x1 = np.exp(df_a1)
    return df_x1


def calc_solubiltiy_with_gamma(prod,db, T,R_sol):
    # calculate gamma and x with iteration
        # load sigma profiles and cosmo
    # print(prod)
    COSMO = cCOSMO.COSMO3(prod, db)
    cConsts = COSMO.get_mutable_COSMO_constants()
    cConsts.fast_Gamma = False
    
    # set List for result and iter variable i
    xL = []
    i = 0
    #  loop over Temperature-Range
    gamma = 0
    x2 = R_sol
    if x2 > 0.99:
        x2 = 0.2
    eps = 1
    iter = 0
    # iterate for solubility with activity coefficients
    while eps > 1e-7:
        # get activity coefficient
        try:
            gamma = np.exp(COSMO.get_lngamma(T, [x2,1-x2]))
        except ValueError:
            print('Solubility Calc',x2,gamma,R_sol )
            xn = 1337
            break;
            
        # calc solubility
        xn = R_sol/gamma[0]
        # check boundary conditions of loop can be aborted
        eps = np.abs(xn-x2) 
        # print(eps)
        x2 =  xn
        iter = iter +1
        if iter > 150:
            # flag if calculation has been aborted
            xn = 1337
            break
        # print(i)
        #save result
        # xL.append(xn)
        # sum result in dataframe
    
    return xn


if __name__ == '__main__':
    
    #Daten einlesen
    path = 'C:/Users/X202722/Desktop/Prozesssynthese/12lm_sigma/COSMOSAC-master/CodeSLE/Daten4.xlsx'
    df_data = pd.read_excel(path, sheet_name = 'Tabelle1', header = 0)
    df_Gmehling = pd.read_excel(path, sheet_name = 'Tabelle3', header = 0)
    solvents = pd.read_csv("../profiles2/UD/complist.txt", sep = ' ')
    solvents = solvents.iloc[1:,3].to_numpy()
    # Prepare use of COSMO-SAC powered by NIST
    db = cCOSMO.DelawareProfileDatabase("../profiles2/UD/complist.txt", 
                                        "../profiles2/UD/sigma1")
    #Temperaturbereich festlegen
    T =  range(200, 400, 10)
    
    df_T = np.asarray(T)
    
    for ii in [2,3]:
        # Enthalpy of fusion in J/mol
        Hm = df_Gmehling.iloc[ii,2]
        # Melting point in Ks
        Tm = df_Gmehling.iloc[ii,1]
        # calculate right hand of solubilty (or solubility with gamma = 1)
        R_sol = calc_solubility (df_T,Hm, Tm)
    
        all_data = pd.DataFrame()
        
        
        for solvent in solvents:
            # load sigma profiles and cosmo
            prod = [solvent, df_Gmehling.iloc[ii,0]]
            print(solvent)
        
            res_solv = calc_solubiltiy_with_gamma(prod,db, T)
            # save result in dataframe
            all_data = all_data.append(res_solv)
            
        

