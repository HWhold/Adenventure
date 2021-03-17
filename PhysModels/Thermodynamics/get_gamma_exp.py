# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:49:11 2021
@author: Luisa Peterson

This file loads the function for getting the experimental activity coeffcients.

Functions overview:
    gamma_exp = activity_coefficient_from_solubility(Data, Experiments)
"""

# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# builds and contains data structures and operators for accessing numeric tables
import pandas as pd
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np

# %% EXPERIMENTAL ACTIVITY COEFFICIENT
'''
This function calculates the activity coefficients with given solubilites,
melting points and enthalpy of fusion. The calculation is done based on
Gmehling S. 395 Eutectic Systems
Input: 
    - Data: Melting Point and Enthalpy of Fusion for the molecules
    - Experiments: Experiments to be considered
Output:
    - gamma_exp: experimental activity coefficients
'''
def activity_coefficient_from_solubility(Data, Experiments):
    # load data for one solute and all considered solvents
    # Solute Data
    solutes_fit = Experiments['Solute (1)'].unique()
    solutes = Data[0].copy()
    solutes = solutes.loc[solutes_fit]
    # Experimental Data
    Experiments = Experiments.copy()
    # get solute name for experiment
    solute_name = Experiments.loc[:, 'Solute (1)']
    # match solvent table to experiments
    solute_exp = solutes.loc[solute_name,:]
    # Temperature
    Temperature = np.array(Experiments.loc[:,'Temperature'])
    # Solubility date for solutes and solvents
    solubility_solute = Experiments.loc[:,'x_1']
    solubility_solute.index = np.arange(len(solubility_solute))
    # ideal gas constant
    R = 8.31446261815324
    # convert into np.arrays
    MPT = np.array(solute_exp['MPT'])
    HFUS = np.array(solute_exp['HFUS'])
    cp = HFUS/MPT
    # Calculate activity for the one solute
    e_length = len(Experiments)
    gamma_exp_solute = []
    for e in range(e_length):
        activity_cp = -HFUS[e]/(R*Temperature[e])*(1-Temperature[e]/(MPT[e])) + (cp[e]/(R*Temperature[e]))*(MPT[e]-Temperature[e])-(cp[e]/R)*np.log(MPT[e]/Temperature[e])
        # print((cp[e]/(R*Temperature[e]))*(MPT[e]-Temperature[e])-(cp[e]/R)*np.log(MPT[e]/Temperature[e]))
        # if necesarry: add transition enthalphy
        if solute_exp['TRST'].iloc[e] != 'None' and Temperature[e] < (solute_exp['TRST'].iloc[e]):
            TRST = np.array(solute_exp['TRST'])
            HTRS = np.array(solute_exp['HTRS'])
            transition = -HTRS[e]/(R*Temperature[e])*(1-Temperature[e]/(TRST[e]))+ (cp[e]/(R*Temperature[e]))*(TRST[e]-Temperature[e])-(cp[e]/R)*np.log(TRST[e]/Temperature[e])
            activity_cp = activity_cp + transition
        # get activity coefficient by dividing activity by measured solubility
        gamma_exp_solute_e = np.exp(activity_cp)/solubility_solute[e]
        gamma_exp_solute_e = pd.Series(gamma_exp_solute_e)
        gamma_exp_solute.append(gamma_exp_solute_e)
    gamma_exp_solute = pd.concat(gamma_exp_solute)
    # convert into numpy array
    gamma_exp = np.array(gamma_exp_solute)
    return gamma_exp

# GIVEN GAMMA EXP
def load_exp_gamma(Experiments):
    '''
    This function loads activity coefficients
    Input: 
        - Experiments: Experiments to be considered
    Output:
        - gamma_exp: experimental activity coefficients
    '''
    gamma_exp = Experiments.iloc[:,-1]
    gamma_exp = np.array(gamma_exp)
    return gamma_exp