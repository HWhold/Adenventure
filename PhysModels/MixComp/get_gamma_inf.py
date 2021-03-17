# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:03:55 2021
@author: Luisa Peterson

This file calls the function that calculates gamma at infinite dilutions 
beween the main groups of UNIFAC.

Functions overview:
    gamma_inf = UNIFAC_inf(Data, Database, interaction_param_fit, 
                           interaction_param, interaction_param_where)
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
# For mapping function to repeat values
from itertools import repeat

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.PhysModels.MixComp.get_gamma_const import UNIFAC_const
from Adenventure.PhysModels.MixComp.get_gamma_fit import UNIFAC_fit

#%% GAMMA AT INFINITE DILUTION
'''
With the considered interaction parameters, this function calculates the
activity coefficient at infinite dilution. The interaction between two groups
where the solute is 0 and the solvent is 1 is considered. Therefore the DataFrames 
structural information and Experiments are adapted. 
Input:
    - Data: contains structural information
    - Database: UNIFAC databae
    - interaction_param_fit: interaction parameter that are not given and, thus,
    to be fitted
    - interaction_param: matrix of interaction param including Nones
    - interaction_param_where: localizes where the interaction params that are 
    to be fitted are
Output:
    - gamma_inf: activity coefficient at infinite dilutions
'''
def UNIFAC_inf(Data, Database, interaction_param_fit, interaction_param, interaction_param_where):
    # structural information: in order to get the infinity coefficient at infinie dilusion
    # each group is taken as a molecule
    struct_info = Data[2].copy()
    column_name = list(struct_info.columns)
    struct_info_gamma_inf = pd.DataFrame(np.zeros((len(struct_info.T),len(struct_info.T))))
    # every value on the diagonal is set to 1
    np.fill_diagonal(struct_info_gamma_inf.values, 1)
    # set index and columns
    struct_info_gamma_inf = struct_info_gamma_inf.set_axis(column_name, axis='index')
    struct_info_gamma_inf = struct_info_gamma_inf.set_axis(column_name, axis='columns')
    Data[2] = struct_info_gamma_inf
    
    # Build experiments: every group that exists in one of the present components interacts 
    # with all other groups of those components
    Exp = pd.DataFrame(columns=['Temperature','Solute (1)','Solvent (2)', 'x_1', 'x_2'])
    for i in range(len(struct_info.T)):
        for j in range(len(struct_info.T)):
            if i !=j:
                solute = column_name[i]
                solvent = column_name[j]
                Exp_ij = {'Temperature': 373.15, 'Solute (1)': solute, 'Solvent (2)': solvent, 'x_1': 0, 'x_2':1}
                Exp_ij = pd.DataFrame(Exp_ij, index = [0])
                Exp = pd.concat([Exp,Exp_ij])
                
    # call UNIFAC (UNIFAC_const + UNIFAC_fit) to get activity coefficient at infinite dilution
    Exp = Exp.values
    UNIFAC_const_out = list(map(UNIFAC_const, repeat(Data), Exp, repeat(Database), repeat(False))) 
    gamma_inf = list(map(UNIFAC_fit,repeat(interaction_param_fit),repeat(Data),
                         Exp, repeat(Database),repeat(interaction_param),
                         repeat(interaction_param_where),UNIFAC_const_out, repeat(False)))
    gamma_inf = pd.concat(gamma_inf)
    gamma_inf = np.array(gamma_inf[0])             
    return gamma_inf