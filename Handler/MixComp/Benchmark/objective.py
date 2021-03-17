# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:07:53 2021
@author: Luisa Peterson


This file calls the objective function needed for the optimization of the
fit parameter.

Functions overview:
    RSS = objective_UNIFAC(interaction_param_fit, Data_all, Exp, Database_all, 
                           UNIFAC_const_out, gamma_exp, interaction_param, 
                           interaction_param_where)
    RSS = objective_UNISAC(segment_area_param_fit, Data, Experiments, Database, 
                           combinatorial, segment_area_param_where, gamma_exp)
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
from Adenventure.Handler.MixComp.Benchmark.find_parameter import filter_groups
from Adenventure.PhysModels.MixComp.get_gamma_fit import UNIFAC_fit, UNISAC_fit, run_COSMOSAC
from Adenventure.PhysModels.MixComp.get_gamma_inf import UNIFAC_inf

#%% OBJECTIVE FUNCTION UNIFAC
'''
The target function calls the part of UNIFAC that depends on the interaction
parameters. With the given interaction parameters the activity coecient is cal-
culated. The result is compared with an experimentally determined activity coef-
cient. For this purpose, the error square sum is determined. If this function is
called by the optimizer, it tries to find the interaction parameters so that the sum
of squares of errors is minimal.
Input:
        - interaction_param: Parameter to be fitted
        - Data_all: contains structural information
        - Exp: experiments to be considered
        - Database: UNIFAC databse
        - Databse_all: used to reload Data
        - UNIFAC_const_out: reults from the constant part of UNIFAC
        - gamma_exp: activity coefficient from experimental solubilities
        - interaction_param: matrix of interaction param including Nones
        - interaction_param_where: localizes where the interaction params that are 
        missing are
Output:
        - RSS: Residual Sum of Square (Sum over the squared error of all experiments)
'''
def objective_UNIFAC(interaction_param_fit, Data_all, Exp, Database_all, UNIFAC_const_out, 
              gamma_exp, interaction_param, interaction_param_where):
    
    # Get gammas for interaction params that are to be fitted
    Data, Database,_ = filter_groups(Exp, Data_all, Database_all)
    Exp_UNIFAC_fit = Exp.values
    gamma_mod = list(map(UNIFAC_fit,repeat(interaction_param_fit),repeat(Data),
                         Exp_UNIFAC_fit, repeat(Database),repeat(interaction_param),
                         repeat(interaction_param_where),UNIFAC_const_out, repeat(False)))
    gamma_mod = pd.concat(gamma_mod)
    gamma_mod = np.array(gamma_mod[0])
    
    # get gamma_inf from database
    gamma_inf = Database[3]
    # get position of diagonal elements
    dim = len(gamma_inf)
    diagonal = []
    for i in range(dim):
        diagonal_i = i*(dim+1)
        diagonal.append(diagonal_i)
    # flatten matrix
    gamma_inf = gamma_inf.values.flatten()
    # delete former diagonal elements from vector
    gamma_inf = np.delete(gamma_inf, diagonal)
    # create area in which gamma_inf_mod should be
    gamma_inf_max = gamma_inf*1.5+10
    gamma_inf_min = gamma_inf*0.5-10
    
    # get gamma_inf for each set of fitted interaction_param
    gamma_inf_mod = UNIFAC_inf(Data, Database, interaction_param_fit, interaction_param, interaction_param_where)
    
    # Penalty Functions
    # compare gamma_inf from database and gamma_inf with fitted interaction parameter
    # penalty 0: ensure combination of interaction params leads to a gamma within gamma_inf 
    penalty_gamma_inf = np.zeros(len(gamma_inf))
    # punish larger values
    gamma_inf_greater_where = np.where(gamma_inf_mod > gamma_inf_max)
    penalty_gamma_inf[gamma_inf_greater_where] = (gamma_inf_mod-gamma_inf_max)[gamma_inf_greater_where]
    # punish smaller values
    gamma_inf_smaller_where = np.where(gamma_inf_mod < gamma_inf_min)
    penalty_gamma_inf[gamma_inf_smaller_where] = (gamma_inf_min-gamma_inf_mod)[gamma_inf_smaller_where]
    # sum over all penalty values
    penalty_gamma_inf = sum(penalty_gamma_inf)
    
    # penalty 1: prohibit extreme interaction parameter       
    # minimal allowed value for interaction parameter
    interaction_param_min = -2500
    # maximal allowed value for interaction parameter
    interaction_param_max = 10000
    penalty_interaction_param = np.zeros(len(interaction_param_fit))
    # punish larger values
    interaction_param_greater_where = np.where(interaction_param_fit > interaction_param_max)
    penalty_interaction_param[interaction_param_greater_where] = (interaction_param_fit-interaction_param_max)[interaction_param_greater_where]
    # punish smaller values
    interaction_param_smaller_where = np.where(interaction_param_fit < interaction_param_min)
    penalty_interaction_param[interaction_param_smaller_where] = (interaction_param_min-interaction_param_fit)[interaction_param_smaller_where]
    # sum over all penalty values
    penalty_interaction_param = sum(penalty_interaction_param)
    
    # Error between model and experimental gamma
    # add panalty functions
    Err = abs(gamma_mod-gamma_exp) + penalty_gamma_inf**2 + penalty_interaction_param
    # Sum of Squares divided
    No_exp = len(gamma_mod)    
    # get Error for every experiment
    RSS = np.sum(Err)/No_exp
    return RSS

#%% OBJECTIVE FUNCTION UNISAC
'''
The target function calls the part of UNISAC that depends on the segment area
parameters. With the given segment area parameters the activity coefficient is cal-
culated. The result is compared with an experimentally determined activity coef-
cient. For this purpose, the error square sum is determined. If this function is
called by the optimizer, it tries to find the segment area parameters so that the sum
of squared errors is minimal.
Input:
        - segment_area_param_fit: Parameter to be fitted
        - Data_all: contains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC databse
        - combinatorial: reults from the constant part of UNIFAC
        - interaction_param_where: localizes where the segment area parameter
        that are missing are
        - gamma_exp: activity coefficient from experimental solubilities
Output:
        - RSS: Residual Sum of Square (Sum over the squared error of all experiments)
'''
def objective_UNISAC(segment_area_param_fit, Data, Experiments, Database, 
              combinatorial, segment_area_param_where, gamma_exp):
    # Get gammas for interaction params that are to be fitted
    Exp_UNISAC = Experiments.values
    gamma_mod = list(map(UNISAC_fit,repeat(segment_area_param_fit),repeat(Data),
                         Exp_UNISAC, repeat(Database), combinatorial,
                         repeat(segment_area_param_where),repeat(False)))
    gamma_mod = np.array(gamma_mod)
    gamma_mod = gamma_mod[:,0] 
    
    # penalty 1: prohibit extreme interaction parameter       
    # minimal allowed value for interaction parameter
    segment_area_param_min = 0
    # maximal allowed value for interaction parameter
    segment_area_param_max = 3
    penalty = np.zeros(len(segment_area_param_fit))
    # punish larger values
    segment_area_param_greater_where = np.where(segment_area_param_fit > segment_area_param_max)
    penalty[segment_area_param_greater_where] = (segment_area_param_fit-segment_area_param_max)[segment_area_param_greater_where]
    # punish smaller values
    segment_area_param_smaller_where = np.where(segment_area_param_fit < segment_area_param_min)
    penalty[segment_area_param_smaller_where] = (segment_area_param_min-segment_area_param_fit)[segment_area_param_smaller_where]+10**3
    # sum over all penalty values
    penalty= sum(penalty)
    
    # Error between model and experimental gamma
    # add panalty functions
    Err = abs(gamma_mod-gamma_exp)+penalty
    
    # Sum of Squares divided
    No_exp = len(gamma_mod)    
    # get Error for every experiment
    RSS = np.sum(Err)/No_exp
    return RSS


#%% OBJECTIVE FUNCTION COSMOSAC
'''
The objective function calls COSMO-SAC which depends on the effective radius.
With the given molecular interaction params the activity coefficient is calculated. 
The result is compared with an experimentally determined activity coefcient.
For this purpose, the error square sum is determined. If this function is
called by the optimizer, the optimizer tries to find the molecular interaction
parameters so that the sum of squaresd errors is minimal.
Input:
    - Parameter: Effective Radius that is to be fitted
    - Experiments: Experiments to be considered
    - gamma_exp: Experimentally determined activity coefficients
    - db: active Delaware Database with relevant sigma profiles
    - Parameter_0: initial value for the Fit-Parameter
Output:
    - RSS: Residual Sum of Square (Sum over the squared error of all experiments)
'''
def objective_COSMOSAC(Parameter, Experiments, gamma_exp, db):
    # Get gammas for interaction params that are to be fitted
    Exp_to_COSMO = Experiments.values
    gamma_mod = list(map(run_COSMOSAC, repeat(Parameter), Exp_to_COSMO, repeat(db)))
    gamma_mod = np.asarray(gamma_mod)
    #gamma_mod = np.array(gamma_mod[0])
    # Param_min = np.array(Parameter_0/10)
    # Param_max = np.array(Parameter_0*10)
    if len(Parameter)==3:
        # Fit chb
        param_min = np.array([401.378, 92.332, 301.643])
        param_max = np.array([10000.8, 9233.2, 30164.3])
    if len(Parameter)==2:
        # Fit electrostatic
        param_min = np.array([6525.690000,1.4859*10**7])
        param_max = np.array([652569.0000,1.4859*10**9])

    param_fit = Parameter
    penalty = np.zeros(len(param_fit))
    # punish greater values
    param_greater_where = np.where(param_fit > param_max)
    penalty[param_greater_where] = (param_fit-param_max)[param_greater_where]
    # punish smaller values
    param_smaller_where = np.where(param_fit < param_min)
    penalty[param_smaller_where] = (param_min-param_fit)[param_smaller_where]   
    Err = abs(gamma_mod-gamma_exp)+sum(penalty)
    # Sum of Squares divided
    No_exp = len(gamma_mod)    # get Error for every experiment
    RSS = np.sum(Err)/No_exp
    return RSS