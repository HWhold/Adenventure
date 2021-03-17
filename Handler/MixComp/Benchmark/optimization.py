# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:15:12 2021
@author: Luisa Peterson

This file calls the optimizer that minimizes the objective function in order to
find the best fitting fit-parameter.

Functions overview:
    solution = optimization(inputs, interaction_param_fit_0, bnds, solver)
"""

# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
# provides functions for minimizing (or maximizing) objective_UNIFAC functions.
from scipy.optimize import  minimize

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Handler.MixComp.Benchmark.objective import objective_UNIFAC, objective_UNISAC, objective_COSMOSAC

# %% OPTIMIZATION
'''
Set Up of the Optimization Problem:
Initial values and boundary conditions that are given to the optimizer.
Input:
    - Vector of inputs:
        - Data
        - Experiments
        - Database
        - UNIFAC_const_out
        - gamma_exp
        - interaction_param
        - interaction_param_where
    - interaction_param_fit_0: initial values for the solver
    - bnds: boundaries on the fitting parameters
    - solver: selection of the solver
Output:
    - solution: output of the optimizer including fitted interaction parameters,
    final Residual Sum of Squares and stats from the optimizer
'''
def optimization(inputs, fit_0, solver, model):
    print('I was here')
    if model == 'UNIFAC':    
        # RUN OPTIMIZATION
        Data_all, Exp, Database_all, UNIFAC_const_out, gamma_fit_exp, interaction_param, interaction_param_where = inputs
        # vector of arguments
        args = (Data_all, Exp, Database_all, UNIFAC_const_out, 
                  gamma_fit_exp, interaction_param, interaction_param_where)
        objective = objective_UNIFAC
        
    if model == 'UNISAC':
        # RUN OPTIMIZATION
        Data, Exp, Database, combinatorial, segment_area_param_where, gamma_exp = inputs
        # vector of arguments
        args = (Data, Exp, Database, combinatorial, segment_area_param_where, gamma_exp)
        objective = objective_UNISAC
        
    if model == 'COSMOSAC':
        Experiments_fit, gamma_fit_exp, db = inputs
        args = (Experiments_fit, gamma_fit_exp, db)
        objective = objective_COSMOSAC
        
    # different solvers can be chosen by function input
    if solver == 'BFGS':
        ## MINIMIZE PARALLEL
        fit_0 = np.array(fit_0)
        solution = minimize(objective, fit_0, method='BFGS', args = args)
                            #options={'maxiter': 2})
    
    if solver == 'Simplex':
        ## SIMPLEX
        fit_0 = np.array(fit_0)
        solution =  minimize(objective, fit_0, method='nelder-mead', args = args)
                             #options={'maxiter': 5, 'maxfev': 5, 'disp': True})
    print('I was being called')
    
    return solution