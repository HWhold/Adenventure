# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:33:25 2021
@author: M294241

Testing the implementation of the UNIFAC model which is implemented in two parts:
UNIFAC_const+UNIFAC_fit. In the following test, UNIFAC calculates activity
coefficients for the system n-Hexane in 2-Butanone. For comparison the activity
coefficients were obtained by the Dortmund Database for the same experimental
conditions. The activity coefficients calculated from the implemented UNIFAC and
the activity coefficients from the DDBSP are compared.
"""
# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
import pytest
# builds and contains data structures and operators for accessing numeric tables
import pandas as pd
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
# For mapping function to repeat values
from itertools import repeat, starmap
# Factorial Sampling
from pyDOE import lhs
# import COSMO
import cCOSMO

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Test.g_Test_COSMOSAC_Fit.g_set_up_test import load_Experiments
from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import load_exp_gamma
from Adenventure.PhysModels.MixComp.get_gamma_fit import run_COSMOSAC
from Adenventure.Handler.MixComp.Benchmark.optimization import optimization

def test_COSMOSAC_chb_Handler():
    '''
    This function runs all subfunctions in the right order
    Input:
            - seed: check reproducibility at different set of starting values
            - name group: name to load the right excel tables
    Output:
            - out_trainingsset:
                - Experiments_Error_fit: Table of Experiments with activity coefficient
                - RSS_fit: residual sum of squares for original UNISAC and UNISAC with own groups
                - Stats_fit: Table of statistical evaluation (mean value, standard deviation)
                - db: active COSMO databse
                - Parameter_opti: optimized chb parameters
    '''
    # Load Data
    seed = np.random.randint(1e9)
    name_group = 'test'
    model = 'COSMOSAC'
    # Load Experiments
    Experiments_fit = load_Experiments(name_group,'fit')
    Exp_to_COSMO = Experiments_fit.values
    
    # Calculate activity coefficients from experimentally determined solubilities
    gamma_fit_exp = load_exp_gamma(Experiments_fit)
    
    # Load COSMO Database
    #here = os.path.abspath(os.path.dirname(__file__))
    here = "c:/users/m294241/desktop/GitEntwicklung/Adenventure/Handler/MixComp/Benchmark/Data/COSMOSAC"
    db = cCOSMO.DelawareProfileDatabase(
    # table for translation
    here+"/profiles/UD/complist.txt", 
    here+"/profiles/UD/sigma3/")
    
    # add sigma profiles for all solute and solvent molecules to active database
    names = Experiments_fit['Solute (1)'].unique()  
    names = np.append(names, Experiments_fit['Solvent (2)'].unique())    
    for iden in names:
        # Delaware
        db.add_profile(db.normalize_identifier(iden))  
    
    ## OPTIMIZATION
    # Set Up Optimization
    # number of parameters to be optimized (=dimension of the optimization problem)
    dim = 3
    # number of sets of initial values
    num = dim**2
    
    # bounds on interaction parameters (one magnitude within original value)
    # interaction between OH hydroxyl groups
    bnds_min = []
    bnds_max = []
    bnds = []
    bnds_min_c_OH_OH = 401.378
    bnds_min.append(bnds_min_c_OH_OH)
    bnds_max_c_OH_OH = 40137.8
    bnds_max.append(bnds_max_c_OH_OH)
    bnds_new = [(bnds_min_c_OH_OH, bnds_max_c_OH_OH)]
    bnds_new = pd.Series(bnds_new)
    bnds.append(bnds_new)
    
    # interaction between other hydrogen bonding groups (e.g. H bound to O, N and F)
    bnds_min_c_OT_OT = 92.332
    bnds_min.append(bnds_min_c_OT_OT)
    bnds_max_c_OT_OT  = 9233.2
    bnds_max.append(bnds_max_c_OT_OT)
    bnds_new = [(bnds_min_c_OT_OT , bnds_max_c_OT_OT)]
    bnds_new = pd.Series(bnds_new)
    bnds.append(bnds_new)
    
    # interaction between OH hydroxyl group and other hydrogen bonding group
    bnds_min_c_OH_OT = 301.643
    bnds_min.append(bnds_min_c_OH_OT)
    bnds_max_c_OH_OT  = 30164.3
    bnds_max.append(bnds_max_c_OH_OT)
    bnds_new = [(bnds_min_c_OH_OT , bnds_max_c_OH_OT)]
    bnds_new = pd.Series(bnds_new)
    bnds.append(bnds_new)
    
    # convert into numpy array and concatenate
    bnds_min = np.array(bnds_min)
    bnds_max = np.array(bnds_max)
    bnds = pd.concat(bnds)
    bnds = tuple(bnds)
    
    # create semi-random initial values via LHS
    np.random.seed(seed)
    LHS_design = lhs(dim, samples=num)
    slope = bnds_max - bnds_min
    offset = bnds_min
    SLOPE = np.ones((dim,num))
    OFFSET = np.ones((dim,num))
    for i in range(num):
        SLOPE[:,i] = np.ones(dim)*slope
        OFFSET[:,i] = np.ones(dim)*offset
    Parameter = SLOPE.T*LHS_design+OFFSET.T
    obj = 1e10

    # get parameter as inputs for starmap
    inputs = [Experiments_fit, gamma_fit_exp, db]
    params = []
    for Para in Parameter:
        params.append(tuple((inputs, Para, 'BFGS', model)))
    
    # Run optimization with the BFGS solver
    solutions_BFGS_starmap = starmap(optimization, params)
    solutions_BFGS = list(solutions_BFGS_starmap)
    
    # order solutions according to their RSS and choose the best 20%
    Parameter_BFGS = []
    for solution in solutions_BFGS:
        Parameter_BFGS.append([solution.x, solution.fun])
    Parameter_BFGS = np.concatenate(Parameter_BFGS)
    Parameter_BFGS = np.reshape(Parameter_BFGS, (-1, 2))
    Parameter_BFGS = Parameter_BFGS[Parameter_BFGS[:,1].argsort()]
    Parameter_BFGS = np.delete(Parameter_BFGS,1,1)
    num_simplex = np.ceil(0.2*num)
    Parameter_BFGS = Parameter_BFGS[0:int(num_simplex),:]

    # get parameter as inputs for starmap including solution of the BFGS optimization
    # as the input to the simplex optimizer
    params = []
    for Para in Parameter_BFGS:
        Para = np.vstack(Para).astype(np.float)
        inputs = [Experiments_fit, gamma_fit_exp, db]
        params.append(tuple((inputs, Para, 'Simplex', model)))
        
    # Run optimization with the simplex solver
    solutions_Simplex_starmap = starmap(optimization, params)
    solutions_Simplex = list(solutions_Simplex_starmap)

    # choose the solution with the lowest RSS
    for solutions in solutions_Simplex, solutions_BFGS:
        for solution in solutions:
            if solution.fun < obj:
                obj = solution.fun
                solution_save = solution
    solution = solution_save
    solution = solution.x.reshape(len(solution.x),1)
    Parameter_opti = solution
    
    ## FIT - Apply optimization result to trainings set
    # COSMOSAC with optimized Paramters
    gamma_out_COSMO_opti_fit = list(map(run_COSMOSAC, repeat(Parameter_opti), Exp_to_COSMO, repeat(db)))
    gamma_out_COSMO_opti_fit = np.array(gamma_out_COSMO_opti_fit)

    assert list(Parameter_opti) == pytest.approx(list(np.array([4013.78,923.32, 3016.43])), abs=5)
    