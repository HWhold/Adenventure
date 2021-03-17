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
# multiprocessing for parallization of starmap
import multiprocessing as mp
from multiprocessing import Pool

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Test.f_Test_UNISAC_Fit.f_set_up_test import load_Experiments, load_Data, load_Database
from Adenventure.Handler.MixComp.Benchmark.find_parameter import disassemble
from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import load_exp_gamma
from Adenventure.PhysModels.MixComp.get_gamma_const import UNISAC_const
from Adenventure.PhysModels.MixComp.get_gamma_fit import UNISAC_fit
from Adenventure.Handler.MixComp.Benchmark.optimization import optimization

def test_UNISAC_Handler():
    # set up
    # Load Data andExperiments
    seed = np.random.randint(1e9)
    name_group = 'test'
    model = 'UNISAC'
    Exp = load_Experiments(name_group,'fit')
    Exp_fit = Exp.values
    Data = load_Data(name_group, model)
    Database = load_Database(name_group, model)
    
    # disassemble segment area parameters that are unknown from those who are known
    UNISAC_Properties = Database[1]
    segment_area_param_where = disassemble(UNISAC_Properties, model)
    
    # Calculate activity coefficients from solubilites
    gamma_fit_exp = load_exp_gamma(Exp)
    
    # Run constant UNISAC
    # own groups
    combinatorial = list(map(UNISAC_const, repeat(Data), Exp_fit, repeat(Database), repeat(False)))
    
    # Set up optimization problem
    # bounds on interaction parameters
    bnds_min = 0
    bnds_max = 2.5
    bnds_new = [(bnds_min, bnds_max)]
    bnds_new = pd.Series(bnds_new)
    bnds = []
    for i in range (int(len(segment_area_param_where[0])*7)):
        bnds.append(bnds_new)
    bnds = pd.concat(bnds)
    bnds = tuple(bnds)
    
    ## number of sets of intial values
    dim = len(segment_area_param_where[0])*7
    num = dim**2
    
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
    segment_area_param_sampling = SLOPE.T*LHS_design+OFFSET.T
    
    ## Optimization
    # initial value of obj funciton
    obj = 1e1000
    
    # vector of summarized inputs
    inputs = [Data, Exp, Database, combinatorial,
              segment_area_param_where, gamma_fit_exp]
    
    params = []
    for segment_area_param_fit in segment_area_param_sampling:
        segment_area_param_fit = np.vstack(segment_area_param_fit).astype(np.float)
        params.append(tuple((inputs, segment_area_param_fit, 'BFGS', model)))
    # Run optimization with Simplex for the results of the BFGS optimizer
    # set num_workers equal to the number of CPU cores
    #num_workers = mp.cpu_count() 
    # create a pool of "num workers" subprocesses
    #with Pool(num_workers) as p:    
    solutions_BFGS_starmap = starmap(optimization, params)
    solutions_BFGS = list(solutions_BFGS_starmap)
    
    # call BFGS optimizer
    # solutions_BFGS = list(map(optimization,repeat(inputs), segment_area_param_sampling, 
    #                      repeat(bnds), repeat('BFGS'), model))

    # go trough all solution of the optimizer
    for solution in solutions_BFGS:
        # accept solution if better than the best yet existing
        if solution.fun < obj:
            obj = solution.fun
            solution_save = solution
    
    # sort optimization results for best solutions
    segment_area_param_fit_BFGS = []
    for solution in solutions_BFGS:
        segment_area_param_fit_BFGS.append([solution.x, solution.fun])
    segment_area_param_fit_BFGS = np.concatenate(segment_area_param_fit_BFGS)
    segment_area_param_fit_BFGS = np.reshape(segment_area_param_fit_BFGS, (-1, 2))
    segment_area_param_fit_BFGS = segment_area_param_fit_BFGS[segment_area_param_fit_BFGS[:,1].argsort()]
    segment_area_param_fit_BFGS = np.delete(segment_area_param_fit_BFGS,1,1)
   # select best 10% of the solutions obtained by the BFGS solver (ceiled to the next integer)
    num_simplex = np.ceil(0.1*num)
    segment_area_param_fit_BFGS = segment_area_param_fit_BFGS[0:int(num_simplex),:]

    # get parameter as inputs for starmap
    params = []
    for segment_area_param_fit in segment_area_param_fit_BFGS:
        segment_area_param_fit = np.vstack(segment_area_param_fit).astype(np.float)
        params.append(tuple((inputs, segment_area_param_fit, 'Simplex', model),))
        
    # Run optimization with Simplex for the results of the BFGS optimizer
    # set num_workers equal to the number of CPU cores
    #num_workers = mp.cpu_count() 
    # create a pool of "num workers" subprocesses
    #with Pool(num_workers) as p:    
    solutions_Simplex_starmap = starmap(optimization, params)
    solutions_Simplex = list(solutions_Simplex_starmap)
    
    # go through all solutions obtained by the simplex solver
    for solution in solutions_Simplex:
        # accept solution if better than the best yet existing
        if solution.fun < obj:
            obj = solution.fun
            solution_save = solution
        
    #save best fitting interaction parameters
    solution = solution_save
    solution_x = solution.x.reshape(len(solution.x),1)
    
    # save best solution and put it into the database
    Database_groups = Database[1].copy()
    segment_area_param = Database_groups.iloc[:,2:9]
    segment_area_param_missing = segment_area_param_where[0].flatten()
    segment_area_param.iloc[segment_area_param_missing,:] = solution_x.T
    Database_groups.iloc[:,2:9] = segment_area_param
    Database[1] = Database_groups

    ### ALL EXPERIMENTS FROM FIT GORUP
    # Run results for all Experiments in the Fit group
    # Run with own groups
    # run constant part of UNIFAC
    gamma_fit_opt = list(map(UNISAC_fit,repeat('No_fit_values'),repeat(Data),
                         Exp_fit, repeat(Database), combinatorial,
                         repeat(segment_area_param_where),repeat(False)))
        
    # run fit part of UNISAC with optimized parameter
    gamma_fit_opt = np.array(gamma_fit_opt)
    gamma_fit_opt = gamma_fit_opt[:,0] 

    assert list(solution_x) == pytest.approx(list(np.array([0.54,0,0,0,0,0,0])), abs=0.5)
    