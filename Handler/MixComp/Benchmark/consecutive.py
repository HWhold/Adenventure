# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:38:12 2021
@author: Luisa Peterson

This file calls the function that organizes the consecutive fit of the UNIFAC
interaction parameters.

Functions overview:
    Database_all = consecutive(seed, Exp, Data_all, Database_all)
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
# Factorial Sampling
from pyDOE import lhs
# For mapping function to repeat values
from itertools import repeat, starmap
# multiprocessing for parallization of starmap
import multiprocessing as mp
from multiprocessing import Pool

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Handler.MixComp.Benchmark.find_parameter import filter_groups, disassemble
from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import activity_coefficient_from_solubility
from Adenventure.PhysModels.MixComp.get_gamma_const import UNIFAC_const
from Adenventure.Handler.MixComp.Benchmark.optimization import optimization

# %% CONSECUTIVE FITTING
'''
This function is responsible for making the groups fit one by one. The experiments 
are divided into fit groups for this purpose. For each fit group this function
is called. In this function the data necessary for the fit are loaded and not
yet certain interaction parameters are fitted. A variety of initial conditions
are used for the fit. The result, the determined interaction parameters, are
then added to the interaction matrix. The determined interaction parameters do
not need to be adjusted for further fits on the fit list.
Input:
    - seed: to check reproducibility at different set of starting values:
    - Exp: experiments according to the fit list
    - Data_all: Data for all solutes and solvents in the trainings set
    - Database_all: unfiltered database
Output:
    - Database_all: database with updated interaction parameters
'''
def consecutive(seed, Exp, Data_all, Database_all):
    model = 'UNIFAC'
    # load Data
    Data, Database, _ = filter_groups(Exp, Data_all, Database_all)
    # get matrix of interactin parameter
    interaction_param = Database[2]
    interaction_param_index = interaction_param.index
    interaction_param = interaction_param.values.flatten()
    # disassemble interaction parameters into already given and still unknown
    interaction_param_where = disassemble(interaction_param, model)
    
    # Calculate activity coefficients from solubilites
    gamma_fit_exp = activity_coefficient_from_solubility(Data, Exp)
    
    # Run constant UNIFAC with own groups
    Exp_UNIFAC_const = Exp.values
    UNIFAC_const_out = list(map(UNIFAC_const, repeat(Data), Exp_UNIFAC_const, repeat(Database), repeat(False)))
    
    # Set up optimization problem
    # bounds on interaction parameters
    bnds_min = -2500
    bnds_max = 10000
    bnds_new = [(bnds_min, bnds_max)]
    bnds_new = pd.Series(bnds_new)
    bnds = []
    for i in range (int(len(interaction_param_where[0]))):
        bnds.append(bnds_new)
    bnds = pd.concat(bnds)
    bnds = tuple(bnds)
    
    # dimension of the optimization problem: 
    # "how many unknwon interaction parameters do I have?"
    dim = len(interaction_param_where[0])
    # number of sets of intial values
    num = dim**5
    
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
    interaction_param_sampling = SLOPE.T*LHS_design+OFFSET.T
    
    ## Optimization
    # initial value for the objective_UNIFAC function
    obj = 1e1000
    # vector of summarized inputs
    inputs = [Data_all, Exp, Database_all, UNIFAC_const_out, 
              gamma_fit_exp, interaction_param, interaction_param_where]
    
    # call BFGS optimizer
    # get parameter as inputs for starmap
    params = []
    for interaction_param_fit in interaction_param_sampling:
        interaction_param_fit = np.vstack(interaction_param_fit).astype(np.float)
        params.append(tuple((inputs, interaction_param_fit, 'BFGS', model)))
      
    # Run optimization with Simplex for the results of the BFGS optimizer
    # set num_workers equal to the number of CPU cores
    num_workers = 1#mp.cpu_count() 
    # create a pool of "num workers" subprocesses
    #with Pool(num_workers) as p:    
    solutions_BFGS_starmap = starmap(optimization, params)
    solutions_BFGS = list(solutions_BFGS_starmap)

    # go trough all solution of the optimizer
    for solution in solutions_BFGS:
        # accept solution if better than the best yet existing
        if solution.fun < obj:
            obj = solution.fun
            solution_save = solution
    
    # sort optimization results for best solutions
    interaction_param_fit_BFGS = []
    for solution in solutions_BFGS:
        interaction_param_fit_BFGS.append([solution.x, solution.fun])
    interaction_param_fit_BFGS = np.concatenate(interaction_param_fit_BFGS)
    interaction_param_fit_BFGS = np.reshape(interaction_param_fit_BFGS, (-1, 2))
    interaction_param_fit_BFGS = interaction_param_fit_BFGS[interaction_param_fit_BFGS[:,1].argsort()]
    interaction_param_fit_BFGS = np.delete(interaction_param_fit_BFGS,1,1)
    # select best 10% of the solutions obtained by the BFGS solver (ceiled to the next integer)
    num_simplex = np.ceil(0.1*num)
    interaction_param_fit_BFGS = interaction_param_fit_BFGS[0:int(num_simplex),:]

    # get parameter as inputs for starmap
    params = []
    for interaction_param_fit in interaction_param_fit_BFGS:
        interaction_param_fit = np.vstack(interaction_param_fit).astype(np.float)
        params.append(tuple((inputs, interaction_param_fit,'Simplex', model)))
      
    # Run optimization with Simplex for the results of the BFGS optimizer
    # set num_workers equal to the number of CPU cores
    num_workers = mp.cpu_count() 
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
            print("new Simplex solution accepted")
        
    #save best fitting interaction parameters
    solution = solution_save
    
    # save optimized interaction parameters into matrix of interactions
    # resize and annotate optimized interaction parameters
    solution_x = solution.x.reshape(len(solution.x),1)
    interaction_param[interaction_param_where[0]] = solution_x
    interaction_param.resize(int(np.sqrt(len(interaction_param))),int(np.sqrt(len(interaction_param))))
    interaction_param_update = pd.DataFrame(interaction_param)
    interaction_param_update = interaction_param_update.set_axis(interaction_param_index, axis='index')
    interaction_param_update = interaction_param_update.set_axis(interaction_param_index, axis='columns')
    # replace 'none' interactions in matrix of interation
    # result: updated interacton matrix with optimized interaction parameter
    for a in interaction_param_index:
        for b in interaction_param_index:
            Database_all[2].loc[a,b] = interaction_param_update.loc[a,b]
    
    return Database_all


