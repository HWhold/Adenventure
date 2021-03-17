# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 11:04:28 2021
@author: Luisa Peterson

This file calls the function that summerizes all operations on the trainingsset
of the activity-coefficient-model handlers

Functions overview:
    out_trainingsset = run_UNIFAC_testsetrun_UNIFAC_trainingsset(seed, name_group)
    out_trainingsset = run_UNIFAC_testsetrun_UNISAC_trainingsset(seed, name_group)
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
from itertools import repeat, starmap
# Factorial Sampling
from pyDOE import lhs
# multiprocessing for parallization of starmap
import multiprocessing as mp
from multiprocessing import Pool
# import COSMO
import cCOSMO

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Handler.MixComp.Benchmark.set_up import load_Experiments, load_OrderOfFit, load_Data, load_Database
from Adenventure.Handler.MixComp.Benchmark.find_parameter import filter_groups, disassemble
from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import activity_coefficient_from_solubility
from Adenventure.PhysModels.MixComp.get_gamma_const import UNIFAC_const, UNISAC_const
from Adenventure.PhysModels.MixComp.get_gamma_fit import UNIFAC_fit, UNISAC_fit, run_COSMOSAC
from Adenventure.Handler.MixComp.Benchmark.optimization import optimization
from Adenventure.Handler.MixComp.Benchmark.evaluation import evaluation
from Adenventure.Handler.MixComp.Benchmark.consecutive import consecutive

#%% UNIFAC TRAININGSSET
def run_UNIFAC_trainingsset(seed, name_group):
    '''
    This function runs all subfunctions in the right order
    Input:
            - seed: check reproducibility at different set of starting values
            - name group: name to load the right excel tables
    Output:
            - out_trainingsset:
                - Experiments_Error_fit: Table of Experiments with activity coefficient
                - RSS_fit: residual sum of squares for original UNIFAC and UNIFAC with own groups
                - Stats_fit: Table of statistical evaluation (mean value, standard deviation)
                - interaction_param_where: information on where to find the interaction 
                parameters that need to be fit
    '''
    # Load Data and Database
    # Load Experiments
    model = 'UNIFAC'
    Experiments_all = load_Experiments(name_group,'fit')
    Order_Fit = load_OrderOfFit(name_group)
    Exp_all = Experiments_all.values
    Data_all = load_Data(name_group, model)
    Database_all = load_Database(name_group, model)
    
    # Select Experiments from Fit Group k
    for k in range (np.max(Order_Fit.iloc[:,0])+1):
        print('Fitting Group', k+1, 'out of', np.max(Order_Fit.iloc[:,0])+1)
        Order = Order_Fit[(Order_Fit['Fit']==k)]
        Exp = []
        for j in range (len(Order)):
            Exp_j = Experiments_all[(Experiments_all['Solute (1)'] == Order.Molecule.iloc[j]) & (Experiments_all['Solvent (2)'] == Order.Solvent.iloc[j])]
            Exp.append(Exp_j)
        Exp = pd.concat(Exp)
        Exp.index = range(len(Exp))
        # run optimization and obtain the updated database that is the input
        # for the next consecutive fit
        Database_all = consecutive(seed, Exp, Data_all, Database_all)

    ### EXPERIMENTS FROM FIT GORUP
    # Run results for all Experiments in the Fit group
    Data, Database, Database_UNIFAC = filter_groups(Experiments_all, Data_all, Database_all)
    # Calculate activity coefficients from solubilites
    gamma_fit_exp = activity_coefficient_from_solubility(Data, Experiments_all)
    
    ## Run with own groups
    # load updated matrix of interaction parameters from database
    interaction_param = Database[2]
    interaction_param_where = disassemble(interaction_param,model)
    # run constant part of UNIFAC
    UNIFAC_const_out = list(map(UNIFAC_const, repeat(Data), Exp_all, repeat(Database), repeat(False)))
    # run fit part of UNIFAC with optimized parameter
    gamma_fit_opt = list(map(UNIFAC_fit,repeat('None'),repeat(Data),
                         Exp_all, repeat(Database),repeat(interaction_param),
                         repeat(interaction_param_where),UNIFAC_const_out, repeat(False)))
    gamma_fit_opt = np.array(gamma_fit_opt)

    ## Run with own group split to UNIFAC groups
    # load matrix of interaction parameters
    interaction_param = Database_UNIFAC[2]
    interaction_param_where = disassemble(interaction_param, model)
    # run constant UNIFAC
    UNIFAC_const_out_UNIFAC = list(map(UNIFAC_const, repeat(Data), Exp_all, repeat(Database_UNIFAC), repeat(True)))
    # fit UNIFAC, all interaction parameters given by UNIFAC
    gamma_fit_UNIFAC = list(map(UNIFAC_fit,repeat('None'),repeat(Data), 
                                Exp_all, repeat(Database_UNIFAC), repeat(interaction_param),
                                repeat(interaction_param_where),UNIFAC_const_out_UNIFAC, repeat(True)))
    gamma_fit_UNIFAC = np.array(gamma_fit_UNIFAC)
    
    # evaluate the results - if wanted
    gamma_fit = [gamma_fit_opt, gamma_fit_exp, gamma_fit_UNIFAC]
    Experiments_Error_fit, RSS_fit, Stats_fit = evaluation(
    gamma_fit, 'fit', Experiments_all)
        
    out_trainingsset = [interaction_param, Experiments_Error_fit, RSS_fit,
                        Stats_fit, Database_all]
    return out_trainingsset
   
#%% UNISAC TRAININGSSET
def run_UNISAC_trainingsset(seed, name_group):
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
                - segment_area_param_where: information on where to find the segment_area_parameter
                parameters that need to be fit
    '''
    # Load Data andExperiments
    model = 'UNISAC'
    Exp = load_Experiments(name_group,'fit')
    Exp_fit = Exp.values
    Data = load_Data(name_group, model)
    Database = load_Database(name_group, model)
    
    # disassemble segment area parameters that are unknown from those who are known
    UNISAC_Properties = Database[1]
    segment_area_param_where = disassemble(UNISAC_Properties, model)
    
    # Calculate activity coefficients from solubilites
    gamma_fit_exp = activity_coefficient_from_solubility(Data, Exp)
    
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
    num = dim#**3
    
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
    num_workers = mp.cpu_count() 
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

    # Run with UNISAC groups
    # run constant part of UNISAC
    combinatorial_UNISAC = list(map(UNISAC_const, repeat(Data), Exp_fit, repeat(Database), repeat(True)))
    # run fit part of UNISAC with UNISAC parameter
    gamma_fit_UNISAC = list(map(UNISAC_fit,repeat('No_fit_values'),repeat(Data),
                         Exp_fit, repeat(Database), combinatorial_UNISAC,
                         repeat(segment_area_param_where),repeat(True)))
    gamma_fit_UNISAC = np.array(gamma_fit_UNISAC)
    gamma_fit_UNISAC = np.array(gamma_fit_UNISAC)
    gamma_fit_UNISAC = gamma_fit_UNISAC[:,0] 
    
    # evaluate the results - if wanted
    gamma_fit = [gamma_fit_opt, gamma_fit_exp, gamma_fit_UNISAC]
    Experiments_Error_fit, RSS_fit, Stats_fit = evaluation(
    gamma_fit, 'fit', Exp)

        
    out_trainingsset = [segment_area_param, Experiments_Error_fit, RSS_fit,
                        Stats_fit, Database]  
    return out_trainingsset

#%% COSMOSAC chb TRAININGSSET
def run_COSMOSAC_chb_trainingsset(seed, name_group):
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
    model = 'COSMOSAC'
    Data = load_Data(name_group,model)
    # Load Experiments
    Experiments_fit = load_Experiments(name_group,'fit')
    Exp_to_COSMO = Experiments_fit.values
    
    # Calculate activity coefficients from experimentally determined solubilities
    gamma_fit_exp = activity_coefficient_from_solubility(Data, Experiments_fit)
    
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
    num = dim**3
    
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
    
    #COSMOSAC with original Parameters
    Parameter_COSMO = [4013.78, 923.32, 3016.43]
    gamma_out_COSMO_original_fit = list(map(run_COSMOSAC, repeat(Parameter_COSMO), Exp_to_COSMO, repeat(db)))
    gamma_out_COSMO_original_fit = np.array(gamma_out_COSMO_original_fit)

    # Statistics
    gamma_fit = [gamma_out_COSMO_opti_fit, gamma_fit_exp, gamma_out_COSMO_original_fit]
    Experiments_Error_fit, RSS_fit, Stats_fit = evaluation(
    gamma_fit, 'fit', Experiments_fit)
    
    out_trainingsset = [Experiments_Error_fit, RSS_fit, Stats_fit, db, Parameter_opti]  
    return out_trainingsset

#%% COSMOSAC ELECTROSTATIC TRAININGSSET
def run_COSMOSAC_electrostatic_trainingsset(seed, name_group):
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
                - Parameter_opti: optimized electrostatic parameters
    '''
    
    # Load Data
    model = 'COSMOSAC'
    Data = load_Data(name_group,model)
    # Load Experiments
    Experiments_fit = load_Experiments(name_group,'fit')
    Exp_to_COSMO = Experiments_fit.values
    
    # Calculate activity coefficients from experimentally determined solubilities
    gamma_fit_exp = activity_coefficient_from_solubility(Data, Experiments_fit)
    
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
    dim = 2
    # number of sets of initial values
    num = dim**3
    
    # bounds on interaction parameters (one magnitude within original value)
    # C_ES = A_ES+B_ES/T
    bnds_min = []
    bnds_max = []
    bnds = []
    bnds_min_AES = 652.569
    bnds_min.append(bnds_min_AES)
    bnds_max_AES = 65250.69
    bnds_max.append(bnds_max_AES)
    bnds_new = [(bnds_min_AES, bnds_max_AES)]
    bnds_new = pd.Series(bnds_new)
    bnds.append(bnds_new)
    
    bnds_min_BES = 1.4859*10**7
    bnds_min.append(bnds_min_BES)
    bnds_max_BES = 1.4859*10**9
    bnds_max.append(bnds_max_BES)
    bnds_new = [(bnds_min_BES, bnds_max_BES)]
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
    
    #COSMOSAC with original Parameters
    Parameter_COSMO = [6525.690000, 148590000.000000]
    gamma_out_COSMO_original_fit = list(map(run_COSMOSAC, repeat(Parameter_COSMO), Exp_to_COSMO, repeat(db)))
    gamma_out_COSMO_original_fit = np.array(gamma_out_COSMO_original_fit)

    # Statistics
    gamma_fit = [gamma_out_COSMO_opti_fit, gamma_fit_exp, gamma_out_COSMO_original_fit]
    Experiments_Error_fit, RSS_fit, Stats_fit = evaluation(
    gamma_fit, 'fit', Experiments_fit)
    
    out_trainingsset = [Experiments_Error_fit, RSS_fit,Stats_fit, db, Parameter_opti]  
    return out_trainingsset