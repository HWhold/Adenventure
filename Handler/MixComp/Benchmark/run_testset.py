# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:46:34 2021
@author: Luisa Peterson

This file calls the function that summerizes all operations on the testset
of the activity-coefficient-model handlers

Functions overview:
    out_testset = run_UNIFAC_testset(name_group, interaction_param_where)
    out_testset = run_UNISAC_testset(name_group, segment_area_param_where)
"""

# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
# For mapping function to repeat values
from itertools import repeat

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Handler.MixComp.Benchmark.set_up import load_Experiments, load_Data
from Adenventure.Handler.MixComp.Benchmark.find_parameter import filter_groups
from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import activity_coefficient_from_solubility
from Adenventure.PhysModels.MixComp.get_gamma_const import UNIFAC_const, UNISAC_const
from Adenventure.PhysModels.MixComp.get_gamma_fit import UNIFAC_fit, UNISAC_fit, run_COSMOSAC
from Adenventure.Handler.MixComp.Benchmark.evaluation import evaluation

# %% RUN UNIFAC TESTSET
'''
For the testset, this function runs all subfunctions in the right order
Input:
        - name group: name to load the right excel tables
        - interaction_param_where: information on where to find the interaction 
            parameters that need to be fit
Output:
        - out_testset:
            - Experiments_Error_extrapolation: Table of Experiments with activity coefficient
            - RSS_extrapolation: residual sum of squares for original UNIFAC and UNIFAC with own groups
            - Stats_extrapolation: Table of statistical evaluation (mean value, standard deviation)
'''
def run_UNIFAC_testset(name_group, Database_all):
    ### EXPERIMENTS FROM EXTRAPOLATION GROUP
    # Run UNIFAC fit with similar molecules and optimized interaction parameters
    # Load extrapolation Experiments
    model = 'UNIFAC'
    Experiments_all = load_Experiments(name_group,'extrapolation')
    Exp_all = Experiments_all.values
    Data_all = load_Data(name_group, model)
    
    Data, Database, Database_UNIFAC = filter_groups(Experiments_all, Data_all, Database_all)
    Exp_all = Experiments_all.values
    
    # experimentally obtained activity coefficients for the extrapolation group
    gamma_ext_exp = activity_coefficient_from_solubility(Data, Experiments_all)
    gamma_ext_exp = np.array(gamma_ext_exp)
    
    ## Run with own groups
    # load updated matrix of interaction parameters from database
    interaction_param = Database[2]
    # run constant part of UNIFAC
    UNIFAC_const_out = list(map(UNIFAC_const, repeat(Data), Exp_all, repeat(Database), repeat(False)))
    # run fit part of UNIFAC with optimized interaction parameters
    gamma_ext_opt = list(map(UNIFAC_fit,repeat('None'),repeat(Data),
                             Exp_all, repeat(Database),repeat(interaction_param),
                             repeat('None'),UNIFAC_const_out, repeat(False)))
    gamma_ext_opt = np.array(gamma_ext_opt)

    ## Run with own group split to UNIFAC groups
    # load matrix of interaction parameters
    interaction_param = Database_UNIFAC[2]
    # run constant part of UNIFAC
    UNIFAC_const_out_UNIFAC = list(map(UNIFAC_const, repeat(Data), Exp_all, 
                                       repeat(Database_UNIFAC), repeat(True)))
    # run firt part of UNIFAC, interaction parameters are given by UNIFAC
    gamma_ext_UNIFAC = list(map(UNIFAC_fit,repeat('None'),repeat(Data),
                                Exp_all, repeat(Database_UNIFAC),repeat(interaction_param),
                                repeat('None'),UNIFAC_const_out_UNIFAC, repeat(True)))
    gamma_ext_UNIFAC = np.array(gamma_ext_UNIFAC)
    
    # Statistic
    gamma_extrapolation = [gamma_ext_opt, gamma_ext_exp, 
                           gamma_ext_UNIFAC]
    Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation = evaluation(
    gamma_extrapolation, 'extrapolation', Experiments_all)


    out_testset = [Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation]
    return out_testset


# %% RUN UNISAC TESTSET
def run_UNISAC_testset(name_group, Database):
    '''
    For the testset, this function runs all subfunctions in the right order
    Input:
            - name group: name to load the right excel tables
            - segment_area_param_where: information on where to find the segment
              area parameters that need to be fit
    Output:
            - out_testset:
                - Experiments_Error_extrapolation: Table of Experiments with activity coefficient
                - RSS_extrapolation: residual sum of squares for original UNIFAC and UNIFAC with own groups
                - Stats_extrapolation: Table of statistical evaluation (mean value, standard deviation)
    '''
    ### EXPERIMENTS FROM EXTRAPOLATION GROUP
    # Run UNISAC fit with similar molecules and optimized interaction parameters
    # Load extrapolation Experiments
    model = 'UNISAC'
    Exp = load_Experiments(name_group,'extrapolation')
    Exp_ext = Exp.values
    Data = load_Data(name_group, model)
    
    ## experimentally obtained activity coefficients
    gamma_ext_exp = activity_coefficient_from_solubility(Data, Exp)
    gamma_ext_exp = np.array(gamma_ext_exp)
    
    ## Run with own groups
    # run constant part of UNISAC
    combinatorial = list(map(UNISAC_const, repeat(Data), Exp_ext, repeat(Database), repeat(False)))
    # fit UNISAC with optimized interaction paramters
    gamma_ext_opt = list(map(UNISAC_fit,repeat('No_fit_values'),repeat(Data),
                         Exp_ext, repeat(Database), combinatorial,
                         repeat('None'),repeat(False)))
    gamma_ext_opt = np.array(gamma_ext_opt)
    gamma_ext_opt = gamma_ext_opt[:,0] 

    # Run with own group split to UNIFAC groups
    # run combinatorial part of UNISAC
    combinatorial_UNISAC = list(map(UNISAC_const, repeat(Data), Exp_ext, repeat(Database), repeat(True)))
    # fit UNISAC, all interaction parameters given by UNISAC
    gamma_ext_UNISAC = list(map(UNISAC_fit,repeat('No_fit_values'),repeat(Data),
                         Exp_ext, repeat(Database), combinatorial_UNISAC,
                         repeat('None'),repeat(True)))
    gamma_ext_UNISAC = np.array(gamma_ext_UNISAC)
    gamma_ext_UNISAC = gamma_ext_UNISAC[:,0] 
    
    # evaluate the results - if wanted
    gamma_ext = [gamma_ext_opt, gamma_ext_exp, gamma_ext_UNISAC]
    Experiments_Error_ext, RSS_ext, Stats_ext = evaluation(
    gamma_ext, 'fit', Exp)

        
    out_testset = [Experiments_Error_ext, RSS_ext, Stats_ext]
    return out_testset


# %% RUN COSMOSAC TESTSET
def run_COSMOSAC_testset(name_group, db, Parameter_opti):
    '''
    For the testset, this function runs all subfunctions in the right order
    Input:
            - name group: name to load the right excel tables
            - segment_area_param_where: information on where to find the segment
              area parameters that need to be fit
    Output:
            - out_testset:
                - Experiments_Error_extrapolation: Table of Experiments with activity coefficient
                - RSS_extrapolation: residual sum of squares for original UNIFAC and UNIFAC with own groups
                - Stats_extrapolation: Table of statistical evaluation (mean value, standard deviation)
    '''
    # EXTRAPOLATION - Apply optimization result to test set
    # Load Data
    model = 'COSMOSAC'
    Data = load_Data(name_group,model)
    
    # load test set
    # load experiments
    Experiments_extrapolation = load_Experiments(name_group,'extrapolation')
    Exp_to_COSMO = Experiments_extrapolation.values
    # load active database
    names = Experiments_extrapolation['Solute (1)'].unique()  
    names = np.append(names, Experiments_extrapolation['Solvent (2)'].unique())    
    for iden in names:
        # Delaware
        db.add_profile(db.normalize_identifier(iden))    

    # COSMOSAC with optimized Paramters
    gamma_out_COSMO_opti_extrapolation = list(map(run_COSMOSAC, 
                                                  repeat(Parameter_opti), 
                                                  Exp_to_COSMO, repeat(db)))
    gamma_out_COSMO_opti_extrapolation = np.array(gamma_out_COSMO_opti_extrapolation)    
    # Experimental Results
    gamma_ext_exp = activity_coefficient_from_solubility(Data, Experiments_extrapolation)
    gamma_ext_exp = np.array(gamma_ext_exp)
    #COSMOSAC with original Parameters
    Parameter_COSMO = [4013.78, 923.32, 3016.43]
    gamma_out_COSMO_original_extrapolation = list(map(run_COSMOSAC, 
                                                      repeat(Parameter_COSMO), 
                                                      Exp_to_COSMO, repeat(db)))
    gamma_out_COSMO_original_extrapolation = np.array(gamma_out_COSMO_original_extrapolation)
    
    # Statistics
    gamma_extrapolation = [gamma_out_COSMO_opti_extrapolation, gamma_ext_exp, gamma_out_COSMO_original_extrapolation]
    Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation = evaluation(
    gamma_extrapolation, 'extrapolation', Experiments_extrapolation)
    
    out_testset = [Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation]
    return out_testset