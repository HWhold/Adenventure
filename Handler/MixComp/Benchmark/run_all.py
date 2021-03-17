# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:46:34 2021
@author: Luisa Peterson

This file 

Functions overview:
    interaction_param, Experiments_with_Error, RSS, Stats = run_all_UNIFAC(seed, name_group)
    segment_area_param, Experiments_with_Error, RSS, Stats = run_all_UNISAC(seed, name_group)
"""

# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# builds and contains data structures and operators for accessing numeric tables
import pandas as pd

#%% IMPORT FUNCTION
'''
Here, functions are loaded from other Python files
'''
from Adenventure.Handler.MixComp.Benchmark.run_trainingsset import run_UNIFAC_trainingsset, run_UNISAC_trainingsset, run_COSMOSAC_chb_trainingsset, run_COSMOSAC_electrostatic_trainingsset
from Adenventure.Handler.MixComp.Benchmark.run_testset import run_UNIFAC_testset, run_UNISAC_testset, run_COSMOSAC_testset

# %% RUN ALL UNIFAC
'''
This function runs all subfunctions in the right order
Input:
    - seed: check reproducibility at different set of starting values
    - name group: name to load the right excel tables
Output:
    - interaction_param: updated matrix of interaction params with optimized values
    - Experiments_with_Error: Table of Experiments with activity coefficient
    - RSS: residual sum of squares for original UNIFAC and UNIFAC with own groups
    - Stats: Table of statistical evaluation (mean value, standard deviation)
'''
def run_all_UNIFAC(seed, name_group):

    out_trainingsset = run_UNIFAC_trainingsset(seed, name_group)
    interaction_param, Experiments_Error_fit, RSS_fit, Stats_fit, Database = out_trainingsset
    out_testset = run_UNIFAC_testset(name_group, Database)
    Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation = out_testset
    
    # revise output data
    Experiments_with_Error = [Experiments_Error_fit,Experiments_Error_extrapolation]
    RSS = [RSS_fit,RSS_extrapolation]
    RSS = pd.DataFrame(RSS)
    RSS.columns = ['UNIFAC params', 'optimized params']
    RSS.index = ['Fit', 'Extrapolation']
    Stats = [Stats_fit,Stats_extrapolation]
    return interaction_param, Experiments_with_Error, RSS, Stats

# %% RUN ALL UNISAC
'''
This function runs all subfunctions in the right order
inputs:
    - seed: check reproducibility at different set of starting values
    - name group: name to load the right excel tables
outputs:
    - segment_area_param: updated matrix of segment area params with optimized values
    - Experiments_with_Error: Table of Experiments with activity coefficient
    - RSS: residual sum of squares for original UNIFAC and UNIFAC with own groups
    - Stats: Table of statistical evaluation (mean value, standard deviation)
'''
def run_all_UNISAC(seed, name_group):
    
    out_trainingsset = run_UNISAC_trainingsset(seed, name_group)
    segment_area_param, Experiments_Error_fit, RSS_fit, Stats_fit, Database = out_trainingsset
    out_testset = run_UNISAC_testset(name_group, Database)
    Experiments_Error_ext, RSS_ext, Stats_ext = out_testset
    
    # revise output data
    Experiments_with_Error = [Experiments_Error_fit, Experiments_Error_ext]
    RSS = [RSS_fit,RSS_ext]
    RSS = pd.DataFrame(RSS)
    RSS.columns = ['UNISAC params', 'optimized params']
    RSS.index = ['Fit', 'Extrapolation']
    Stats = [Stats_fit,Stats_ext]
    return segment_area_param, Experiments_with_Error, RSS, Stats

#%% RUN ALL COSMOSAC chb
def run_all_COSMOSAC_chb(seed, name_group):
    '''
    This function runs all subfunctions in the right order
    inputs:
        - seed: check reproducibility at different set of starting values
        - name_group: name to load the right excel tables
    outputs:
    - optimizedParameter: optimized molecular interaction parameter
    - Experiments_with_Error: Table of Experiments with activity coefficient
    - RSS: residual sum of squares for original UNIFAC and UNIFAC with own groups
    - Stats: Table of statistical evaluation (mean value, standard deviation)
    '''

    out_trainingsset = run_COSMOSAC_chb_trainingsset(seed, name_group)
    Experiments_Error_fit, RSS_fit, Stats_fit, db, Parameter_opti = out_trainingsset
    out_testset = run_COSMOSAC_testset(name_group, db, Parameter_opti)
    Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation = out_testset
    
    # receive output data
    optimizedParameter = Parameter_opti
    Experiments_with_Error = [Experiments_Error_fit,Experiments_Error_extrapolation]
    RSS = [RSS_fit,RSS_extrapolation]
    RSS = pd.DataFrame(RSS)
    RSS.columns = ['COSMO_SAC params', 'optimized params']
    RSS.index = ['Fit', 'Extrapolation']
    Stats = [Stats_fit,Stats_extrapolation]
    return optimizedParameter, Experiments_with_Error, RSS, Stats

#%% RUN ALL ELECTROSTATIC
def run_all_COSMOSAC_electrostatic(seed, name_group):
    '''
    This function runs all subfunctions in the right order
    inputs:
        - seed: check reproducibility at different set of starting values
        - name_group: name to load the right excel tables
    outputs:
    - optimizedParameter: optimized molecular interaction parameter
    - Experiments_with_Error: Table of Experiments with activity coefficient
    - RSS: residual sum of squares for original UNIFAC and UNIFAC with own groups
    - Stats: Table of statistical evaluation (mean value, standard deviation)
    '''

    out_trainingsset = run_COSMOSAC_electrostatic_trainingsset(seed, name_group)
    Experiments_Error_fit, RSS_fit, Stats_fit, db, Parameter_opti = out_trainingsset
    out_testset = run_COSMOSAC_testset(name_group, db, Parameter_opti)
    Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation = out_testset
    
    # receive output data
    optimizedParameter = Parameter_opti
    Experiments_with_Error = [Experiments_Error_fit,Experiments_Error_extrapolation]
    RSS = [RSS_fit,RSS_extrapolation]
    RSS = pd.DataFrame(RSS)
    RSS.columns = ['COSMO_SAC params', 'optimized params']
    RSS.index = ['Fit', 'Extrapolation']
    Stats = [Stats_fit,Stats_extrapolation]
    return optimizedParameter, Experiments_with_Error, RSS, Stats