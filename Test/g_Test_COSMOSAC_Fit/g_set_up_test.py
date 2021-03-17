# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:24:53 2021
@author: Luisa Peterson

This file summarizes all the functions for setting up the various 
activity-coefficient-model handlers.

Functions overview:
    Experiments = load_Experiments(name_group, fit_extrapolation)
    Order_Fit = load_OrderOfFit(name_group)
    Data = load_Data(name_group)
    Database = load_Database(name_group)
"""

# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# builds and contains data structures and operators for accessing numeric tables
import pandas as pd

'''
In this function, experiments are imported from Excel tables
Input:
        - name_group: group of solutes so that the right file can be loaded
    - fit_extrapolation: choose wheter the test or trainings set should be considered
Output:
        - Experiments: Experiments for the group and either fit or extrapolation
'''    
def load_Experiments(name_group, fit_extrapolation):
    # load experimental data for a given solute and different solvents
    # Columns: Temperature in °C, solvent name, solute name, solubility
    path_empty = 'Data_test/07_Experiments_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    Experiments = pd.read_excel(path, sheet_name = fit_extrapolation, delimiter = ' ', header = 0)
    
    return Experiments
