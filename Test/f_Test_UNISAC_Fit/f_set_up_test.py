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
    path_empty = 'Data_test/06_Experiments_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    Experiments = pd.read_excel(path, sheet_name = fit_extrapolation, delimiter = ' ', header = 0)
    
    return Experiments

'''
In this function, solubility data of the molecules and their assignment to
structural UNIFAC groups are loaded.
Input:
        - name_group: group of solutes so that the right file can be loaded
Output:
        - Data: solubility data and structural information for solutes and solvents
''' 
def load_Data(name_group, model):
    # load information about the solutes
    # Columns: Material, FileName, Formula Smiles, CAS,	MPT °C,	HFUS J/mol
    Solutes = pd.read_excel('Data_test/06_Molecules.xlsx', sheet_name= name_group, index_col ="Material")
    
    # load information about the solvents
    # Columns: Material, Formula Smiles, CAS, Source MPT, MPT °C, Source HFUS, HFUS J/mol
    Solvents = pd.read_excel('Data_test/06_Molecules.xlsx', sheet_name= 'Solvents', index_col ="Material")
    
    if model == 'COSMOSAC':
        Data = [Solutes,Solvents]
        
    if model == 'UNIFAC':
        # load Structural Information
        # subgroup of each component and their counts for all species
        # structural information with own group
        path_empty = 'Data_test/06_StructuralInformation_'
        path_to_ctc2 = '-UNIFAC.xlsx'
        path = path_empty + str(name_group) + path_to_ctc2
        structural_information = pd.read_excel(path, delimiter = ' ', header = 0, index_col ="Name")
        # structural information where the own group is split to existing UNIFAC groups
        # used for comparison
        path_to_ctc2 = '-UNIFAC_original.xlsx'
        path = path_empty + str(name_group) + path_to_ctc2
        structural_information_UNIFAC = pd.read_excel(path, delimiter = ' ', header = 0, index_col ="Name")
        # combine output
        Data = [Solutes,Solvents,structural_information,structural_information_UNIFAC]  
    
    elif model == 'UNISAC':
        # load Structural Information
        # subgroup of each component and their counts for all species
        # structural information with own group
        path_empty = 'Data_test/06_StructuralInformation_'
        path_to_ctc2 = '-UNISAC.xlsx'
        path = path_empty + str(name_group) + path_to_ctc2
        structural_information_SG = pd.read_excel(path, sheet_name= 'StructInfo_SG', delimiter = ' ', header = 0, index_col ="Name")
        structural_information_seg = pd.read_excel(path, sheet_name= 'StructInfo_Seg', delimiter = ' ', header = 0, index_col ="Name")
        structural_information = [structural_information_SG, structural_information_seg]
        # structural information where the own group is split to existing UNISAC groups
        # used for comparison
        path_to_ctc2 = '-UNISAC_original.xlsx'
        path = path_empty + str(name_group) + path_to_ctc2
        structural_information_SG = pd.read_excel(path, sheet_name= 'StructInfo_SG', delimiter = ' ', header = 0, index_col ="Name")
        structural_information_seg = pd.read_excel(path, sheet_name= 'StructInfo_Seg', delimiter = ' ', header = 0, index_col ="Name")
        structural_information_UNISAC = [structural_information_SG, structural_information_seg]
        # combine output
        Data = [Solutes,Solvents,structural_information,structural_information_UNISAC]  

    return Data

'''
In this function, the UNIFAC database including the vdW properties for the groups,
the translation table from main groups to groups and the Group interaction 
parameters between the main groups is loaded.
Input:
        - name_group: group of solutes so that the right file can be loaded
Output:
        - Database: vdW Parameter, group translation table and interaction parameter
'''
def load_Database(name_group, model):
    # source: http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
    # Access data from database
    database = {}
    path_empty = 'Data_test/06_database_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(model) + path_to_ctc2
    
    excel_file = pd.ExcelFile(path)
    # load all sheets in excel file
    for sheet_name in excel_file.sheet_names:
        database[sheet_name] = excel_file.parse(sheet_name)
        
    if model == 'UNIFAC':
        # Subset of van der Waals Properties
        # vdW Properties are reindexed again
        vdW_Properties = database.get("UNIFAC_vanDerWaalsProperties")
        #vdW_Properties.set_index("Group_ID", inplace=True)
        
        # Translates Subgroups to Maingroups
        # Subset ist used for groups and main groups and numeration is changed
        group_to_mainGroup = database.get("group_to_mainGroup")
        group_to_mainGroup.set_index("Substructure", inplace=True)
        group_to_mainGroup.sort_values(by=['ID'], inplace=True) # sort by unique subgroup ID
        
        # Subset of interaction parameter
        # This values are fitted but used for the comparison of the fitted values
        interaction_param = database.get("UNIFAC_GroupInteractionParam")
        interaction_param.set_index("a_nm", inplace=True)
        # gamma at infinite dilution
        gamma_inf = database.get("gamma_inf")   
        gamma_inf.set_index("y_inf (373.15K)", inplace=True)
        
        # combine output
        Database = [vdW_Properties, group_to_mainGroup, interaction_param, gamma_inf]
        
    elif model == 'UNISAC':
         # van der Waals Properties and segment Area parameters
        vdW_Properties = database.get("vdW_Properties")
        UNISAC_Properties = database.get("UNISAC_Properties")
        # interaction parameter between different segments
        interaction_param = database.get("GroupInteractionParam")
        interaction_param.set_index("a_nm", inplace=True)
        # combine output
        Database = [vdW_Properties, UNISAC_Properties, interaction_param]
        
    return Database