# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:35:16 2021
@author: LuisaPeterson

This file summarizes all the functions for considering only the relevant
UNIFAC-groups and finding the interacion parameters to be fit.

Functions overview:
    Data, Database, Database_UNIFAC = filter_groups(Experiments, Data_all, Database_all)
    out = disassemble(interaction_param)
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

#%% FILTER GROUPS
'''
In this function the previously loaded data is filtered for relevant information.
For the activity coefficient determination via UNIFAC, only a few group parameters
are needed for the molecule-solvent mixtures of the fit groups. Only the groups 
relevant to the fit are considered within the fit.
Input: 
        - Experiments: Experiments in consecutive fitting
        - Data_all: Data before filtering
        - Database_all: Database before filtering
Output:
        - Data: Data after filtering
        - Database: Database after filtering (to run with own groups)
        - Database_UNIFAC: Database after filtering (to run with only original UNIFAC groups)
'''
def filter_groups(Experiments, Data_all, Database_all):
    # Filter for Solutes in Experiments; dublicates are not considered
    Solutes = Data_all[0].copy()
    solutes_fit = Experiments['Solute (1)'].unique()
    Solutes = Solutes.loc[solutes_fit]
    
    # Filter for Solvents in Experiments; dublicates are not considered
    Solvents = Data_all[1].copy()
    solvents_fit = Experiments['Solvent (2)'].unique()
    Solvents = Solvents.loc[solvents_fit]
    
    # get structural information for all solutes and solvents in Experiments
    # considering own groups
    structural_information_all = Data_all[2].copy()
    structural_information_solutes = structural_information_all.loc[solutes_fit]
    structural_information_solvents = structural_information_all.loc[solvents_fit]
    structural_information = [structural_information_solutes,structural_information_solvents]
    structural_information = pd.concat(structural_information)
    # only consider the subgroups that are part of used solutes and solvents
    # get considered subgroups
    FilterSubgroups = structural_information.loc[:, (structural_information != 0).any(axis=0)]
    FilterSubgroups = np.array(list(FilterSubgroups))
    FilterSubgroups = np.subtract(FilterSubgroups,1.0)
    # apply subset of subgroups to structural information
    structural_information = structural_information.iloc[:,FilterSubgroups]
    
    # repeat  line 142-155 for when the own group is split to UNIFAC groups
    structural_information_UNIFAC = Data_all[3].copy()
    structural_information_solutes_UNIFAC = structural_information_UNIFAC.loc[solutes_fit]
    structural_information_solvents_UNIFAC = structural_information_UNIFAC.loc[solvents_fit]
    structural_information_UNIFAC = [structural_information_solutes_UNIFAC,structural_information_solvents_UNIFAC]
    structural_information_UNIFAC = pd.concat(structural_information_UNIFAC)
    FilterSubgroups_UNIFAC = structural_information_UNIFAC.loc[:, (structural_information_UNIFAC != 0).any(axis=0)]
    FilterSubgroups_UNIFAC = np.array(list(FilterSubgroups_UNIFAC))
    FilterSubgroups_UNIFAC = np.subtract(FilterSubgroups_UNIFAC,1.0)
    # apply subset of subgroups to structural information
    structural_information_UNIFAC = structural_information_UNIFAC.iloc[:,FilterSubgroups_UNIFAC]
    
    # combine output
    Data = [Solutes, Solvents, structural_information, structural_information_UNIFAC]
   
    ## Database for UNIFAC with own groups   
    # Filter for relevant van der waals properties; reindex Group IDs
    vdW_Properties = Database_all[0].copy()
    vdW_Properties = vdW_Properties.iloc[FilterSubgroups,:]
    Group_ID = vdW_Properties.Group_ID.unique()
    Group_ID_new = np.arange(len(Group_ID))
    vdW_Properties_Group_ID = vdW_Properties.loc[:,'Group_ID'].replace(to_replace = Group_ID, value = Group_ID_new)
    vdW_Properties.loc[:,'Group_ID'] = np.array(vdW_Properties_Group_ID)
    vdW_Properties.set_index(["Group_ID","Subgroup"], inplace=True)
    
    # filtered translation table: get the main group when having the group number
    group_to_mainGroup = Database_all[1].copy()
    group_to_mainGroup = group_to_mainGroup.iloc[FilterSubgroups,:]
    group_to_mainGroup.loc[:,'ID'] = np.arange(len(group_to_mainGroup))
    MainGroups = group_to_mainGroup.Main_Group_ID.unique()
    MainGroups_new = np.arange(len(MainGroups))
    group_to_mainGroup.loc[:,'Main_Group_ID'] = group_to_mainGroup['Main_Group_ID'].replace(MainGroups, MainGroups_new)
    
    # filtered interaction parameters between main groups
    interaction_param = Database_all[2].copy()
    interaction_param = interaction_param.iloc[Group_ID-1,Group_ID-1]
    
    # activity coefficient at infinite dilution
    gamma_inf = Database_all[3].copy()
    gamma_inf = gamma_inf.iloc[FilterSubgroups,FilterSubgroups]
    
    # combine output
    Database = [vdW_Properties, group_to_mainGroup, interaction_param, gamma_inf]
    
    ## Database for UNIFAC with existing UNIFAC groups
    # repeat lines 174-199
    vdW_Properties = Database_all[0].copy()
    vdW_Properties = vdW_Properties.iloc[FilterSubgroups_UNIFAC,:]
    Group_ID = vdW_Properties.Group_ID.unique()
    Group_ID_new = np.arange(len(Group_ID))
    vdW_Properties_Group_ID = vdW_Properties.loc[:,'Group_ID'].replace(to_replace = Group_ID, value = Group_ID_new)
    vdW_Properties.loc[:,'Group_ID'] = np.array(vdW_Properties_Group_ID)
    vdW_Properties.set_index(["Group_ID","Subgroup"], inplace=True)
    
    group_to_mainGroup = Database_all[1].copy()
    group_to_mainGroup = group_to_mainGroup.iloc[FilterSubgroups_UNIFAC,:]
    group_to_mainGroup.loc[:,'ID'] = np.arange(len(group_to_mainGroup))
    MainGroups = group_to_mainGroup.Main_Group_ID.unique()
    MainGroups_new = np.arange(len(MainGroups))
    group_to_mainGroup.loc[:,'Main_Group_ID'] = group_to_mainGroup['Main_Group_ID'].replace(MainGroups, MainGroups_new)
    
    interaction_param = Database_all[2].copy()
    interaction_param = interaction_param.iloc[Group_ID-1,Group_ID-1]
    
    # combine output
    Database_UNIFAC = [vdW_Properties, group_to_mainGroup, interaction_param]
    return Data, Database, Database_UNIFAC
# %% DISASSEMBLE
'''
This function dissembles the interaction parameter for the main groups into 
two parts: Firstly, interaction parameters given by UNIFAC and secondly, interaction
parameters not given by UNIFAC. Those can be interaction parameter that are not
available for the considered main groups or interaction parameters for additional
groups.
Input: 
        - param: params of the model including Nones
Output: 
        if model == UNIFAC
        - out: interaction parameters disassmbled to whether they
        are given or not
        if model == UNISAC
        - out: segment area parameter disassmbled to whether they
        are given or not
'''
def disassemble(param, model):
    if model == 'UNIFAC':
        # get the position of interaction parameter that are not given
        interaction_param_None_where = np.array(np.where(param == 'None')).T
        # get the position of interaction parameter that are given
        interaction_param_UNIFAC_where = np.array(np.where(param != 'None')).T
        out = [interaction_param_None_where, interaction_param_UNIFAC_where]
    if model == 'UNISAC':
        # get the position of interaction parameter that are not given
        segment_area_param_None_where = np.array(np.where(param['A'] == 'None')).T
        # get the position of interaction parameter that are given
        segment_area_param_UNISAC_where = np.array(np.where(param['A']  != 'None')).T
        out = [segment_area_param_None_where, segment_area_param_UNISAC_where]
    return out