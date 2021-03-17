# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:42:53 2021
This file contains the stage construction for crystallization
@author: X202722
"""

import pandas as pd
import numpy as np
from Adenventure.PhysModels.Thermodynamics.Solubility_rev1 import calc_solubiltiy_with_gamma, calc_solubility
import cCOSMO

def UO_crystall_cooling(data,table,components2,Ts,db,Z,K):
    # data: VU-dictionary as dataframe [pd.DataFrame()]
    # table: database for solvent search, containing ID, pS,   [pd.DataFrame()]
    # components are the components used in this extraction evaluation: list[int]
    # Ts temperatures in K to evaluate given by solubility: [list]
    # db delaware Database containing sigma profiles [object]
    # Z amount of purification to be done Z = x_out/x_in [float]
    # Distribution coefficient y = K*x for crystallizsation [float]
    R_sol = []
    cristallization_save = pd.DataFrame()
    # calculate solubility for fit molecueles
    # first without gamma
    all_data = pd.DataFrame()
    
    for t in Ts:
        data['R_soll'+str(t)]=data.apply(lambda x: calc_solubility(t,abs(x['HFUS_FIT']), x['MPT_FIT']), axis= 1)
        # calc_solubiltiy_with_gamma(prods,db, Ts,R_sol)
        # get nammes of columns already produced
        columns_names = cristallization_save.columns
        # calculate new cristallization_save
        cristallization_save2=data.apply(lambda x: table['CAS#'].apply(lambda y: calc_solubiltiy_with_gamma([str(x['Filename']),y],db,t, x['R_soll'+str(t)])), axis= 1).T
        cristallization_save2.columns = data['Number of molecule'].values
        # concat the results
        cristallization_save = pd.concat([cristallization_save, cristallization_save2], ignore_index =True, axis =1)
        # calculate stage number
        temp = data.apply(lambda x: cristallization_save2.apply(lambda y: crystallization(x['Mol_perc'],Z,K, y[x['Number of molecule']]),axis=1), axis= 1).T
        # concat the results
        cristallization_save = pd.concat([cristallization_save, temp], ignore_index =True, axis =1)
        # rename columns
        cristallization_save.columns = list(columns_names) +list(str(t)+'K ' + str(x ) for x in components2)+list(str(t)+'K ' + 'StageNumber' + str(x ) for x in components2)
        data['Krist stage Number'] = cristallization_save.max()
    return data, cristallization_save

def crystallization(x_Vu,Z,K,x_sle):
    #this function calculates stage number according to Hannes
    part0 = np.log((Z*(x_Vu-1))/(Z*x_Vu-1 ))
    part1 = np.log(K*x_sle/(1-x_sle))
    return part0/part1