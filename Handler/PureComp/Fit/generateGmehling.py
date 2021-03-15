# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 18:40:49 2021

@author: X202722
"""


import itertools
from itertools import repeat
from Adenventure.PhysModels.Thermodynamics.BinaryActivityCoefficients import CalcBinary_single, CalcBinary_single
import cCOSMO
import numpy as np
import pandas as pd
from Adenventure.PhysModels.KPIs.calc_KPI_gmehling import calc_KPI_gmehling

def generate_gmehling2(table, T,prod, name_prod1, components,db):
 # this function functionalizes the gmehling funciton which is descirbed in paper --> google gmeling solvent choice lle
    # 
    # T=293
    # define weighting in EVAL_SOL
    dif_vec = np.zeros(4)
    # capacity
    dif_vec[0]  =1 
    #distillation
    dif_vec[1] = 1
    #selectivity
    dif_vec[2] = 3
    # solvent loss
    dif_vec[3] =  0.75
    
    size = len(table)
    k = 0
    combinations = list(itertools.product(table['NAME'],repeat = 2))
    
    #calculate activity coeffients at infinite dilution

    for comps in components:    
        # print(comps)
        table['comb'+str(comps)] = list( CalcBinary_single([str(comps),names],db,T) for names in table['CAS#'])
    table['comb_prod'] = list( CalcBinary_single([str(int(prod)),names],db,T) for names in table['CAS#'])
    # calculate binary interaction (inf dil act coeff) coefficients for all components with defined solvents
    #get product index
    
    # now create distribution coefficents out of table these are defined as K = X_lm1/X_lm2 = gamma_LM2/gamma_LM1
    distribution_coeff = pd.DataFrame()
    distribution_coeff['solvent_combs'] = combinations
    for comps in components:    
        distribution_coeff['combination'+str(comps)+ ' '+ str(T)+'K'] = list( values[1][0]/values[0][0] for values in itertools.product(table['comb'+str(comps)],repeat = 2))
    distribution_coeff['comb_prod'+ ' '+ str(T)+'K'] = list( values[1][0]/values[0][0] for values in itertools.product(table['comb_prod'],repeat = 2))        
    
    # calculate necessary alphas
    ps = table['ps@25KhPa'].values
    names = table['NAME'].values
    ps_prod = ps[np.where(names == str(prod))]
    table['alpha1'] = table.apply(lambda x: x['comb_prod'][0]*ps_prod/x['ps@25KhPa'], axis = 1)
    
    # make list for distribution coefficients
    value_list = []   
    for subs,i in zip(distribution_coeff.columns, range(len(distribution_coeff.columns))):
        a = subs.find('combination')
        if a != -1:
            value_list.append(i)
    
    #calculate gmehling
    prod_pos = len(distribution_coeff.columns)-1
    for comps in itertools.product(value_list ,repeat = 2):
        # calc_KPI_gmehling(K_E_i, K_R_j, x_SR, alpha_dis, vec)
        # to do - implement creation of alphas according to solvent choice and gamma inf
        distribution_coeff[str(comps)+ ' '+ str(T)+'K'] = distribution_coeff.apply(lambda x: calc_KPI_gmehling(np.array(x.iloc[list(comps)]), x.iloc[prod_pos], 0.01, 10, dif_vec ), axis = 1)

       
    return distribution_coeff, table