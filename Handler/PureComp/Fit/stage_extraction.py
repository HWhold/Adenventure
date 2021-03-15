# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:45:16 2021
this contains stage construction for extraction via kremser equation
1) create new x_out for impurity under guessed yield of product
2) calculate new concentration considering dilution with solvent to max solubility @ temp for product
3) calc kremser for product -> stage number
4) calc product retention
5) calc solvent free composition
6) 
@author: X202722
"""

from Adenventure.ProcessModels.Extraction.kremser import calc_kremser, calc_comp_kremser
from Adenventure.Handler.PureComp.Fit.generateGmehling import generate_gmehling2
import pandas as pd
import itertools
import numpy as np 


def UO_extraction(data, table, process_data, components2, prod, db, solubilities,Ts,S,Z):
    # data: VU-dictionary as dataframe
    # Process data: fitted physical properties
    # components are the components used in this extraction evaluation: list[int]
    # prod is the product of this extraction evaluation it should stay in extract: int 
    # db is the object of the cosmo DB for activity coefficient calculations: db database object
    # solubilities are the product solubilities at given temperatures: dataframe with solubilities
    # Ts temperatures in K to evaluate given by solubility: [list]
    # S yield of product in process x_out/x_in
    # Z purification of impurities done in process x_out/x_int
    
    
    name_prod1 = components2
    prod2 = process_data['Number of molecule'][int(prod)]
    extraction_ress = []
    # calculate mol flow rate for outlet stream if maximum loss is defined by product. n_total_in  = 1 mol where S is max Loss Product and Z is Purification rate impurity
    
    n_1Prod =  data['Mol_perc'][prod]
    # this calculates residue x for guessed maximal loss of prodcut S and cut Z of VU
    data['n_out'] = data['Mol_perc'].apply(lambda x: (-S*Z*n_1Prod*x/((Z - 1)*x - n_1Prod)))
    for T,i in zip(Ts,range(len(Ts))):          
        # extraction_res.append(generate_gmehling(table, T,prod, name_prod1, components,db))
        distribution_coefficient, table = (generate_gmehling2(table, T,prod2, name_prod1, components2,db))
   
        column_names = distribution_coefficient.columns
        #spread solubilities to distribution coefficent length
        
        solubility_list = pd.DataFrame(list(itertools.product(solubilities.iloc[:,i].values,repeat = 2))).iloc[:,0]
        # calculate true mol perc of VU according to dilution of LM -> max-solubility of product @ T
        # solubilities_comp =  data['Mol_perc'].apply(lambda x : solubilities.iloc[:,i].apply(lambda value: (x*value)/(x*value+(1-x))))
        solubilities_comp =  data['Mol_perc'].apply(lambda x : solubilities.iloc[:,i].apply(lambda value: (x)/(1+(n_1Prod)/value)))
        solubilities_comp_list = (solubilities_comp.apply(lambda x: pd.DataFrame(list(itertools.product(x,repeat = 2))).iloc[:,0], axis = 1)).T
        # calculate the outgoing impurity fraction
        solubilities_comp2 =  data['n_out'].apply(lambda x : solubilities.iloc[:,i].apply(lambda value: (x*value)/(x*value+(1-x))))
        solubilities_comp2_list = (solubilities_comp2.apply(lambda x: pd.DataFrame(list(itertools.product(x,repeat = 2))).iloc[:,0], axis = 1)).T
        #name every column
        # concat everything
        distribution_coefficient = pd.concat([distribution_coefficient,solubilities_comp_list,solubilities_comp2_list],ignore_index = True, axis = 1) 
        
        distribution_coefficient.columns = list(column_names )+ list('x_in'+str(x ) for x in components2)+list('x_out'+str(x ) for x in components2)
        # calc_kremser(K, x_in, x_out, V_L) - necessary stages for purification
        stage_number = data['Number of molecule'].apply(lambda x: distribution_coefficient.apply(lambda K: calc_kremser(K['combination'+str(x)+ ' '+ str(T)+'K'],K['x_in'+str(x )],K['x_out'+str(x)],0.1), axis = 1)).T
        #sometimes stage construction becomes negative -> fail -> parsing to nan
        stage_number[stage_number<0] = np.nan
        # concat arrays to one big result array
        extraction_res = pd.concat([stage_number, distribution_coefficient],ignore_index = True, axis = 1).copy()
        # rename columns
        extraction_res.columns = [list('Kremser'+str(x ) for x in components2)+list(distribution_coefficient.columns)]
        
        #calc retained product in solvent
        #via  calc_comp_kremser(K, x_in, n, V_L)
        extraction_res2 = data['Number of molecule'].apply(lambda x: extraction_res.apply(lambda K: calc_comp_kremser(float(K['comb_prod'+ ' '+ str(T)+'K']),float(K['x_in'+str(int(prod2))]),float(K['Kremser'+str(x)]),0.1), axis = 1)).T
        
        # back to solvent free view (n_prod = x_prod2/x_prod_solubility/(1-x_prod2)) with n_prod_start = x_prod_start * 1 mol and n_ges = n_lm+x_prod*n_ges
        ResidueProd = pd.DataFrame(list((a/(1-a)*1/solubility_list for a in extraction_res2.T.values)))
        ResidueProd = ResidueProd.T
        
        # concat arrays to one big result array
        extraction_res = pd.concat([extraction_res, ResidueProd],ignore_index = True, axis = 1)
        
        #rename columns
        extraction_res.columns = [list('Kremser'+str(x ) for x in components2)+list(distribution_coefficient.columns)+list('ResidueProd'+str(x ) for x in components2)]
        
        # now rate options
        ResidueProd.columns = list('ResidueProd'+str(x ) for x in components2)
        idx = ResidueProd.idxmax()
        
        
        data['extr stage Number'+ str(T)] = pd.DataFrame(list(stage_number.iloc[idx.iloc[n], n] for n in range(len(components2))))
        data['extra Product extraction_residue'+str(T)] = ResidueProd.max().values
        extraction_ress.append(extraction_res)
    return data, ResidueProd, extraction_ress