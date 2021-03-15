# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:40:12 2021
This function contains the wrapper for the offset calculation of
1) Vapor pressure
2) Discrete properties

@author: X202722
"""
import numpy as np
from Adenventure.PhysModels.PureComp.testMPT2 import make_parameter_MPT_fit
from Adenventure.PhysModels.PureComp.testHFUS import make_parameter_HFUS_fit
from Adenventure.PhysModels.PureComp.testBPT import make_parameter_BPT_fit
from Adenventure.Handler.PureComp.Fit.make_parameter_VAP_fit import make_parameter_VAP_fit
import pandas as pd
from Adenventure.Handler.PureComp.Fit.load_prereq2 import load_prereq2

def create_column_names(num_method,properti):
    #this function creates column names for panda array where create estimator fit stuff is saved so this monster array is kept readable
    names = []
    
    for counter in range(num_method):
        names.append(properti + ' Offsets' + str(counter))
        
    for counter in range(num_method):
        names.append(properti + ' Rel Dev result' + str(counter))
        
    names.append(properti + ' Best estimator')
    return names


# Test Functions
def create_estimator_fit(num_method, rest_val, Tb, na, M, properti):
    # this is the wrapper to create the fit of the experimental values with each chosen estimator
#   num_method
    # rest_val: simulated values from ddbst from vap1m
    # Tb: experimental value of evaluated property and molecule 
    # na: number of atoms heavier than hydrogen in evaluated molecule [integer]
    # M: molar weight of evaluated molecule in g/mol [float]
    # properti: Property to evaluate [string]
    
    #create result array
    results = np.ones(num_method*2+1)
    offsets=[]
    T_s = []
    for method in range(num_method):
        # get ddbst data
        if properti == 'BPT':
            T_sim = rest_val.iloc[method*3,4]
            # create offset
            if T_sim != 0:
                offset = make_parameter_BPT_fit(T_sim, Tb, method, na, M)
            else:
                offset = 0
            
        elif properti == 'MPT':
            T_sim = rest_val.iloc[method,5]
            # create offset
            if T_sim != 0:
                offset = make_parameter_MPT_fit(T_sim, Tb, method, na, M)
            else:
                offset = 0
        elif properti == 'HFUS':
            T_sim = rest_val.iloc[method,6]
            # create offset
            if T_sim != 0:
                offset = make_parameter_HFUS_fit(T_sim, Tb, method, na, M)
            else:
                offset = 0            
        
        #save offset
        offsets.append(offset)
        # save simulation temperature
        T_s.append(T_sim)
    # now load all offsets into db
    results[-2*num_method-1:-1*num_method-1]=offsets 
    
    # 3) Benchmark performance and choose best estimator
    # and results of fit
    results[-1*num_method-1:-1]=abs(np.array(T_s)-Tb)/(Tb) 
    
    # save best estimator 
    results[-1]=np.argmin(results[-num_method-1:-1])
    return results





def create_estimator_fit_VAP(num_method, VAP_val, VAP_exp, T_exp, Tb_exp, Tb_sim):
    # this is offset function for vapor pressure. This is a bit more complicated due the need of fitted Tb and starting Tb_nannoolal
     #create result array
     #output results is offsets of each method, then benchmark values and then the enumerator of the best method
    # results = np.ones(num_method*2+1)
    offsets=[]
    P_sim = []
    # load p and T from experiments (mbar and °C in Excel)
    p = pd.to_numeric(VAP_exp)
    p = p[~np.isnan(p)]
    #in K
    T = pd.to_numeric(T_exp)
    T = T[~np.isnan(T)]    
    
    #generate test sets with low pressures            
    T_sim = np.ones(14)
    
    #generate data for GCT reverse engineering
    for i in VAP_val.index:
        T_sim[i] = 345 + i*15
        
    # uses parameters from DDBST and creates VAP_values from experimental temperatures, then normalizes whith exp 
    # VAP data and plots results
    for method in range(num_method):
        # get ddbst data
        # from kPa to mbar
        p_sim = VAP_val.iloc[:,2+method]*10
        # create offset
        if any(VAP_val != 0):
            offset, p_sim_scal = make_parameter_VAP_fit (Tb_sim, Tb_exp,T_sim, p_sim, p.to_numpy().astype(float),T.to_numpy().astype(float), method)
        else:
            offset = 0  
        #save offset
        offsets.append(offset)
        # save simulation temperature
        P_sim.append(p_sim_scal)    
    #unload offset these brackets are a wrapper to keep list together
    results = pd.DataFrame([[offsets[0]]])
    results = results.append(pd.Series([offsets[1]]),ignore_index = True).copy()
    # save results here: first double parameter moller then nannoolal    
    # 3) Benchmark performance and choose best estimator before fit
    # and results of fit
    for p_sim,i in zip(P_sim, range(num_method)):
        results = results.append(pd.Series(np.mean(abs(np.array(p_sim_scal)-p)/(p) )).astype(object), ignore_index = True).copy()
    # save best estimator 
    results = results.append(pd.Series(np.argmin(results.iloc[-num_method:,0])), ignore_index = True).copy()
    return results    

def fit_wrapper_main (exp_data, ges, num_method_BPT, num_method_HFUS, num_method_MPT, num_method_VAP, path_ctc):
        #create sapce for fit molecules into exp_data2 BPT
        exp_data2 = exp_data.join(pd.DataFrame(np.zeros((len(exp_data.index),num_method_BPT*2+1)), columns=create_column_names(num_method_BPT,'BPT') ))
        #write down indexes where BPT starts and ends
        BPT_offset_start = len(exp_data.columns)
        BPT_offset_end = BPT_offset_start +num_method_BPT*2+1
        
        #create sapce for fit molecules into exp_data2 MPT
        exp_data2 = exp_data2.join(pd.DataFrame(np.zeros((len(exp_data.index),num_method_MPT*2+1)), columns=create_column_names(num_method_MPT,'MPT') ))
        #write down indexes where MPT starts and ends
        MPT_offset_start = BPT_offset_end
        MPT_offset_end = MPT_offset_start +num_method_MPT*2+1
        
        #create sapce for fit molecules into exp_data2 MPT
        exp_data2 = exp_data2.join(pd.DataFrame(np.zeros((len(exp_data.index),num_method_HFUS*2+1)), columns=create_column_names(num_method_HFUS,'HFUS') ))
        #write down indexes where HFUS starts and ends
        HFUS_offset_start = MPT_offset_end
        HFUS_offset_end = HFUS_offset_start +num_method_HFUS*2+1
        
        #create sapce for fit molecules into exp_data2 MPT
        exp_data2 = exp_data2.join(pd.DataFrame(np.zeros((len(exp_data.index),num_method_VAP*2+1)), columns=create_column_names(num_method_VAP,'VAP') ))
        #write down indexes where VAP starts and ends
        VAP_offset_start = HFUS_offset_end
        VAP_offset_end = VAP_offset_start +num_method_VAP*2+1
    
        #load all necessary data
        offsets = []
        # create DataFrame which can contain all offsets and add to exp_data2
        OffsetList = pd.DataFrame(columns=())
        for iii in exp_data2.index:
            material = iii   
            print(iii)
            name_material = ges['Filename'][material]
            para = 0
            MPT = ges.iloc[material,7]+273.15
            Tb = ges.iloc[material,3]+273.15
            na = ges.iloc[material,4]
            M = ges.iloc[material,5]
            HFUS = ges.iloc[material,8]
            VAP = ges.iloc[material,9:]
            _, vaps_data, rest_val =  load_prereq2 (name_material,path_ctc)
            
        
            #get offset for BPT
            exp_data2.iloc[material,BPT_offset_start :BPT_offset_end] = create_estimator_fit(num_method_BPT, rest_val, Tb, na, M, 'BPT').copy()
            #get offset for MPT
            exp_data2.iloc[material,MPT_offset_start :MPT_offset_end] = create_estimator_fit(num_method_MPT, rest_val, MPT, na, M, 'MPT').copy()
            #get offset for HFUS
            exp_data2.iloc[material,HFUS_offset_start :HFUS_offset_end] = create_estimator_fit(num_method_HFUS, rest_val, HFUS, na, M, 'HFUS').copy()
            #get offset  for VAP
            # first initiate necessary data for vap calc
            VAP_exp = VAP.iloc[:6]
            T_exp = VAP.iloc[6:]+273
            Tb_sim = rest_val.iloc[0,4]
            
            temp = create_estimator_fit_VAP(num_method_VAP, vaps_data, VAP_exp, T_exp, Tb, Tb_sim).copy() 
            # exp_data2 = exp_data2.iloc[:,52].astype(object).copy()
            exp_data2 = exp_data2.astype(object).copy()
            exp_data2.iloc[material,VAP_offset_start] = temp.iloc[0,0].copy()
            exp_data2.iloc[material,VAP_offset_start +1 :] = temp.iloc[1:,0].values.copy()
        return exp_data2 , [BPT_offset_start, BPT_offset_end, MPT_offset_start, MPT_offset_end, HFUS_offset_start, HFUS_offset_end, VAP_offset_start, VAP_offset_end]