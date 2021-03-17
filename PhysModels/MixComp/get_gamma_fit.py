# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:57:33 2021
@author: Luisa Peterson

This file contains the functions to calculate the part of the activity-
coefficient-model that depends on the fit-parameter

Functions overview:
    gamma = UNIFAC_fit(interaction_param_fit, Data, Experiments, Database, 
                       interaction_param, interaction_param_where, 
                       UNIFAC_const_out, UNIFAC)
    gamma = UNISAC_fit(segment_area_param_fit, Data, Experiments, Database, 
               combinatorial, segment_area_param_where, UNISAC)
    gamma = run_COSMOSAC(Parameter, Experiments, db)
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
import itertools
# supress warning
import warnings
# import COSMO
import cCOSMO

# %% UNIFAC FIT
'''
This function calculates the second part of UNIFAC. This part is dependent from
the Interaction Parameter and, thus, needs to be considered by the objective func-
tion.
Input:
        - interaction_param_fit: interaction parameter that are not given and, thus,
        to be fitted
        - Data: contains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC database
        - interaction_param: matrix of interaction param including Nones
        - interaction_param_where: localizes where the interaction params that 
        are missing are
        - UNIFAC_const_out: part of UNIFAC that is not depentend from interaction
        parameters
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
Output:
        - gamma: Activity coefficients of all species in the mixture, [-]
'''
def UNIFAC_fit(interaction_param_fit, Data, Experiments, Database, 
               interaction_param, interaction_param_where, UNIFAC_const_out, UNIFAC):
    
    # small number to add to divisions
    eps = np.finfo(float).eps
    # Allocate interaction param
    warnings.simplefilter(action='ignore', category=FutureWarning)
    if interaction_param_fit != 'None':
        interaction_param[interaction_param_where[0]] = interaction_param_fit.reshape(-1,1)
    # Reshape Interaction Param
    dim = int(np.sqrt(interaction_param.size))
    interaction_param = np.reshape(interaction_param,(dim, dim))
    interaction_param = pd.DataFrame(interaction_param)   
    
    # Unpack Input
    if UNIFAC == False:
        structInfo = Data[2].copy()
    else:
        structInfo = Data[3].copy()    
    Experiments = Experiments.copy()
    vdW_Properties = Database[0].copy()
    group_to_mainGroup = Database[1].copy() 
    
    # get Temperature    
    temperature = Experiments[0]
    
    # length for loops
    k_length = len((structInfo.columns))
    m_length = len(group_to_mainGroup)
    n_length = len(group_to_mainGroup)
        
    # Solute and Solvent for experiment e
    solute_name = Experiments[1]
    solvent_name = Experiments[2]
    # Structural information for experiment e
    structural_information = structInfo.loc[[solute_name,solvent_name],:]
    # UNIFAC const output for experiment e
    combinatorial = UNIFAC_const_out[0]
    Theta = UNIFAC_const_out[1]
    Theta_i = UNIFAC_const_out[2]
    Theta_i = Theta_i.squeeze()
    
    ## MIX
    # temperature depsendence of the activity coefficients
    interaction_coeff_temp = np.zeros((m_length, n_length))
    interaction_coeff_temp = list(np.exp(-interaction_param.iloc[m,n]/(temperature+eps)) for m,n in itertools.product(group_to_mainGroup.iloc[:,1],repeat = 2))        
    interaction_coeff_temp = np.reshape(interaction_coeff_temp,(m_length,n_length))   
                  
    # Now all data are avalable to calculate the group activity coefficients in
    # the binary system
    #group_activity_coeff = np.zeros((len(structural_information.columns)))
    #sum1 = np.tensordot(Theta,interaction_coeff_temp,[0,0])
    # sum2 = np.multiply(Theta_i,interaction_coeff_temp)
    # sum3 = np.sum(np.divide(sum2, sum1),axis=1)
    # vdW_Prop = np.array(vdW_Properties.iloc[:,1])
    # group_activity_coeff = vdW_Prop*(1-np.log(sum1)-sum3)
    
    sum1 = np.tensordot(Theta,interaction_coeff_temp,[0,0])
    sum2 = np.zeros(k_length)
    # Loop over all components
    for k in range(k_length):
            for m in range(m_length):
                sum2[k] += [(Theta[m]*interaction_coeff_temp[k,m])
                            /(sum1[m]+eps)]           
    
    group_activity_coeff = np.zeros((len(structural_information.columns)))
    # Group activity coefficients
    vdW_Prop = np.array(vdW_Properties.iloc[:,1])
    group_activity_coeff = vdW_Prop*(1-np.log(sum1+eps)-sum2) 
            
    # group activity coefficient
    sum1_i = np.tensordot(Theta_i,interaction_coeff_temp,[0,0])
    # Loop over all components       
    sum2_i = np.zeros(len(structural_information.columns))
    for k in range(k_length):
        for m in range (m_length):
            sum2_i[k] += [(Theta_i.iloc[m]*interaction_coeff_temp[k,m])/(sum1_i[m]+eps)]
           
    group_activity_coeff_i = 0
    group_activity_coeff_i = vdW_Properties.iloc[:,1]*(1-np.log(sum1_i+eps)-sum2_i)
            
    # Herewith all data are available to calculate the residual part of the
    # activity coefficients following the solutions of group concepst and finally
    # to calculate the required activity coefficient
    # residual part
    residual = 0
    residual += np.dot(structural_information.iloc[0,:],(group_activity_coeff-group_activity_coeff_i))
    ### ACTIVITY COEFFICIENT
    # activity coefficient
    gamma = np.exp(combinatorial+residual)
    return gamma

# %% UNISAC FIT
'''
This function calculates the second part of UNISAC. This part is dependent from
the segment area parameters and, thus, needs to be considered by the 
objective function.
Input:
        - segment_area_param_fit: interaction parameter that are not given and, thus,
        to be fitted
        - Data: contains structural information
        - Experiments: experiments to be considered
        - Database: UNISAC database
        - combinatorial: combinatorial part of UNISAC
        - segment_area_param_where: localizes where the segment area parameter
        that are missing are
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
Output:
        - gamma: Activity coefficients of all species in the mixture, [-]
'''

def UNISAC_fit(segment_area_param_fit, Data, Experiments, Database, 
               combinatorial, segment_area_param_where, UNISAC):
    
    # small number to add to divisions
    eps = np.finfo(float).eps
    # length for loops
    # number of segments
    n_seg = 7
    # number of molecules in the mixture
    n_length = 2    
    
    # get already exsisting segment area params
    Database_groups = Database[1].copy()
    segment_area_param = Database_groups.iloc[:,2:9]
    # Allocate segment area parameter
    warnings.simplefilter(action='ignore', category=FutureWarning)
    if segment_area_param_fit != 'No_fit_values':
        segment_area_param_missing = segment_area_param_where[0].flatten()
        segment_area_param.iloc[int(segment_area_param_missing),:] = segment_area_param_fit.reshape(1,-1)
        segment_area_param = np.array(segment_area_param, dtype=np.float)
    else:
        segment_area_param = np.array(segment_area_param, dtype=np.float)

    # Unpack Input
    if UNISAC == False:
        structInfo = Data[2].copy()
        structInfo = structInfo[1]
    else:
        structInfo = Data[3].copy()    
        structInfo = structInfo[1]
    # get Temperature    
    Experiments = Experiments.copy()
    temperature = Experiments[0]
    # get Interaction Param
    interaction_param = Database[2]
    # get Mole Fractions
    mole_fraction = np.array((Experiments[3],Experiments[4]))
    
    # Solute and Solvent for experiment e
    solute_name = Experiments[1]
    solvent_name = Experiments[2]
    # Structural information for experiment e
    structural_information = structInfo.loc[[solute_name,solvent_name],:]
    structural_information = np.array(structural_information)
    
    # total segment area of segment k in component i
    # Equation 3.18    
    v_s = np.zeros((n_seg,n_length))
    for i in range(n_length):
        for k in range(n_seg):
            v_s[k,i] =  np.dot(structural_information.T[:,i],segment_area_param.T[k,:])
    
    # Pure Component: segment area fraction of segment m in pure component i
    # Equation 3.20
    Theta = v_s/(sum(v_s)+eps)
    # Mixture: segment area fraction of segment m in the mixture  
    # Equation 3.22
    Theta_i = np.sum(v_s*mole_fraction,axis = 1)/(np.sum(np.sum(v_s*mole_fraction,axis = 1),axis = 0)+eps)    
    
    # temperature dependen segment interaction between u and v
    # interaction_coeff_temp = gabel
    interaction_coeff_temp = np.exp(-interaction_param/temperature+eps)
    
    # nuatrual logarithm of the segment activity coefficient for the mixture
    # Equation 3.21
    ln_gamma_mix = np.zeros(n_seg) 
    first_p = np.tensordot(Theta_i,interaction_coeff_temp,[0,0])
    second_p = np.multiply(Theta_i,interaction_coeff_temp)
    third_p=np.sum(np.divide(second_p,(first_p+eps)),axis=1)
    first_p[np.where(first_p <0)]=0
    ln_gamma_mix = 1-np.log(first_p+eps)-third_p
    
    # for i in range(n_length):
    # natural logarithm of the segment activity coefficient of segment k in
    # the pure component i
    # Equation 3.19
    ln_gamma_pure = np.zeros((n_length,n_seg)).T    
    for i in range(n_length):
            first_p = np.tensordot(Theta[:,i],interaction_coeff_temp,[0,0])
            second_p = np.multiply(Theta[:,i],interaction_coeff_temp)
            third_p=np.sum(np.divide(second_p,(first_p+eps)),axis=1)
            first_p[np.where(first_p <0)]=0
            ln_gamma_pure[:,i] = 1-np.log(first_p+eps)-third_p
        
    residual = np.tensordot(v_s,(np.array([ln_gamma_mix], ndmin=2).T-ln_gamma_pure),[0,0])
    residual[:,0]
    
    # Equation 1.49
    np.warnings.filterwarnings('ignore', 'overflow')
    ln_gamma = combinatorial + residual
    gamma = np.exp(ln_gamma)
    gamma = gamma[0]
    return gamma

#%% RUN COSMO SAC
def run_COSMOSAC(Parameter, Experiments, db):
    '''
    Run COSMO-SAC C++ source code (from Ian H. Bell, Erik Mickoleit, 
    Chieh-Ming Hsieh, Shiang-Tai Lin, Jadran Vrabec, Cornelia Breitkopf, 
    and Andreas Jäger, Journal of Chemical Theory and Computation 2020 16 (4),
    2635-2646, DOI: 10.1021/acs.jctc.9b01016) to get activity coefficients here.
    While calculating, the molecular interaction params are set and therefore 
    not always equal to the original COSMO-SAC value. The molecular interaction 
    params arealso called by the optimizer that optimizes this value by minimizing 
    the error between experimental and COSMO-SAC obtained interaction parameter.
    Input:
    - Parameter: values for the Parameter that needs to be optimized
    - Experiments: Experiments to be considered
    - db: active Delaware Database with relevant sigma profiles
    Output:
        - gamma: Activity Coefficient
    '''
    # load solute and solvent name of the required Experiment
    solute_name = Experiments[1]
    solvent_name = Experiments[2]
    names = [solute_name, solvent_name]
    
    #  temperature
    T = Experiments[0]
    # mole fraction
    x_1 = Experiments[5]
    z = np.array([x_1, 1-x_1])
    
    # COSMO modification by Hsieh etal.
    COSMO = cCOSMO.COSMO3(names, db)
    
    # change molecular interaction parameter
    comb_consts_COSMO = COSMO.get_mutable_COSMO_constants()
    if len(Parameter) == 3:
        comb_consts_COSMO.c_OH_OH = Parameter[0]
        comb_consts_COSMO.c_OT_OT = Parameter[1]
        comb_consts_COSMO.c_OH_OT = Parameter[2]
    if len(Parameter) == 2:
        comb_consts_COSMO.A_ES = Parameter[0]
        comb_consts_COSMO.B_ES = Parameter[1]

    try:
        # calculate activity coefficent via COSMO
        lngamma = COSMO.get_lngamma(T, z)
    except:
        lngamma = np.array([0,0])
    # activity coefficient
    gamma = np.exp(lngamma[0])
    return gamma