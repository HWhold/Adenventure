# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 16:41:24 2021
@author: Luisa Peterson

Functions overview:
    gamma_UNIFAC = run_UNIFAC(seed, name_group)
    gamma_UNISAC = run_UNIFAC(seed, name_group)
    gamma_COSMOSAC = run_UNIFAC(seed, name_group)
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
# import COSMO
import cCOSMO

# %% UNIFAC
def run_UNIFAC(Data, Experiments, Database, UNIFAC_groups):
    '''
    This function runs UNIFAC. Alle equations are taken from:
    J. Gmehling u. a. Chemical Thermodynamics for Process Simulation. 2nd.
    Weinheim, Germany: Wiley-VCH Verlag und Co. KGaA, 2019. isbn:
    9783527809479.

    Input:
            - Data: constains structural information
            - Experiments: experiments to be considered
            - Database: UNIFAC database
            - UNIFAC_groups:
                - True: run only with original UNIFAC groups
                - False: run with own groups
    Output:
            - gamma: Activity coefficients of all species in the mixture, [-]
    '''
    eps = np.finfo(float).eps
    # Unpack Input
    
    # Whenever UNIFAC with our own group is used, a different set of structural
    # information needs to be loaded (UNIFAC_groups == FALSE) than when UNIFAC
    # is used with original UNIFAC groups (else)
    if UNIFAC_groups == False:
        structInfo = Data[2].copy()
    else:
        structInfo = Data[3].copy() 
    vdW_Properties = Database[0].copy()
    group_to_mainGroup = Database[1].copy() 
    interaction_param = Database[2].copy()
    # Mole fraction of solute and solvent
    mole_fraction = np.array((Experiments[3],Experiments[4]))
    
    # length for loops
    k_length = len(structInfo.columns)
    m_length = len(vdW_Properties)
    n_length = len(group_to_mainGroup)
    
    # get Temperature
    temperature = Experiments[0]
    # Solute name for experiment e
    solute_name = Experiments[1]
    # Solvent name for experiment e
    solvent_name = Experiments[2]
    structural_information = structInfo.loc[[solute_name,solvent_name],:]
    
    ## COMBINATORIAL PART
    # The van der Waals properties of the two compounds can be calculated with the
    # help of the van der Waals properties of the groups
    vdW_Properties_r = vdW_Properties.iloc[:,0]
    r = np.tensordot(structural_information,vdW_Properties_r,[1,0])
    vdW_Properties_q = vdW_Properties.iloc[:,1]
    q = np.tensordot(structural_information,vdW_Properties_q,[1,0]) 
            
    # Using these van der Waals properties for the given mole fractions, the
    # following iloc are obtained for the volume-mole-ratio (V) and the surface-
    # mole-ratio (F)
    # mol volume / mol surface * sum over all mole fractions
    V_sum = np.sum(r*np.array(mole_fraction))
    F_sum = np.sum(q*np.array(mole_fraction))
    # Volume-mole-ratio
    V = r[0]/(V_sum+eps)
    # Surface-mole-ratio
    F = q[0]/(F_sum+eps)
        
    # With these data, the combinatorial part can be calculated
    # Loop over every component
    combinatorial = 1-V+np.log(V+eps)-5*q[0]*(1-(V/(F +eps))+np.log(V/(F+eps)+eps))
    
    # RESIDUAL PART
    # MIX
    # The following group mole fractions and sruface area fractions are obtained
    # for the considered binary system at given mole fractions
    
    X = np.zeros(m_length)
    # Loop over every component and every group
    X_num = np.tensordot(structural_information,mole_fraction,[0,0])
    for m in range(m_length):
        # mole fraction of group m in mix
        X = X_num/(sum(X_num)+eps)
    
    Theta_num = np.array(vdW_Properties_q*X)
    # surface fraction of group m in mix
    Theta = np.zeros(m_length)
    Theta = Theta_num/(sum(Theta_num)+eps)
    
    # PURE COMPOUNDS 
    # group mole fractions
    X_i = np.zeros((m_length))
    for k in range(k_length): # for every group
        # mole fraction of group m and pure component i
        X_i[k] = structural_information.iloc[0,k]/(sum(structural_information.iloc[0,:])+eps)

    # surface area fractions
    Theta_i_num = np.zeros(m_length)
    Theta_i = np.zeros(m_length)
    Theta_i_num += vdW_Properties.iloc[:,1]*X_i
    Theta_i_sum = sum(Theta_i_num.T)
    Theta_i = Theta_i_num/Theta_i_sum

    # get output
    combinatorial = pd.Series(combinatorial)
    Theta = pd.Series(Theta)
    Theta_i = pd.DataFrame(Theta_i)
    UNIFAC_const_out = []
    UNIFAC_const_out.append(combinatorial)
    UNIFAC_const_out.append(Theta)
    UNIFAC_const_out.append(Theta_i)
    
    # Allocate interaction param
    # Reshape Interaction Param
    dim = int(np.sqrt(interaction_param.size))
    interaction_param = np.reshape(interaction_param,(dim, dim))
    interaction_param = pd.DataFrame(interaction_param)   
        
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
    sum1_i = sum1_i.flatten()
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
    gamma_UNIFAC = np.exp(combinatorial+residual)
    return gamma_UNIFAC

# %% UNISAC
def run_UNISAC(Data, Experiments, Database, UNISAC_groups):
    '''
    Input:
        - Data: constains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC database
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
    '''
    eps = np.finfo(float).eps
    # number of segments
    n_seg = 7
    # number of molecules in the mixture
    n_length = 2    
    # Unpack Input
    
    vdW_Properties = Database[0].copy()
    vdW_Properties = vdW_Properties.iloc[:,0:4]
    # get Interaction Param
    interaction_param = Database[2]
    
    # Whenever UNISAC with our own group is used, a different set of structural
    # information needs to be loaded (UNISAC_groups == FALSE) than when UNISAC
    # is used with original UNISAC groups (else)
    if UNISAC_groups == False:
        structInfo = Data[2].copy()
        structInfo = structInfo[1]
    else:
        structInfo = Data[3].copy()    
        structInfo = structInfo[1]
    
    Experiments = Experiments.copy()
    # get Temperature    
    temperature = Experiments[0]
    # Mole fraction of solute and solvent
    # Transpose?
    mole_fraction = np.array((Experiments[3],Experiments[4]))
    
    # Solute name for experiment e
    solute_name = Experiments[1]
    # Solvent name for experiment e
    solvent_name = Experiments[2]
    structural_information = structInfo.loc[[solute_name,solvent_name],:]

    # The van der Waals properties of the two compounds can be calculated with the
    # help of the van der Waals properties of the groups
    vdW_Properties_r = vdW_Properties.loc[:,'Rk']
    r = np.tensordot(structural_information,vdW_Properties_r,[1,0])
    vdW_Properties_q = vdW_Properties.loc[:,'Qk']
    q = np.tensordot(structural_information,vdW_Properties_q,[1,0]) 
            
    # Using these van der Waals properties for the given mole fractions, the
    # following iloc are obtained for the volume-mole-ratio (V) and the surface-
    # mole-ratio (F)
    # mol volume / mol surface * sum over all mole fractions
    V_sum = np.sum(r*np.array(mole_fraction))
    F_sum = np.sum(q*np.array(mole_fraction))
    # Volume-mole-ratio
    V = r[0]/(V_sum+eps)
    # Surface-mole-ratio
    F = q[0]/(F_sum+eps)
        
    # With these data, the combinatorial part can be calculated
    # Loop over every component
    combinatorial = 1-V+np.log(V+eps)-5*q[0]*(1-(V/(F +eps))+np.log(V/(F+eps)+eps))
    
    # get already exsisting segment area params
    Database_groups = Database[1].copy()
    segment_area_param = Database_groups.iloc[:,2:9]
    # Allocate segment area parameter
    segment_area_param = np.array(segment_area_param, dtype=np.float)
    
    # total segment area of segment k in component i
    # Equation 3.18    
    v_s = np.zeros((n_seg,n_length))
    for i in range(n_length):
        for k in range(n_seg):
            v_s[k,i] =  np.dot(structural_information.T.iloc[:,i],segment_area_param.T[k,:])
    
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
    gamma_UNISAC = gamma[0,0]
    return gamma_UNISAC

#%% COSMOSAC

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