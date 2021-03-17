# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:54:06 2021
@author: Luisa Peterson

This file summarizes all the functions to calculate the part of the
activity coefficient model that is constant regardless of the fitting-parameter.

Functions overview:
    UNIFAC_const_out = UNIFAC_const(Data, Experiments, Database, UNIFAC)
    UNISAC_const_out = UNIFAC_const(Data, Experiments, Database, UNISAC)
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

# %% UNIFAC CONSTANT
'''
This function calculates the part of UNIFAC that remains constant regardless of
the interaction parameters. The function does not need to be called by the opti-
mizer.
Input:
        - Data: constains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC database
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
Output:
        - UNIFAC_const_out: combinatorial part of UNIFAC and surface fraction of
        pure substances and substances in the mix (Theta_i, Theta)
'''
def UNIFAC_const(Data, Experiments, Database, UNIFAC):
    eps = np.finfo(float).eps
    # Unpack Input
        # Unpack Input
    if UNIFAC == False:
        structInfo = Data[2].copy()
    else:
        structInfo = Data[3].copy()  
    vdW_Properties = Database[0].copy()
    # Mole fraction of solute and solvent
    # Transpose?
    mole_fraction = np.array((Experiments[3],Experiments[4]))
    
    # length for loops
    k_length = len(structInfo.columns)
    m_length = len(vdW_Properties)
    
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
    return UNIFAC_const_out

# %% UNISAC CONSTANT
'''
This function calculates the part of UNISAC (= Stavermann Guggenheim Equation)
that remains constant regardless of
the segment area parameters. The function does not need to be called by the opti-
mizer.
Input:
        - Data: constains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC database
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
Output:
        - combinatorial: combinatorial part of UNIFAC and surface fraction of
'''
def UNISAC_const(Data, Experiments, Database, UNISAC):
    eps = np.finfo(float).eps
    # Unpack Input
        # Unpack Input
    if UNISAC == False:
        structInfo = Data[2].copy()
        structInfo = structInfo[0]
    else:
        structInfo = Data[3].copy()
        structInfo = structInfo[0]
    vdW_Properties = Database[0].copy()
    vdW_Properties = vdW_Properties.iloc[:,0:4]
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
    return combinatorial