# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 11:20:02 2020
@author: X202722

Only commented, no modifications in the code

Equations based on Moddley 2016
"""

# build UNISAC model after email rarey 
import numpy as np
import time
import pandas as pd
#%% Set UP
# mashine einteraction_coeff_templon
eps = np.finfo(float).eps
# binary System
n_length = 2
# temperature
temperature = 298.15
# mole fraction
mole_fraction = np.array([.0448, 1-.0448])
# frequency of group j in component i
# =^structural information; CH3, CH2, ACH, AC, CC3COO, COOH
v2 = np.array([[0,0,4,2,1,1],[1,1,0,0,1,0]])
# 49 = extended UNIFAC
Q_segs = np.zeros(49)
# group contribution values of modified UNIFAC
# aspirin in ethyl acetate
# equal to UNIFAC parameters
# CH3, CH2, ACH, AC, CC3COO, COOH
R = np.array([0.9011, 0.6744, 0.5313, 0.3652, 1.9031, 1.3013])
Q = np.array([0.8480, 0.54, 0.4, 0.12, 1.7280, 1.2240])
vdW_Properties = pd.DataFrame((R,Q))

# first load data here from rarey email segments and values 
# index gibt nicht 0 Werte an
ind = np.array([0,7,14,15,16,17,18,20,21,35,36,37,42,43,44])
val = [0.848,0.54,0.2012,0.0162,0.0198,0.007,0.058,0.0041,0.228,0.008,0.1985,0.6914,0.008,0.1985,0.6914]

# molecule segments use
# UNISAC groups
v = np.array([[1,0,4,1,1,1,1,1],[2,1,0,0,0,0,1,1]])
# number of segments -> EXTENDED  UNIFAC
n_seg = 7
# now put values into matrix form for segment matrix
Q_segs[ind] = np.array(val )
# squared Table with UNISAC parameters
Q_seqs_rs = np.reshape(Q_segs,(7,7))

#%% Not dependend on a_s
# COMBINATORIAL
## calc combinatorical part
# molecule vdW surface area
# Equation 3.15
r = np.sum(v2*R,axis = 1)
# Equation 3.13
V= r/sum(mole_fraction*r)

# Equation 3.16
q = np.sum(v2*Q,axis = 1)
# molecule vdW volume
# Equation 3.14
F = q/sum(mole_fraction*q)

# Equation 5.2
combinatorial = 1-V+np.log(V)-5*q*(1-V/F+np.log(V/F))

## RESIDUAL - NEWS START HERE
zeitanfang = time.time()
# calculate v_s (= Omega)
# total segment area of segment k in component i
# Equation 3.18
v_s = np.zeros(v.T.shape)
for k in range(7): # nicht segmentanzahl sondern anzahl an UNISAC gruppen
    for i in range(2):
            v_s[k,i] =  np.dot(v.T[:,i],Q_seqs_rs.T[k,:])
print(v_s)

# Pure Component
# segment area fraction of segment m in pure component i
# Equation 3.20
Theta = v_s/(sum(v_s)+eps)
        
# Mixture       
# segment area fraction of segment m in the mixture  
# Equation 3.22
Theta_i = np.sum(v_s*mole_fraction,axis = 1)/np.sum(np.sum(v_s*mole_fraction,axis = 1),axis = 0)    

# segment interaction matrix
# Binary interaction parameters for base segments (B-MRR1) the Extended UNISAC model
# To Be Fitted
interaction_param = np.array([[0,597,104.3, 24.9, 476.4, 1318, 0],
              [24.82,0,491.95, -15.62, -287.5, 242.8,0],
              [-78.45, -54.86,0,51.9,372, 1201,0],
              [36.7,74.04, -30.1, 0, 552.1, 826.76, 0],
              [26.76, 481.7, -39.2,-354.55, 0,472.5,0],
              [300,112.6, 497.54, 353.68, -195.4,0,0],
              [0,0,0,0,0,0,0]])
# temperature dependen segment interaction between u and v
# interaction_coeff_temp = gabel
interaction_coeff_temp = np.exp(-interaction_param/temperature)

# nuatrual logarithm of the segment activity coefficient for the mixture
# Equation 3.21
ln_gamma_mix = np.zeros(n_seg)   
for i in range(n_seg):
    first_p = 0
    second_p = 0
    first_p += np.dot(Theta_i,interaction_coeff_temp[:,i])
    for j in range(n_seg):
        second_p += Theta_i[j]*interaction_coeff_temp[i,j]/np.sum(np.dot(Theta_i,interaction_coeff_temp[:,j]),axis=0)
    ln_gamma_mix[i] = 1-np.log(first_p)-second_p

# for i in range(n_length):
# natural logarithm of the segment activity coefficient of segment k in
# the pure component i
# Equation 3.19
ln_gamma_pure = np.zeros((n_length,n_seg)).T    
for i in range(n_length):
    for j in range(n_seg):
        first_p = 0
        second_p = 0
        first_p += np.dot(Theta[:,i],interaction_coeff_temp[:,j])
        for k in range(n_seg):
            second_p += Theta[k,i] * interaction_coeff_temp[j,k] / np.sum(np.dot(Theta[:,i],interaction_coeff_temp[:,k]),axis=0)
        ln_gamma_pure[j,i] = 1-np.log(first_p)-second_p
    
# ln_gamma_mix = 1-np.log(np.sum(Theta_i*psi),axis=1)-
#                 +-np.sum(Theta_i*psi/(np.sum(Theta_i*psi),axis =0),axis=0)
residual = np.ones(2)
for i in range(n_length):
    helps = 0
    helps += np.dot(v_s[:,i],(ln_gamma_mix-ln_gamma_pure[:,i]))
    residual[i] = helps

# Equation 1.49
ln_gamma = combinatorial + residual
gamma = np.exp(ln_gamma)
print(gamma)