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

#%% Set UP
# binary System
n_comp = 2
# temperature
T_system = 298.15
# mole fraction
composition = np.array([.0448, 1-.0448])
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

# first load data here from rarey email segments and values 
# index gibt nicht 0 Werte an
ind = np.array([0,7,14,15,16,17,18,20,21,35,36,37,42,43,44])
val = [0.848,0.54,0.2012,0.0162,0.0198,0.007,0.058,0.0041,0.228,0.008,0.1985,0.6914,0.008,0.1985,0.6914]

# molecule segments use
# UNISAC groups
v = np.array([[1,0,4,1,1,1,1],[2,1,0,0,0,0,1]])
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
V= r/sum(composition*r)

# Equation 3.16
q = np.sum(v2*Q,axis = 1)
# molecule vdW volume
# Equation 3.14
F = q/sum(composition*q)

# Equation 5.2
ln_gamma_comb = 1-V+np.log(V)-5*q*(1-V/F+np.log(V/F))

## RESIDUAL - NEWS START HERE
zeitanfang = time.time()
# calculate v_s (= Omega)
# total segment area of segment k in component i
# Equation 3.18
iterator = 0
v_s = np.ones(v.T.shape)
for i in v:
    v_s[:,iterator] = v[iterator,:].T@Q_seqs_rs
    iterator = iterator +1

# Pure Component
# segment area fraction of segment m in pure component i
# Equation 3.20
theta = v_s/np.sum(v_s,axis = 0)
for i in range(n_comp):
    help1 = 0
    for j in range(n_seg):
        helps = 0
        for k in range(n_seg):
            helps = helps + v_s[k,i]
        theta[j,i] = v_s[j,i]/helps
        
# Mixture       
# segment area fraction of segment m in the mixture  
# Equation 3.22
theta_m = np.sum(v_s*composition,axis = 1)/np.sum(np.sum(v_s*composition,axis = 1),axis = 0)    
for i in range(n_seg):
    helps = 0
    help1 = 0
    for j in range(n_comp):
        helps = helps+ v_s[i,j]*composition[j] 
        for k in range(n_seg):
            help1 = help1 + v_s[k,j]*composition[j]
    theta_m[i]= helps/help1

#%% Dependent on a_s

# segment interaction matrix
# Binary interaction parameters for base segments (B-MRR1) the Extended UNISAC model
# To Be Fitted
a_s= np.array([[0,597,104.3, 24.9, 476.4, 1318, 0],
              [24.82,0,491.95, -15.62, -287.5, 242.8,0],
              [-78.45, -54.86,0,51.9,372, 1201,0],
              [36.7,74.04, -30.1, 0, 552.1, 826.76, 0],
              [26.76, 481.7, -39.2,-354.55, 0,472.5,0],
              [300,112.6, 497.54, 353.68, -195.4,0,0],
              [0,0,0,0,0,0,0]])
# temperature dependen segment interaction between u and v
# psi = gabel
psi = np.exp(-a_s/T_system)

       
# for i in range(n_comp):
# natural logarithm of the segment activity coefficient of segment k in
# the pure component i
# Equation 3.19
ln_gamma_pure = np.zeros((n_comp,n_seg)).T    
for i in range(n_comp):
    for j in range(n_seg):
        first_p = 0
        second_p = 0
        for k in range(n_seg):
            first_p = first_p+ theta[k,i]*psi[k,j]
            nenner = 0
            for l in range(n_seg):
                nenner = nenner+ theta[l,i] * psi[l,k]
            second_p = second_p + theta[k,i] * psi[j,k] / nenner
        ln_gamma_pure[j,i] = 1-np.log(first_p)-second_p

# nuatrual logarithm of the segment activity coefficient for the mixture
# Equation 3.21
ln_gamma_mix = np.zeros(n_seg)   
for i in range(n_seg):
    first_p = 0
    second_p = 0
    for j in range(n_seg):
        first_p = first_p+ theta_m[j]*psi[j,i]
        nenner = 0
        for k in range(n_seg):
            nenner = nenner + theta_m[k]*psi[k,j]
        second_p = second_p + theta_m[j]*psi[i,j]/nenner
    ln_gamma_mix[i] = 1-np.log(first_p)-second_p
    
# ln_gamma_mix = 1-np.log(np.sum(theta_m*psi),axis=1)-
#                 +-np.sum(theta_m*psi/(np.sum(theta_m*psi),axis =0),axis=0)
ln_gamma_res = np.ones(2)
iterator = 0
    
for i in range(n_comp):
    helps = 0
    for j in range(n_seg):
        helps = helps + v_s[j,i]*(ln_gamma_mix[j]-ln_gamma_pure[j,i])
    ln_gamma_res[i] = helps

# Equation 1.49
ln_gamma = ln_gamma_comb + ln_gamma_res
gamma = np.exp(ln_gamma)
print(gamma)
zeitende = time.time()
print("Dauer Programmausführung:",)
print(zeitende-zeitanfang)