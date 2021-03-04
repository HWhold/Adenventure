# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:05:51 2021

@author: M302212
"""
import numpy as np

def calc_kremser(K, x_in, x_out, V_L):
    # K distribution coefficient [-]
    # x_in, xout compositions start and end in [mol/mol]
    #  V_L feed to solvent ratio [-]
    # Number of Stages
    return np.log((1-K*V_L)*(x_in/x_out)+K*V_L)/np.log(1/V_L/K)


def calc_comp_kremser(K, x_in, n, V_L):
    # calculate residual composition via kremser equation when stage numbers are given
    # ATTENTION if necessary purity is high, high solvent rations are needed to get necessary
    # purification - if not enough solvent-> NAN
    # taken from 
    # K distribution coefficient [-]
    #  x_in =n_in feed mole rate [mol/mol]
    # n stages [-]
    # V_L feed to solvent rate [-]
    
    part1 = np.log(1/V_L/K)
    part2 = (np.exp(part1*n)-V_L*K)/(1-V_L*K)
    
    return x_in/part2

def get_stage_number(solubility, K, feed_purity, product_purity, V_L):
    # Calculate the solvent fraction of the product at feed
    x_in = feed_purity / ( 1 + (product_purity / solubility) )
    # Calculate the solvent fraction of the product stream
    x_out = product_purity / ( 1 + (feed_purity / solubility) )
    
    return calc_kremser(K, x_in, x_out, V_L)

def multicomponent_composition(substance, solvent1, product, temperature, product_purity, feed_purity, V_L, stages):
    # feed purity of *substance*
    # product purity of *product*
    # Calculate the solubility
    solubility = get_solubility(product, solvent, temperature)
    # Get the Equilibrium Constant
    K = get_K_value(substance, solvent, temperature)
    # Calculate the solvent fraction of the substance at feed
    x_in = feed_purity / ( 1 + (product_purity / solubility) )
    return calc_comp_kremser(K, x_in, stages, V_L)