# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:05:51 2021

@author: M302212
"""
import numpy as np

__name__ = "synthesis.equations.kremser"

def calc_kremser(K, x_in, x_out, V_L):
    """
    Compute the number of stages needed for a given extraction, where
    the an in and output solvent composition of a product is given.
    
    The equation is taken from:
        https://www.youtube.com/watch?v=dByYrj7-tYQ.

    Parameters
    ----------
    K : float
        The distribution coefficient of the product between the two
        solvents.
    x_in : float
        The solvent composition of the product at the beginning of the
        extraction.
    x_out : float
        The target solvent compostion of the product at the end of the
        extraction.
    V_L : float
        The feed to solvent ratio.

    Returns
    -------
    N : float
        The number of stages needed for the separation.
    """
    return np.log((1-K*V_L)*(x_in/x_out)+K*V_L)/np.log(1/V_L/K)


def calc_comp_kremser(K, x_in, n, V_L):
    """
    Compute the solvent output composition for a given substance and
    extraction operation.
    
    The equation is taken from:
        https://www.youtube.com/watch?v=dByYrj7-tYQ.

    Parameters
    ----------
    K : float
        The distribution coefficient of the substance between the two
        solvents.
    x_in : float
        The solvent composition of the substance at the beginning of the
        extraction.
    N : float
        The number of stages.
    V_L : float
        The feed to solvent ratio.

    Returns
    -------
    x_out : float
        The solvent compostion of the substance at the end of the
        extraction.
    """
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