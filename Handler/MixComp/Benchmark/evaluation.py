# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 10:21:12 2021
@author: Luisa Peterson

This file loads the function that statistically (including graphically) evaluates
the results of the optimization to the orginal activity-coefficient-model. This
function is applied to the trainings and the test set.

Functions overview:
    Experiments_with_Error, RSS, Stats = evaluation(gamma, fit_extrapolation,Experiments)
"""

# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
# Plotting the results
import matplotlib.pyplot as plt
# implemented statistic functions
import scipy.stats as stats

# %% EVALUTATION (Statistics and Plot)
'''
This function evaluates the quality of the model gammas in comparison with the
experimental gammas
Input:
    - gamma: experimental and model activity coefficient
    - fit_extrapolation: choose wheter the test or trainings set should be considered
    - Experiments: List of experiments that are to be considered
Return:
    - Experiments_with_Error: Experiments with model activity coeefficients and their
                              error added
    - RSS: residual sum of squares between experimental and model activity coefficients
    - Stats: Mean and standard deviation
'''
def evaluation(gamma, fit_extrapolation,Experiments):
    
    # unpack vector of gammas
    gamma_opt = gamma[0]
    gamma_opt = gamma_opt.flatten()
    gamma_exp = gamma[1]
    gamma_original = gamma[2]
    gamma_original = gamma_original.flatten()
    
    # PARITY PLOT
    plt.plot(gamma_opt,gamma_exp,'r*',label='gamma fit') # x vs y
    plt.plot(gamma_original,gamma_exp,'b*', label = 'gamma original')
    plt.plot(gamma_exp,gamma_exp+0.2*gamma_exp,'k-')  
    plt.plot(gamma_exp,gamma_exp-0.2*gamma_exp,'k-')  
    axes = plt.gca()
    gamma_min = np.amin(gamma_exp)
    gamma_max = np.amax(gamma_exp)
    axes.set_xlim([gamma_min, gamma_max])
    axes.set_ylim([gamma_min, gamma_max+0.5*gamma_max])
    plt.xlabel('model gamma')
    plt.ylabel('experimental gamma')
    plt.legend(loc='lower right')
    title = "Parityplot - "  + fit_extrapolation
    plt.title(title)
    plt.show()

    # calculate absolute Error between model and experimental activity coefficients
    Err_opt = abs(gamma_opt-gamma_exp)
    Experiments['Squared Error Fit'] = Err_opt
    Err_original = abs(gamma_original-gamma_exp)
    Experiments['Squared Error original'] = Err_original
    RSS = [sum(Err_original), sum(Err_opt)]
    
    # Output File
    Experiments['exp activity coef'] = gamma_exp
    Experiments['optimized activity coef'] = gamma_opt
    Experiments['original activity coef'] = gamma_original
    Experiments_with_Error = Experiments
    
    # Comparison to original
    eps = np.finfo(float).eps
    comparison = np.sort(np.log10(Err_opt/Err_original)+eps)
    percentage = np.arange(len(comparison))/len(comparison)*100
    plt.plot(percentage, 0*comparison, 'k-')
    plt.plot(percentage, comparison, 'g-')
    plt.xlabel('Percentage of Data Sets')
    plt.ylabel('log10(Squared Error fit / Squared Error original)')
    title = "Comparison fitted results vs original groups results - " + fit_extrapolation
    plt.title(title)
    plt.show()
    
    # Statistics (Mean, Standarddeviation)
    Stats_original = stats.describe(Err_original)
    Stats_fit = stats.describe(Err_opt)
    Stats = [Stats_original, Stats_fit]
    
    # Boxplot
    plt.boxplot((Err_opt,Err_original))
    plt.title('Deviation of the squared error in a Box Plot for '+fit_extrapolation)
    plt.xticks([1, 2], ['gamma fit - ', 'gamma original'])
    axes = plt.gca()
    axes.set_ylim([0, 2])
    plt.show()
    
    return Experiments_with_Error, RSS, Stats