# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 11:33:25 2021
@author: M294241

Testing the implementation of the UNIFAC model which is implemented in two parts:
UNIFAC_const+UNIFAC_fit. In the following test, UNIFAC calculates activity
coefficients for the system n-Hexane in 2-Butanone. For comparison the activity
coefficients were obtained by the Dortmund Database for the same experimental
conditions. The activity coefficients calculated from the implemented UNIFAC and
the activity coefficients from the DDBSP are compared.
"""
# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
from itertools import repeat
import pytest
# import functions
from Adenventure.PhysModels.MixComp.get_gamma_const import UNIFAC_const
from Adenventure.PhysModels.MixComp.get_gamma_fit import UNIFAC_fit
from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import load_exp_gamma
from Adenventure.Test.a_Test_UNIFAC.a_set_up_test import load_Experiments, load_Data, load_Database
from Adenventure.Handler.MixComp.Benchmark.find_parameter import filter_groups

def test_UNIFAC():
    # set up
    Experiments = load_Experiments('test','fit')
    Data = load_Data('test', 'UNIFAC')
    Database = load_Database('test', 'UNIFAC')
    Data, Database, _ = filter_groups(Experiments, Data, Database)
    interaction_param = Database[2]
    
    # get activity coefficients from Dortmund Database
    gamma_test = load_exp_gamma(Experiments)
    
    # get activity coefficients from UNIFAC
    Experiments = Experiments.values
    UNIFAC_const_out = list(map(UNIFAC_const, repeat(Data), Experiments, repeat(Database), repeat(True)))
    gamma_UNIFAC = list(map(UNIFAC_fit, repeat('None'), repeat(Data), Experiments, 
                            repeat(Database), repeat(interaction_param), 
                            repeat('None'), UNIFAC_const_out, repeat(True)))
    gamma_UNIFAC = np.array(gamma_UNIFAC)
    gamma_UNIFAC = gamma_UNIFAC.flatten()
    
    assert list(gamma_UNIFAC) == pytest.approx(list(gamma_test), abs=0.1)
    