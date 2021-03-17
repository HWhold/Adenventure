# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 09:00:51 2021
Luisa Peterson
This program is used to validate the equation for the calculation of SLE as in
Gmehling Chemical Thermodynamics for Process simulation page 394 (equation 8.7).
The data for validation is from Neau et al. 1977. In their paper "differential
molar heat capacities to test ideal solubility estimations", they compared,
amongst other things, the ideal solubilty of Acetominophin. All data is taken
from the paper.
"""
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
import pytest

def test_SLE():
    R = 8.31446261815324
    HFUS = 27*1000
    MPT = 441.7
    cp_exp = 99.8
    cp_approx = 60.9
    cp_null = 0
    Temperature = 298.2
    
    cp = cp_exp
    activity_cp_exp = -HFUS/(R*Temperature)*(1-Temperature/(MPT)) + (cp/(R*Temperature))*(MPT-Temperature)-(cp/R)*np.log(MPT/Temperature)
    
    cp = cp_approx
    activity_cp_approx = -HFUS/(R*Temperature)*(1-Temperature/(MPT)) + (cp/(R*Temperature))*(MPT-Temperature)-(cp/R)*np.log(MPT/Temperature)
    
    cp = cp_null
    activity_cp_null = -HFUS/(R*Temperature)*(1-Temperature/(MPT)) + (cp/(R*Temperature))*(MPT-Temperature)-(cp/R)*np.log(MPT/Temperature)
    
    activity_test = [activity_cp_exp, activity_cp_approx, activity_cp_null]
    activity_literature = [-2.464, -2.878, -3.525]
    assert list(activity_test) == pytest.approx(list(activity_literature), abs=0.1)