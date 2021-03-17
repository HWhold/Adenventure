# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 15:46:03 2021

@author: M302212
"""
import pytest
import pandas as pd
from ..kremser import calc_kremser, calc_comp_kremser

def test_calc_kremser():
    """
    Test the calculation of the stage number according to the following
    video:
        https://www.youtube.com/watch?v=dByYrj7-tYQ.
    """
    df = pd.read_excel("data\\kremser\\kremser_solvent.xlsx").dropna()
    for i, row in df.iterrows():
        K = row['K']
        x_in = row['x_in']
        x_out = row['x_out']
        V_L = row['V_L']
        n_true = row['n']
        
        n_calc = calc_kremser(K, x_in, x_out, V_L)
        assert n_calc == pytest.approx(n_true, rel=0.01)
        
def test_calc_comp_kremser():
    """
    Use the stage number to calculate the x_out values. Test data is
    taken from:
        https://www.youtube.com/watch?v=dByYrj7-tYQ
    """
    df = pd.read_excel("data\\kremser\\kremser_solvent.xlsx").dropna()
    for i, row in df.iterrows():
        K = row['K']
        x_in = row['x_in']
        x_out_true = row['x_out']
        V_L = row['V_L']
        n = row['n']
        
        x_out_calc = calc_comp_kremser(K, x_in, n, V_L)
        assert x_out_calc == pytest.approx(x_out_true, rel=0.01)