# -*- coding: utf-8 -*-
import pytest
import pandas as pd
from .. import fug
from ..substance import substance

def test_fug_minimum_parameters():
    """
    Test the calculation of N_min and R_min values using reference
    calculations from CHEMCAD and [1].
    
    References
    ----------
    
    ..  [1] Halvorsen, Ivar & Skogestad, Sigurd. (2000). 
        Distillation Theory. Encyclopedia of Separation Science. p.49 
        10.1016/B0-12-226770-2/00631-1. 
    """
    df = pd.read_excel("data\\fug\\binary.xlsx").dropna()
    for i, row in df.iterrows():
        light_substance = substance(
            lambda t: row['ps_light_bottom'] if t == row['temp_bottom'] else row['ps_light_top']
        )
        heavy_substance = substance(
            lambda t: row['ps_heavy_bottom'] if t == row['temp_bottom'] else row['ps_heavy_top']
        )
        N_min, R_min = fug.fug_minimum_parameters(
            light_substance, heavy_substance,
            row['temp_bottom'], row['temp_top'],
            row['feed_purity'], 
            row['product_yield'], row['product_purity'],
            product_is_distillate = bool(row['product_is_distillate'])
        )
        assert N_min == pytest.approx(row['N_min'], rel=0.03)
        assert R_min == pytest.approx(row['R_min'], rel=0.085)
        
def test_feed_tray_position():
    """
    Test the calculation of the feed tray location using reference
    calculations from CHEMCAD.
    """
    df = pd.read_excel("data\\fug\\binary.xlsx").dropna()
    for i, row in df.iterrows():
        feed_tray = fug.feed_tray_position(
            row['feed_purity'],
            row['product_yield'],
            row['product_purity'],
            row['N'],
            product_is_distillate = bool(row['product_is_distillate'])
        )
        assert feed_tray == pytest.approx(row['tray_position'], rel=0.035)
        
def test_exact_feed_tray_position():
    """
    Test the calculation of the feed tray location with exact input
    values generated from reference calculations using CHEMCAD.
    """
    feed_purity = 0.8/1
    product_yield = ((0.792/0.794)*0.794)/(0.8*1)
    product_purity = 0.792/0.794
    N = 22.9998
    feed_tray = fug.feed_tray_position(feed_purity, product_yield, product_purity, N)
    assert feed_tray == pytest.approx(14.6569, abs=1e-4)
        

def test_multicomponent_composition():
    """
    Test the calculation of the t compostion using reference
    calculations from CHEMCAD.
    """
    df = pd.read_excel("data\\fug\\multicomponent.xlsx").dropna()            
    for i, row in df.iterrows():
        nonkey_substance = substance(
            lambda t: row['ps_nonkey_bottom'] if t == row['temp_bottom'] else row['ps_nonkey_top']
        )
        reference_substance = substance(
            lambda t: row['ps_reference_bottom'] if t == row['temp_bottom'] else row['ps_reference_top']
        )
        bn, dn = fug.multicomponent_composition(
            nonkey_substance, reference_substance,
            row['temp_bottom'], row['temp_top'],
            row['feed_purity'],
            row['product_yield'],
            row['product_purity'],
            row['fi'], row['fn'], row['N_min'],
            reference_is_distillate = row['reference_is_distillate']
        )
        assert row['bi'] == pytest.approx(bn, rel=0.08)
        assert dn == pytest.approx(row['di'], rel=0.2)