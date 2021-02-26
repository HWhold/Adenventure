# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:59:59 2021

@author: M302212
"""

import numpy as np
from fug import fug_minimum_parameters

class UnitOperation():
    def __init__(self, connections, flow, substances):
        self.connections = connections
        self.flow = flow
        self.substances = substances
    
class rectification(UnitOperation):
    def __init__(self, feed, product):
        fug_minimum_parameters(light_boiler, heavy_boiler, temp_bottom, temp_top, feed_purity, product_yield, product_purity)

        
