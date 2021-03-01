# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:59:59 2021

@author: M302212
"""

import numpy as np
# Import Rectification Operations

# Filepath to the table with possible rectification columns
RECTIFICATION_COLUMNS = ""
SUBSTANCES = ""

class Flow():
    SUBSTANCES = load_substances()
    
    def load_substances():
        """
        Load Substances in the mixture

        Returns
        -------
        Ndarray dtype object with substances

        """
        pass
    
    def __init__(self, molar_flows, solvent=None):
        if molar_flows == 0:
            molar_flows = np.zeros(SUBSTANCES.shape, dtype=float)
        self.molar_flows = molar_flows
        self.total_flow = np.sum(molar_flows)
        self.molar_fractions = np.zeros(SUBSTANCES.shape, dtype=float)
        for i, flow in enumerate(molar_flows):
            molar_fractions[i] = flow/self.total_flow
        self.solvent = solvent
    

class UnitOperation():
    def __init__(self, connections, flow):
        self.connections = connections
        self.flow = flow
        self.inputs = flow[connections[False]]
        self.outputs = flow[connections[True]]
    
    def calculate_cost(self):
        return 0
    
    def calculate_output(self):
        # By default split inputs to outputs
        self.flow[self.connections[True]] = np.array(
            [np.sum(self.inputs, axis=0)/len(self.outputs)]*len(self.outputs)
        )
    
    def change_input(self, flow):
        self.flow[self.connections[False]] = flow
        self.calculate_output()
    
    
    
class Rectification(UnitOperation):
    def load_availabe_colums():
        """
        Load Available Rectification Columns (Stage Number) from File
    
        Returns
        -------
        Ndarray with Stage number
    
        """
        pass
    # Real Available Columns
    AVAILABLE_COLUMNS = load_availabe_colums()
    # Actively Used Columns
    ACTIVELY_USED = np.zeros(AVAILABLE_COLUMNS.shape, dtype=bool)

    def __init__(self, feed, product_purity, product_yield, pressure=1e5):
        super.__init__(np.ndarray([False, True, True]), [feed, Flow(0), Flow(0)])
        feed_purity = feed.mole_fractions[0]
        vapor_pressure = np.zeros(feed.SUBSTANCES.shape, dtype=float)
        boiling_point = np.zeros(feed.SUBSTANCES.shape, dtype=float)
        # Assume boiling point of product for vapor pressure estimation
        boiling_point_product = feed.SUBSTANCES[0].boiling_point(pressure)
        vapor_pressure_product = feed.SUBSTANCES[0].vapor_pressure(boiling_point_product)
        # Assume Product is the most common substance, other key is second most common
        other_key_index = np.argsort(feed.molar_flows)[-2]
        other_key = feed.SUBSTANCES[other_key_index]
        for i, substance in enumerate(feed.SUBSTANCES):
            vapor_pressure[i] = substance.vapor_pressure(boiling_point_product)
            boiling_point[i] = substance.boiling_point(pressure)
        if vapor_pressure_product > vapor_pressure[other_key_index]:
            product_is_distillate = True
            light_boiler = feed.SUBSTANCES[0]
            heavy_boiler = feed.SUBSTANCES[other_key_index]
            temp_top = boiling_point_product
            
            heavier_than_product = np.where(vapor_pressure < vapor_pressure_product)[0]
            temp_bottom = np.dot(
                feed.mole_fractions[heavier_than_product],
                boiling_point[heavier_than_product]
            )
            
        else:
            self.product_is_distillate = False
            light_boiler = feed.SUBSTANCES[other_key_index]
            heavy_boiler = feed.SUBSTANCES[0]
            temp_bottom = boiling_point_product
            
            lighter_than_product = np.where(vapor_pressure < vapor_pressure_product)[0]
            temp_bottom = np.dot(
                feed.mole_fractions[lighter_than_product],
                boiling_point[lighter_than_product]
            )
        
                
        self.N_min, self.R_min = fug.fug_minimum_parameters(
            light_boiler, heavy_boiler, temp_bottom, temp_top, feed_purity, 
            product_yield, product_purity
        )
        possible_columns = np.where(
            (AVAILABLE_COLUMNS > self.N_min) & (ACTIVELY_USED == False) 
        )[0]
        column = np.random.choice(possible_columns)
        self.ACTIVELY_USED[column] = True
        self.N = self.AVAILABLE_COLUMNS[column]
    
    def _compute_minimum_parameters():
        feed = self.inputs[0]
    
    
    def calculate_output(self):
        if se
        fug.multicomponent_composition()
        
        
        
        
