# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:59:59 2021

@author: M302212
"""

import numpy as np
# Import Rectification Operations

# Filepath to the Excel-Sheet with possible rectification columns
RECTIFICATION_COLUMNS = ""
# Filepath to the Excel-Sheet containing the substances in the mixture
SUBSTANCES = ""

class Flow():
    """
    This class represents the molar flow of substances.
    
    Attributes
    ----------
    SUBSTANCES : ndarray, dtype=Substance, shape=(n,)
        The substances in the mixture.
    molar_flow : ndarray, dtype=float, shape=(n,) (mol/s)
        The molar flow of the individual substances in the mixture.
    molar_fractions: ndarray, dtype=float, shape=(n,)
        The molar fractions of the substances in respect to the total
        molar flow.
    total_flow : float (mol/s)
        The total molar substance flow.
    """
       
    def _load_substances():
        """
        Load substances in the mixture. The substance data is located
        in an Excel-Sheet at the location provided in the global 
        variable `SUBSTANCES`.

        Returns
        -------
        substances: ndarray, dtype=object, shape=(n,)
            The substances in the mixture.
        """
        
        pass
    
    # On load, get the substances in the mixture.
    SUBSTANCES = _load_substances()
    
    def __init__(self, molar_flows, solvent=None):
        """
        A flow is made up of the individual molar flows of each 
        substance and the solvent.

        Parameters
        ----------
        molar_flows : ndarray, dtype=flaot, shape=(n,) or 0
            The molar flow of the individual substances.
        solvent : Substance, optional (default: None)
            The solvent the substances are dissolved in.
        """
        
        # If the molar flow is zero, create an `ndarray`, representing
        # each individual flow as zero.
        if molar_flows == 0:
            molar_flows = np.zeros(SUBSTANCES.shape, dtype=float)
        
        # Set the molar flow and solvent attributes
        self.solvent = solvent
        self.molar_flows = molar_flows
        
        # Calculate the total flow
        self.total_flow = np.sum(molar_flows)
        # Calculate the molar fractions
        self.molar_fractions = np.zeros(SUBSTANCES.shape, dtype=float)
        for i, flow in enumerate(molar_flows):
            self.molar_fractions[i] = flow/self.total_flow

class UnitOperation():
    """
    This class is the base class for all Unit Operations.
    
    Attributes
    ----------
    connections : ndarray, dtype=bool, shape=(n,)
        The inputs (False) and outputs (True) of a given unit operation.
    flow : ndarray, dtype=Flow, shape=(n,)
        The flow for each connection.
    inputs : ndarray, dtype=Flow
        The inputs of the unit operation.
    outputs : ndarray, dtype=Flow
        The outputs of the unit operation.
    """
    
    def __init__(self, connections, flow):
        """
        Initialize the unit operation, given its connections and flows.

        Parameters
        ----------
        connections : ndarray, dtype=bool, shape=(n,)
            The inputs (False) and outputs (True) of a given unit 
            operation.
        flow : ndarray, dtype=Flow, shape=(n,)
            The flow for each connection.
        """
        self.connections = connections
        self.flow = flow
        # False connections correspond to inputs
        self.inputs = flow[connections[False]]
        # True connections correspond to outputs
        self.outputs = flow[connections[True]]
    
    def calculate_cost(self):
        """
        Calculate the cost of the unit operation. The base unit 
        operation is free.

        Returns
        -------
        cost : flaot
            The cost of the unit operation.
        """
        return 0
    
    def calculate_output(self):
        """
        Calculate the output of the unit operation. By default the 
        input molar stream is distributed evenly across the outputs.
        """
        self.flow[self.connections[True]] = np.array(
            [np.sum(self.inputs, axis=0)/len(self.outputs)]*len(self.outputs)
        )
        
        
    def change_input(self, flow):
        """
        Change the input of a given unit operation. The output is not
        automatically updated and has to be recalculated manually.
        
        Parameters
        ----------
        flow : ndarray, dtype=Flow, shape=(n,)
            The flow for each connection.
        """
        self.flow[self.connections[False]] = flow
    
class Feed(UnitOperation):
    """
    A unit operation describing the feed of a system of unit operations.
    """
    def __init__(self, flow):
        """
        Initialize the feed.
        
        Parameters
        ----------
        flow : ndarray, dtype=Flow, shape=(1,)
            The feed flow.
        """
        super.__init__(np.ndarray([True]), np.ndarray(flow))

class Product(UnitOperation):
    """
    A unit operation describing the product endpoint of a system of 
    unit operations.
    """
    def __init__(self, feed):
        """
        Initialize the product end point.
        
        Parameters
        ----------
        flow : ndarray, dtype=Flow, shape=(1,)
            The product flow.
        """
        super.__init__(np.ndarray([False]), np.ndarray(feed))

class Waste(UnitOperation):
    """
    A unit operation describing the waste endpoint of a system of 
    unit operations.
    """
    def __init__(self, feed):
        """
        Initialize the waste end point.
        
        Parameters
        ----------
        flow : ndarray, dtype=Flow, shape=(1,)
            The waste flow.
        """
        super.__init__(np.ndarray([False]), np.ndarray(feed))          
        
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
        super.__init__(np.ndarray([False, True, True]), np.ndarray([feed, Flow(0), Flow(0)]))
        self.product_purity = product_purity
        self.product_yield = product_yield
        self.pressure = pressure
        possible_columns = np.where(
            (self.AVAILABLE_COLUMNS > self.N_min) & (self.ACTIVELY_USED == False) 
        )[0]
        column = np.random.choice(possible_columns)
        self.ACTIVELY_USED[column] = True
        self.N = self.AVAILABLE_COLUMNS[column]
    
    def _compute_minimum_parameters(self):
        feed = self.inputs[0]
        feed_purity = feed.mole_fractions[0]
        vapor_pressure = np.zeros(feed.SUBSTANCES.shape, dtype=float)
        boiling_point = np.zeros(feed.SUBSTANCES.shape, dtype=float)
        # Assume boiling point of product for vapor pressure estimation
        boiling_point_product = feed.SUBSTANCES[0].boiling_point(self.pressure)
        vapor_pressure_product = feed.SUBSTANCES[0].vapor_pressure(boiling_point_product)
        # Assume Product is the most common substance, other key is second most common
        other_key_index = np.argsort(feed.molar_flows)[-2]
        for i, substance in enumerate(feed.SUBSTANCES):
            vapor_pressure[i] = substance.vapor_pressure(boiling_point_product)
            boiling_point[i] = substance.boiling_point(self.pressure)
        if vapor_pressure_product > vapor_pressure[other_key_index]:
            self.product_is_distillate = True
            light_boiler = feed.SUBSTANCES[0]
            heavy_boiler = feed.SUBSTANCES[other_key_index]
            self.temp_top = boiling_point_product
            
            heavier_than_product = np.where(vapor_pressure < vapor_pressure_product)[0]
            self.temp_bottom = np.dot(
                feed.mole_fractions[heavier_than_product],
                boiling_point[heavier_than_product]
            )
            
        else:
            self.product_is_distillate = False
            light_boiler = feed.SUBSTANCES[other_key_index]
            heavy_boiler = feed.SUBSTANCES[0]
            self.temp_bottom = boiling_point_product
            
            lighter_than_product = np.where(vapor_pressure < vapor_pressure_product)[0]
            self.temp_top = np.dot(
                feed.mole_fractions[lighter_than_product],
                boiling_point[lighter_than_product]
            )
        
                
        self.N_min, self.R_min = fug.fug_minimum_parameters(
            light_boiler, heavy_boiler, self.temp_bottom, self.temp_top, feed_purity, 
            self.product_yield, self.product_purity
        )
    
    def calculate_output(self):
        for i, (nonkey_substance, nonkey_feed) in enumerate(zip(self.inputs[0].SUBSTANCES, self.inputs[0].molar_flows)):
            bn, dn = fug.multicomponent_composition(
                nonkey_substance, self.inputs[0].SUBSTANCES[0], self.temp_bottom, self.temp_top, 
                self.inputs[0].molar_fractions[0], self.product_yield, self.product_purity, nonkey_feed, 
                self.inputs[0].molar_flows[0], self.N, reference_is_distillate=self.product_is_distillate
            )
            self.outputs[0][i] = dn
            self.outputs[1][i] = bn
            
    def change_input(self, flow):
        super.change_input(flow)
        self._compute_minimum_parameters()
        
        
        
        
