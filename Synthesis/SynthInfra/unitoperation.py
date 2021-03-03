# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:59:59 2021

@author: M302212
"""

import numpy as np
# Import Rectification Operations
from ...ProcessModels.Distillation.fug import fug_minimum_parameters, multicomponent_composition

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
    """
    A unit operation describing a rectification column. The FUG method
    is employed for all calculations.
    
    Attributes
    ----------
    AVAILABLE_COLUMNS : ndarray, dype=float, shape=(n,)
        The tray number of the available columns.
    ACTIVELY_USED : ndarray, dtype=bool, shape=(n,)
        A boolean ndarray indicating whether each available column is 
        actively used.
    product_purity : float
        The mole fraction of the product.
    product_yield : float
        The desired molar product yield.
    pressure : float (Pa)
        The pressure inside the rectification column.
    N_min : float
        The minimum number of stages to separate the product for the
        current feed.
    R_min : float
        The minimum reflux ratio to separate the product for the 
        current feed.
    N : int
        The number of stages of the currently selected column.
    """
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
        """
        Initialize the rectification column.
        
        Parameters
        ----------
        feed : flow
            The feed flow.
        product_purity : float
            The mole fraction of the product.
        product_yield : float
            The desired molar product yield.
        pressure : float (Pa)
            The pressure inside the rectification column.

        """
        # Initialize the inputs and outputs
        super.__init__(np.ndarray([False, True, True]), np.ndarray([feed, Flow(0), Flow(0)]))
        # Set the unit operations specific attributes
        self.product_purity = product_purity
        self.product_yield = product_yield
        self.pressure = pressure
        # Compute the minimum parameters of the unit operation
        self.compute_minimum_parameters()
        # Use the minimum stage number to randomly select a valid column
        possible_columns = np.where(
            (self.AVAILABLE_COLUMNS > self.N_min) & (self.ACTIVELY_USED == False) 
        )[0]
        column = np.random.choice(possible_columns)
        # Set the selected column as active, such that it is not used
        # multiple times
        self.ACTIVELY_USED[column] = True
        # Safe the number of stage of the currently selected column
        self.N = self.AVAILABLE_COLUMNS[column]
    
    def _compute_minimum_parameters(self):
        """
        Compute the minimum parameters for the given rectification 
        column.
        """
        # Rename feed for easy usage
        feed = self.inputs[0]
        # Get the feed purity
        feed_purity = feed.mole_fractions[0]
        # Arrays that store the respective vapor pressures and boiling
        # points of the substances at the given pressures and 
        # compositions.
        vapor_pressure = np.zeros(feed.SUBSTANCES.shape, dtype=float)
        boiling_point = np.zeros(feed.SUBSTANCES.shape, dtype=float)
        # Assume boiling point of product for vapor pressure estimation.
        boiling_point_product = feed.SUBSTANCES[0].boiling_point(self.pressure)
        vapor_pressure_product = feed.SUBSTANCES[0].vapor_pressure(boiling_point_product)
        # Assume Product is the most common substance, 
        # other key is second most common
        other_key_index = np.argsort(feed.molar_flows)[-2]
        
        # Compute the vapor pressure and boiling point for each 
        # substance
        for i, substance in enumerate(feed.SUBSTANCES):
            vapor_pressure[i] = substance.vapor_pressure(boiling_point_product)
            boiling_point[i] = substance.boiling_point(self.pressure)
            
        # if the product is the light boiler, calculate top and bottom
        # temperatures accordingly
        if vapor_pressure_product > vapor_pressure[other_key_index]:
            # The product is the light boiler and thus distillate
            self.product_is_distillate = True
            
            # Get the light and heavy boilder Substance objects.
            light_boiler = feed.SUBSTANCES[0]
            heavy_boiler = feed.SUBSTANCES[other_key_index]
            
            # The temperature at the top of the column is equivalent
            # to the boiling point of the product.
            self.temp_top = boiling_point_product
            
            # All substances with a lower vapor pressure are heavier 
            # than the product
            heavier_than_product = np.where(vapor_pressure < vapor_pressure_product)[0]
            # The boiling point at the bottom of the rectification 
            # columnn is assumed to be according to the fractions of
            # the heavier substances.
            self.temp_bottom = np.dot(
                feed.mole_fractions[heavier_than_product],
                boiling_point[heavier_than_product]
            )
            
        # if the product is the heavy boiler, calculate top and bottom
        # temperatures accordingly    
        else:
            # The product is the heavy boiler and thus distillate.
            self.product_is_distillate = False
            
            # Get the light and heavy boilder Substance objects.
            light_boiler = feed.SUBSTANCES[other_key_index]
            heavy_boiler = feed.SUBSTANCES[0]
            
            # The temperature at the bottom of the column is equivalent
            # to the boiling point of the product.
            self.temp_bottom = boiling_point_product
            
            # All substances with a higher vapor pressure are lighter 
            # than the product
            lighter_than_product = np.where(vapor_pressure > vapor_pressure_product)[0]
            # The boiling point at the top of the rectification 
            # columnn is assumed to be according to the fractions of
            # the lighter substances.
            self.temp_top = np.dot(
                feed.mole_fractions[lighter_than_product],
                boiling_point[lighter_than_product]
            )
        
        # Calculate the minimum values        
        self.N_min, self.R_min = fug_minimum_parameters(
            light_boiler, heavy_boiler, self.temp_bottom, self.temp_top, feed_purity, 
            self.product_yield, self.product_purity
        )
    
    def calculate_output(self):
        """
        Calculate the output of the rectification based on the input
        and column type.
        """
        # Calculate the output for each component in the feed, according
        # to the required product yield and purity.
        for i, (nonkey_substance, nonkey_feed) in enumerate(zip(self.inputs[0].SUBSTANCES, self.inputs[0].molar_flows)):
            bn, dn = multicomponent_composition(
                nonkey_substance, self.inputs[0].SUBSTANCES[0], self.temp_bottom, self.temp_top, 
                self.inputs[0].molar_fractions[0], self.product_yield, self.product_purity, nonkey_feed, 
                self.inputs[0].molar_flows[0], self.N, reference_is_distillate=self.product_is_distillate
            )
            self.outputs[0][i] = dn
            self.outputs[1][i] = bn
            
    def change_input(self, flow):
        """
        Change the input of the rectification column.
        """
        # Adjust the input parameters
        super.change_input(flow)
        # Recalculate the minimum parameters. This is neccessary as the
        # reflux has to be adjusted in order to match the product yield
        # if the feed composition changes.
        self._compute_minimum_parameters()
        
        
        
        
