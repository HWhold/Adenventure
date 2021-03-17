# -*- coding: utf-8 -*-
import numpy as np

__name__ = "synthesis.equations.fug"
def fug_minimum_parameters(
        light_boiler, heavy_boiler, temp_bottom, temp_top, 
        feed_purity, product_yield, product_purity, product_is_distillate=True
    ):
    """
    Compute the minimum number of trays and minimum reflux ratio for 
    given feed and product streams of binary mixtures.
    
    The FUG method is employed to determine the values. Fenske´s
    equation is used to determine the minimum number of trays and
    Underwood´s equation is used to determine the minimum reflux ratio.

    Parameters
    ----------
    light_boiler : substance
        The properties of the light boiler.
    heavy_boiler : substance
        The properties of the heavy boiler.
    temp_bottom : float (K)
        The temperature at the bottom of the rectification column.
    temp_top : float (K)
        The temperature at the top of the rectification column.
    feed_purity : float
        The mole fraction of the product in the feed.
    product_yield : float
        The desired molar product yield.
    product_purity : float
        The mole fraction of the product.
    product_is_distillate : bool (default: True)
        If set to `True`, the product will assumed to be the light 
        substance, otherwise the product will assumed to be the heavy 
        substance.

    Returns
    -------
    N_min : float
        The minimum number of trays needed for the separation.
    R_min : float
        The minimum reflux ratio needed for the separation.
    """
    
    # Rename input variables according to equation symbols
    xF = feed_purity
    Y = product_yield
    
    # Calculate the relative volatility at the bottom of the 
    # rectification column
    alpha_bottom = light_boiler.vapor_pressure(temp_bottom)/heavy_boiler.vapor_pressure(temp_bottom)
    # Calculate the relative volatility at the top of the rectification 
    # column
    alpha_top = light_boiler.vapor_pressure(temp_top)/heavy_boiler.vapor_pressure(temp_top)
    # For calculation use the geometric mean of top and bottom alpha
    # values
    alpha = np.sqrt(alpha_bottom * alpha_top)
    
    # if the product is the distillate, compute the bottom composition
    if product_is_distillate:        
        xD = product_purity
        
        # Calculate the bottom mole fraction
        xB = ( xF * (1 - Y) ) / ( 1 - ( ( xF * Y) / xD ) )
    
    # Otherwise compute the top compostion
    else:
        xB = product_purity

        # Calculate the bottom mole fraction
        xD = ( xF * (1 - Y) ) / ( 1 - ( ( xF * Y) / xB ) )
        
        # If the product is the bottom, the relative volatility has to
        # be inverted
        alpha = 1/alpha

    # Calculate the minimum number of trays
    N_min = np.log( ( xD / ( 1 - xD ) ) * ( ( 1 - xB ) / xB ) ) / np.log(alpha)

    # Calculate the minimum reflux ratio
    R_min = ( 1 / ( alpha - 1) ) * ( ( xD / xF ) - ( (alpha * (1 - xD) ) / (1 - xF) ) )
    
    return N_min, R_min
    
    
def feed_tray_position(
        feed_purity, product_yield, product_purity, tray_number, 
        product_is_distillate=True
    ):
    """
    Calculate the optimal feed tray position based on Kirkbride´s 
    equation.

    Parameters
    ----------
    feed_purity : float
        The mole fraction of the product in the feed.
    product_yield : float
        The desired molar product yield.
    product_purity : float
        The mole fraction of the product.
    tray_number : int
        The number of trays in the rectification column.
    product_is_distillate : bool (default: True)
        If set to `True`, the product will assumed to be the light 
        substance, otherwise the product will assumed to be the heavy 
        substance.

    Returns
    -------
    feed_tray_position : float
        The optimal position of the feed tray according to Kirkbride´s
        equation.
    """
    # Rename input variables according to equation symbols
    xF = feed_purity
    Y = product_yield
    N = tray_number
    
    # if the product is the distillate, compute the bottom composition
    if product_is_distillate:        
        xD = product_purity
        
        # Calculate the bottom mole fraction
        xB = ( xF * (1 - Y) ) / ( 1 - ( ( xF * Y) / xD ) )
        
        # Calculate the B/D ratio
        BtoD = (xD / (xF * Y)) - 1
    
    # Otherwise compute the top compostion
    else:
        xB = product_purity
        
        # Calculate the bottom mole fraction
        xD = ( xF * (1 - Y) ) / ( 1 - ( ( xF * Y) / xB ) )
        
        # Calculate the B/D ratio
        BtoD = (xB / (xF * Y)) - 1

    # Calculate the Kirkbride Ratio (NR/NS)
    ratio = ( ( ( 1 - xF) / xF ) * np.power(xB / ( 1- xD ), 2  ) * ( BtoD ) ) ** 0.206

    # Calculate optimal feed postion (NR)
    return ( N * ratio ) / ( ratio + 1 )

def multicomponent_composition(
        nonkey_substance, reference_substance, temp_bottom, temp_top, 
        reference_feed_purity, reference_yield, reference_purity, nonkey_feed, 
        reference_feed, N_min, reference_is_distillate=True
    ):
    """
    Compute the the composition at the top and bottom of a nonkey 
    substance for a rectification column designed for a reference light
    or heavy key substance.

    Parameters
    ----------
    nonkey_substance : substance
        The properties of the nonkey substance.
    reference_substance : substance
        The properties of the reference substance.
    temp_bottom : float (K)
        The temperature at the bottom of the rectification column.
    temp_top : float (K)
        The temperature at the top of the rectification column.
    reference_feed_purity : float
        The mole fraction of the reference substance in the feed.
    reference_yield : float
        The yield of the reference substance.
    reference_purity : float
        The purity (mole fraction) of the reference substance after 
        rectification.
    nonkey_feed : float (mol/s)
        The molar feed of the nonkey substance
    reference_feed : float (mol/s)
        The molar feed of the reference substance
    N_min : float
        The minimum number of trays needed for the separation.
    product_is_distillate : bool (default: True)
        If set to `True`, the product will assumed to be the light 
        substance, otherwise the product will assumed to be the heavy 
        substance.

    Returns
    -------
    non_key_bottom_molar_flow_rate : float (mol/s)
        The flow rate of the nonkey substance at the bottom of the
        column.
    non_key_distillate_molar_flow_rate : float (mol/s)
        The flow rate of the nonkey substance at the top of the column.
    """
    
    # Rename input variables according to equation symbols
    xF = reference_feed_purity
    Y = reference_yield
    fn = nonkey_feed
    fr = reference_feed
    
    # if the reference substance is the distillate, compute the bottom 
    # composition
    if reference_is_distillate:        
        xD = reference_purity
        
        # Calculate the total distillate molar flow
        D = (fr / xD) * Y
        
        # Calculate the distillate molar flow of the reference substance
        dr = xD * D
        
        # Calculate the bottom molar flow of the reference substance
        br = fr - dr
    
    # Otherwise compute the top compostion
    else:
        xB = reference_purity
                
        # Calculate the total distillate molar flow
        B = (fr / xB) * Y
        
        # Calculate the distillate molar flow of the reference substance
        br = xB * B
        
        # Calculate the bottom molar flow of the reference substance
        dr = fr - br
      
    # Calculate the relative volatility at the bottom of the 
    # rectification column
    alpha_bottom = nonkey_substance.vapor_pressure(temp_bottom)/reference_substance.vapor_pressure(temp_bottom)
    # Calculate the relative volatility at the top of the 
    # rectification column
    alpha_top = nonkey_substance.vapor_pressure(temp_top)/reference_substance.vapor_pressure(temp_top)
    # For calculation use the geometric mean of top and bottom alpha
    # values
    alpha = np.sqrt(alpha_bottom * alpha_top)
    
    # Calculate the  molar flow of the nonkey substance
    bn = fn / ( 1 + ( ( dr / br ) * ( alpha ** N_min ) ) )
    
    # Calculate the distillate molar flow of the nonkey substance
    dn = ( (dr / br ) * ( alpha ** N_min ) * fn ) / ( 1 + ( ( dr / br ) * ( alpha ** N_min ) ) )
    
    return bn, dn