# -*- coding: utf-8 -*-

class substance():
    """
    This class represents the thermodynamic properties of atoms.
    
    Attributes
    ----------
    vapor_pressure_func : callable (Pa)
        A function that returns the vapor pressure for a given
        Temperature (K)
    """
    def __init__(self, vapor_pressure_func):
        self.vapor_pressure_func = vapor_pressure_func
    
    def vapor_pressure(self, temperature):
        """
        Return the vapor pressure for a given temperature.

        Parameters
        ----------
        temperature : float (K)
            The temperature for the vapor pressure to be returned

        Returns
        -------
        vapor_pressure : float (Pa)
            The vapor pressure at the given temperature.

        """
        return self.vapor_pressure_func(temperature)