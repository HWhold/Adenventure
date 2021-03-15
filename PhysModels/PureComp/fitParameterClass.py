#this class is used throughout the VAP Benchmark
from dataclasses import dataclass

@dataclass
class parameters:
    #Boiling point experimental or fit from experimental data externally
    Tb: float = 0
    # boiling point from DDBST sim
    Tb2: float = 0
    # boiling point fromo nannolal from DDBST calculation
    Tb3: float = 0
    # boiling point extrapolated from moller with smoothed vapor pressure
    Tb4: float = 0
    #boiling point extrapolated from moller with exp values
    Tb5: float = 0
    #boiling point extrapolated from nannoolal with experimental values
    Tb6: float = 0
    # ds parameter for moller straight from DDBST
    Moller_Bs: float = 0
    # bs parameter for moller straight from DDBST
    Moller_Ds: float = 0
    # ds parameter for Nannoolal straight from DDBST
    para_Nannoolal: float = 0
    # ds parameter for nanoolal straight from optimization
    Nann_opt: float = 0
    # # ds parameter for moller after using analytic solution with Tb2
    Ds2: float = 0
    # bs parameter for moller after using analytic solution  with Tb2
    Bs2: float = 0
    # ds parameter for moller after using analytic solution   experimental BPT
    Ds3: float = 0
     # bs parameter for moller after using analytic solution   experimental BPT
    Bs3: float = 0
    # ds parameter for Nannoolal from experimental values
    Ds4: float = 0
    # ds parameter for Nannoolal from experimental values with BPT fit
    Ds5: float = 0
    # ds parameter for moller after using analytic solution with Tb3
    Ds7: float = 0
    # ds parameter for moller after using analytic solution with Tb3
    Bs8: float = 0
    # parameter ds from fit with optimizer for moller
    para1: float = 0
    # parameters bs from fit with optimizer for moller
    para2: float = 0
    # antoine parameters
    antoine1: float = 0
    antoine2: float = 0
    antoine3: float = 0