# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 16:16:59 21

@author: X2722
"""


__docformat__ = "restructuredtext"
################ THIS IMPORTS HANDLER MIX COMP BENCHMARK ####################

from Adenventure.Handler.MixComp.Benchmark.consecutive import (
    consecutive
)

from Adenventure.Handler.MixComp.Benchmark.evaluation import (
    evaluation
)

from Adenventure.Handler.MixComp.Benchmark.find_parameter import (
    filter_groups,
    disassemble
)

from Adenventure.Handler.MixComp.Benchmark.find_parameter import (
    filter_groups,
    disassemble
)

from Adenventure.Handler.MixComp.Benchmark.objective import (
    objective_UNIFAC,
    objective_UNISAC,
    objective_COSMOSAC
)

from Adenventure.Handler.MixComp.Benchmark.optimization import (
     optimization
)

from Adenventure.Handler.MixComp.Benchmark.run_all import (
    run_all_UNIFAC,
    run_all_UNISAC,
    run_all_COSMOSAC_chb,
    run_all_COSMOSAC_electrostatic
)

from Adenventure.Handler.MixComp.Benchmark.run_testset import (
    run_UNIFAC_testset,
    run_UNISAC_testset,
    run_COSMOSAC_testset
)
from Adenventure.Handler.MixComp.Benchmark.run_trainingsset import (
    run_UNIFAC_trainingsset,
    run_UNISAC_trainingsset,
    run_COSMOSAC_chb_trainingsset,
    run_COSMOSAC_electrostatic_trainingsset
)
    
from Adenventure.Handler.MixComp.Benchmark.set_up import (
    load_Experiments,
    load_OrderOfFit,
    load_Data,
    load_Database
)

################ PHYS MODELS ############################

from Adenventure.PhysModels.MixComp.get_gamma_const import (
    UNIFAC_const,
    UNISAC_const
)

from Adenventure.PhysModels.MixComp.get_gamma_fit import (
    UNIFAC_fit,
    UNISAC_fit,
    run_COSMOSAC
)

from Adenventure.PhysModels.MixComp.get_gamma_inf import (
    UNIFAC_inf
)

from Adenventure.PhysModels.Thermodynamics.get_gamma_exp import (
    activity_coefficient_from_solubility
)
