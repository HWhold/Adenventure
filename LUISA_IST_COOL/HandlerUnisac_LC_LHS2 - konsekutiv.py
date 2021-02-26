"""
Created: 2021/02/23
@author: Luisa Peterson

Handler for UNISAC

The goal of the handler is to get segment area parameters for UNISAC.

UNISAC is a model to computationally obtain interaction parameter. To use UNISAC,
segment area parameters are needed. However, not always are all segment area
parameters given. The missing parameters can be found by minimizing the 
weighted squared error between the model and the experimental activity 
coefficients for a trainings set. Experimental activity coefficients are
obtaines by measured solubilities. Model activity coefficients  are obtained 
via UNISAC.

Based on the optimized parameters, a test set with models similar
to the trainings set should be fitted. For these molecules the activity 
coefficient are determined with the previously optimized parameters via UNISAC. 

The structure of the handler is based  on the following theory: 
Similarity hypothesis: Similar molecules have similar solubilities
"""
# %% IMPORT PACKAGES
'''
Program libraries with modules that are used within the script are shown below.
Either the entire library or only individual modules are loaded.
'''
# builds and contains data structures and operators for accessing numeric tables
import pandas as pd
# handling of vectors, matrices or generally large multidimensional arrays
import numpy as np
# provides functions for minimizing (or maximizing) objective functions.
from scipy.optimize import dual_annealing, minimize
# use different kernels to get faster
from optimparallel import minimize_parallel
# Plotting the results
import matplotlib.pyplot as plt
# implemented statistic functions
import scipy.stats as stats
# Factorial Sampling
from pyDOE import lhs
# For mapping function to repeat values
from itertools import repeat
# multiprocessing for parallization of starmap
import multiprocessing as mp
from multiprocessing import Pool
# supress warning
import warnings

# %% LOAD DATA
'''
In this function, experiments are imported from Excel tables
Input:
        - name_group: group of solutes so that the right file can be loaded
        - fit_extrapolation: choose wheter the test or trainings set should be considered
Output:
        - Experiments: Experiments for the group and either fit or extrapolation
'''   
def load_Experiments(name_group, fit_extrapolation):
    # load experimental data for a given solute and different solvents
    # Columns: Temperature in °C, solvent name, solute name, solubility
    path_empty = 'Experiments_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    Experiments = pd.read_excel(path, sheet_name = fit_extrapolation, delimiter = ' ', header = 0)
    
    # determine the order of fit for consecutive fitting
    path_empty = 'OrderOfFit_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    Order_Fit = pd.read_excel(path, delimiter = ' ', header = 0)
    return Experiments, Order_Fit

'''
In this function, solubility data of the molecules and their assignment to
structural UNISAC groups are loaded.
Input:
        - name_group: group of solutes so that the right file can be loaded
Output:
        - Data: solubility data and structural information for solutes and solvents
''' 
def load_Data(name_group):
    # load information about the solutes
    # Columns: Material, FileName, Formula Smiles, CAS,	MPT °C,	HFUS J/mol
    Solutes = pd.read_excel('Molecules.xlsx', sheet_name= name_group, index_col ="Material")
    
    # load information about the solvents
    # Columns: Material, Formula Smiles, CAS, Source MPT, MPT °C, Source HFUS, HFUS J/mol
    Solvents = pd.read_excel('Molecules.xlsx', sheet_name= 'Solvents', index_col ="Material")
    
    # load Structural Information
    # subgroup of each component and their counts for all species
    # structural information with own group
    path_empty = 'StructuralInformation_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    structural_information = pd.read_excel(path, delimiter = ' ', header = 0, index_col ="Name")

    # structural information where the own group is split to existing UNIFAC groups
    # used for comparison
    path_to_ctc2 = '-UNISAC.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    structural_information_UNISAC = pd.read_excel(path, delimiter = ' ', header = 0, index_col ="Name")
    
    # combine output
    Data = [Solutes,Solvents,structural_information,structural_information_UNISAC]
    return Data

'''
In this function, the UNISAC database including the UNISAC properties (van der
Waals values and segment area information) as well as the Group interaction
parameters between the different segment areas is loaded.
Input:
        - name_group: group of solutes so that the right file can be loaded
Output:
        - Database: UNISAC properties and group interaction parameter
'''
def load_Database(name_group):
    # source: http://www.aim.env.uea.ac.uk/aim/info/UNISACgroups.html
    # Access data from database
    database = {}
    excel_file = pd.ExcelFile("database.xlsx")
    # load all sheets in excel file
    for sheet_name in excel_file.sheet_names:
        database[sheet_name] = excel_file.parse(sheet_name)
    
    # van der Waals Properties and segment Area parameters
    Properties = database.get("UNISAC_Properties")
   
    # interaction parameter between different segments
    interaction_param = database.get("GroupInteractionParam")
    interaction_param.set_index("a_nm", inplace=True)
    
    # combine output
    Database = [Properties, interaction_param]
    return Database

# %% DISASSEMBLE
'''
This function dissembles the segment area parameter into 
two parts: Firstly, segment area parameter given by UNISAC and secondly, 
segment area parameter not given by UNISAC.
Input: 
        - Properties: vdW Properties and segment area paraneters
Output: 
        - segment_area_param_where: segment area parameters disassmbled to whether they
        are given or not
'''
def disassemble(Properties):
    # get the position of interaction parameter that are not given
    segment_area_param_None_where = np.array(np.where(Properties['A'] == 'None')).T
    # get the position of interaction parameter that are given
    segment_area_param_UNISAC_where = np.array(np.where(Properties['A']  != 'None')).T
    segment_area_param_where = [segment_area_param_None_where, segment_area_param_UNISAC_where]
    return segment_area_param_where

# %% EXPERIMENTAL ACTIVITY COEFFICIENT
'''
This function calculates the activity coefficients with given solubilites,
melting points and enthalpy of fusion. The calculation is done based on
Gmehling S. 395 Eutectic Systems
Input: 
        - Data: Melting Point and Enthalpy of Fusion for the molecules
        - Experiments: Experiments to be considered
Output:
        - gamma_exp: experimental activity coefficients
'''
def activity_coefficient_from_solubility (Data, Experiments):
    # load data for one solute and all considered solvents
    # Solute Data
    solutes_fit = Experiments['Solute (1)'].unique()
    solutes = Data[0].copy()
    solutes = solutes.loc[solutes_fit]
    # Experimental Data
    Experiments = Experiments.copy()
    # get solute name for experiment
    solute_name = Experiments.loc[:, 'Solute (1)']
    # match solvent table to experiments
    solute_exp = solutes.loc[solute_name,:]
    # Temperature
    Temperature = np.array(Experiments.loc[:,'Temperature'])
    # Solubility date for solutes and solvents
    solubility_solute = Experiments.loc[:,'x_1']
    solubility_solute.index = np.arange(len(solubility_solute))
    # ideal gas constant
    R = 8.31446261815324
    # convert into np.arrays
    MPT = np.array(solute_exp['MPT'])
    HFUS = np.array(solute_exp['HFUS'])
    cp = HFUS/MPT
    # Calculate activity for the one solute
    e_length = len(Experiments)
    gamma_exp_solute = []
    for e in range(e_length):
        activity_cp = -HFUS[e]/(R*Temperature[e])*(1-Temperature[e]/(MPT[e])) + (cp[e]/(R*Temperature[e]))*(MPT[e]-Temperature[e])-(cp[e]/R)*np.log(MPT[e]/Temperature[e])
        # print((cp[e]/(R*Temperature[e]))*(MPT[e]-Temperature[e])-(cp[e]/R)*np.log(MPT[e]/Temperature[e]))
        # if necesarry: add transition enthalphy
        if solute_exp['TRST'].iloc[e] != 'None' and Temperature[e] < (solute_exp['TRST'].iloc[e]):
            TRST = np.array(solute_exp['TRST'])
            HTRS = np.array(solute_exp['HTRS'])
            transition = -HTRS[e]/(R*Temperature[e])*(1-Temperature[e]/(TRST[e]))+ (cp[e]/(R*Temperature[e]))*(TRST[e]-Temperature[e])-(cp[e]/R)*np.log(TRST[e]/Temperature[e])
            activity_cp = activity_cp + transition
        gamma_exp_solute_e = np.exp(activity_cp)/solubility_solute[e]
        gamma_exp_solute_e = pd.Series(gamma_exp_solute_e)
        gamma_exp_solute.append(gamma_exp_solute_e)
    gamma_exp_solute = pd.concat(gamma_exp_solute)
    # convert into numpy array
    gamma_exp = np.array(gamma_exp_solute)
    return gamma_exp

# %% UNISAC CONSTANT
'''
This function calculates the part of UNISAC that remains constant regardless of
the segment area parameters. The function does not need to be called by the opti-
mizer.
Input:
        - Data: constains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC database
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
Output:
        - combinatorial: combinatorial part of UNIFAC and surface fraction of
'''
def UNISAC_const(Data, Experiments, Database, UNISAC):
    eps = np.finfo(float).eps
    # Unpack Input
        # Unpack Input
    if UNISAC == False:
        structInfo = Data[2].copy()
    else:
        structInfo = Data[3].copy()  
    vdW_Properties = Database[0].copy()
    vdW_Properties = vdW_Properties.iloc[:,0:4]
    # Mole fraction of solute and solvent
    # Transpose?
    mole_fraction = np.array((Experiments[3],Experiments[4]))
    
    # Solute name for experiment e
    solute_name = Experiments[1]
    # Solvent name for experiment e
    solvent_name = Experiments[2]
    structural_information = structInfo.loc[[solute_name,solvent_name],:]

    # The van der Waals properties of the two compounds can be calculated with the
    # help of the van der Waals properties of the groups
    vdW_Properties_r = vdW_Properties.loc[:,'Rk']
    r = np.tensordot(structural_information,vdW_Properties_r,[1,0])
    vdW_Properties_q = vdW_Properties.loc[:,'Qk']
    q = np.tensordot(structural_information,vdW_Properties_q,[1,0]) 
            
    # Using these van der Waals properties for the given mole fractions, the
    # following iloc are obtained for the volume-mole-ratio (V) and the surface-
    # mole-ratio (F)
    # mol volume / mol surface * sum over all mole fractions
    V_sum = np.sum(r*np.array(mole_fraction))
    F_sum = np.sum(q*np.array(mole_fraction))
    # Volume-mole-ratio
    V = r[0]/(V_sum+eps)
    # Surface-mole-ratio
    F = q[0]/(F_sum+eps)
        
    # With these data, the combinatorial part can be calculated
    # Loop over every component
    combinatorial = 1-V+np.log(V+eps)-5*q[0]*(1-(V/(F +eps))+np.log(V/(F+eps)+eps))
    return combinatorial

# %% UNISAC FIT
'''
This function calculates the second part of UNISAC. This part is dependent from
the segment area parameters and, thus, needs to be considered by the 
objective function.
Input:
        - segment_area_param_fit: interaction parameter that are not given and, thus,
        to be fitted
        - Data: contains structural information
        - Experiments: experiments to be considered
        - Database: UNISAC database
        - combinatorial: combinatorial part of UNISAC
        - segment_area_param_where: localizes where the segment area parameter
        that are missing are
        - UNIFAC:
            - True: run only with original UNIFAC groups
            - False: run with own groups
Output:
        - gamma: Activity coefficients of all species in the mixture, [-]
'''
def UNISAC_fit(segment_area_param_fit, Data, Experiments, Database, 
               combinatorial, segment_area_param_where, UNISAC):
    
    # small number to add to divisions
    eps = np.finfo(float).eps
    # length for loops
    # number of segments
    n_seg = 7
    # number of molecules in the mixture
    n_length = 2    
    
    # get already exsisting segment area params
    Database_groups = Database[0].copy()
    segment_area_param = Database_groups.iloc[:,4:11]
    # Allocate segment area parameter
    warnings.simplefilter(action='ignore', category=FutureWarning)
    if segment_area_param_fit != 'None':
        segment_area_param_missing = segment_area_param_where[0].flatten()
        segment_area_param[int(segment_area_param_missing):] = segment_area_param_fit.reshape(1,-1)
        segment_area_param = np.array(segment_area_param, dtype=np.float)
    else:
        segment_area_param = np.array(Database_groups.iloc[:,4:11], dtype=np.float)

    
    # Unpack Input
    if UNISAC == False:
        structInfo = Data[2].copy()
    else:
        structInfo = Data[3].copy()    
    
    # get Temperature    
    Experiments = Experiments.copy()
    temperature = Experiments[0]
    # get Interaction Param
    interaction_param = Database[1]
    # get Mole Fractions
    mole_fraction = np.array((Experiments[3],Experiments[4]))
    
    # Solute and Solvent for experiment e
    solute_name = Experiments[1]
    solvent_name = Experiments[2]
    # Structural information for experiment e
    structural_information = structInfo.loc[[solute_name,solvent_name],:]
    structural_information = np.array(structural_information)
    
    # total segment area of segment k in component i
    # Equation 3.18    
    v_s = np.zeros((n_seg,n_length))
    for i in range(n_length):
        for k in range(n_seg):
            v_s[k,i] =  np.dot(structural_information.T[:,i],segment_area_param.T[k,:])
    
    # Pure Component: segment area fraction of segment m in pure component i
    # Equation 3.20
    Theta = v_s/(sum(v_s)+eps)
    # Mixture: segment area fraction of segment m in the mixture  
    # Equation 3.22
    Theta_i = np.sum(v_s*mole_fraction,axis = 1)/(np.sum(np.sum(v_s*mole_fraction,axis = 1),axis = 0)+eps)    
    
    # temperature dependen segment interaction between u and v
    # interaction_coeff_temp = gabel
    interaction_coeff_temp = np.exp(-interaction_param/temperature+eps)
    
    # nuatrual logarithm of the segment activity coefficient for the mixture
    # Equation 3.21
    ln_gamma_mix = np.zeros(n_seg) 
    first_p = np.tensordot(Theta_i,interaction_coeff_temp,[0,0])
    for i in range(n_seg):
        second_p = 0
        for j in range(n_seg):
            second_p += Theta_i[j]*interaction_coeff_temp.iloc[i,j]/(np.sum(np.dot(Theta_i,interaction_coeff_temp.iloc[:,j]),axis=0)+eps)
        ln_gamma_mix[i] = 1-np.log(first_p[i]+eps)-second_p
    
    # for i in range(n_length):
    # natural logarithm of the segment activity coefficient of segment k in
    # the pure component i
    # Equation 3.19
    ln_gamma_pure = np.zeros((n_length,n_seg)).T    
    first_p = np.tensordot(Theta,interaction_coeff_temp,[0,0])
    for i in range(n_length):
        for j in range(n_seg):
            second_p = 0
            for k in range(n_seg):
                second_p += Theta[k,i] * interaction_coeff_temp.iloc[j,k]/(np.sum(np.dot(Theta[:,i],interaction_coeff_temp.iloc[:,k]),axis=0)+eps)
            ln_gamma_pure[j,i] = 1-np.log(first_p[i,j]+eps)-second_p    
        
    residual = np.tensordot(v_s,(np.array([ln_gamma_mix], ndmin=2).T-ln_gamma_pure),[0,0])
    residual[:,0]
    
    # Equation 1.49
    ln_gamma = combinatorial + residual
    gamma = np.exp(ln_gamma)
    gamma = gamma[0]
    return gamma


#%% OBJECTIVE FUNCTION
'''
The target function calls the part of UNISAC that depends on the segment area
parameters. With the given segment area parameters the activity coefficient is cal-
culated. The result is compared with an experimentally determined activity coef-
cient. For this purpose, the error square sum is determined. If this function is
called by the optimizer, it tries to find the segment area parameters so that the sum
of squared errors is minimal.
Input:
        - segment_area_param_fit: Parameter to be fitted
        - Data_all: contains structural information
        - Experiments: experiments to be considered
        - Database: UNIFAC databse
        - combinatorial: reults from the constant part of UNIFAC
        - interaction_param_where: localizes where the segment area parameter
        that are missing are
        - gamma_exp: activity coefficient from experimental solubilities
Output:
        - RSS: Residual Sum of Square (Sum over the squared error of all experiments)
'''
def objective(segment_area_param_fit, Data, Experiments, Database, 
              combinatorial, segment_area_param_where, gamma_exp):
    
    # Get gammas for interaction params that are to be fitted
    Exp_UNISAC = Experiments.values
    gamma_mod = list(map(UNISAC_fit,repeat(segment_area_param_fit),repeat(Data),
                         Exp_UNISAC, repeat(Database), combinatorial,
                         repeat(segment_area_param_where),repeat(False)))
    gamma_mod = np.array(gamma_mod)
    gamma_mod = gamma_mod[:,0] 
    
    # penalty 1: prohibit extreme interaction parameter       
    # minimal allowed value for interaction parameter
    segment_area_param_min = 0
    # maximal allowed value for interaction parameter
    segment_area_param_max = 5
    penalty = np.zeros(len(segment_area_param_fit))
    # punish larger values
    segment_area_param_greater_where = np.where(segment_area_param_fit > segment_area_param_max)
    penalty[segment_area_param_greater_where] = (segment_area_param_fit-segment_area_param_max)[segment_area_param_greater_where]
    # punish smaller values
    segment_area_param_smaller_where = np.where(segment_area_param_fit < segment_area_param_min)
    penalty[segment_area_param_smaller_where] = (segment_area_param_min-segment_area_param_fit)[segment_area_param_smaller_where]
    # sum over all penalty values
    penalty= sum(penalty)
    
    # Error between model and experimental gamma
    # add panalty functions
    Err = abs(gamma_mod-gamma_exp)+penalty
    
    # Sum of Squares divided
    No_exp = len(gamma_mod)    
    # get Error for every experiment
    RSS = np.sum(Err)/No_exp
    return RSS

# %% OPTIMIZATION
'''
Set Up of the Optimization Problem:
Initial values and boundary conditions that are given to the optimizer.
Input:
    - Vector of inputs:
        - Data
        - Experiments
        - Database
        - UNIFAC_const_out
        - gamma_exp
        - interaction_param
        - interaction_param_where
    - interaction_param_fit_0: initial values for the solver
    - bnds: boundaries on the fitting parameters
    - solver: selection of the solver
Output:
    - solution: output of the optimizer including fitted interaction parameters,
    final Residual Sum of Squares and stats from the optimizer
'''
def optimization(inputs, segment_area_param_fit_0, bnds, solver):    
    # RUN OPTIMIZATION
    Data, Exp, Database, combinatorial, segment_area_param_where, gamma_exp = inputs
    # vector of arguments
    args = (Data, Exp, Database, combinatorial, segment_area_param_where, gamma_exp)
    # different solvers can be chosen by function input
    if solver == 'L-BFGS-B':
        ## MINIMIZE PARALLEL
        solution = minimize_parallel(objective, segment_area_param_fit_0, args = args, bounds=bnds,
                                     options={'maxfun': 150, 'maxiter': 150},
                                     parallel={'loginfo': True, 'time': True})
    if solver == 'Simplex':
        ## SIMPLEX
        interaction_param_fit_0 = np.array(segment_area_param_fit_0)
        solution =  minimize(objective, interaction_param_fit_0, method='nelder-mead', args = args,
                             options={'maxiter': 100, 'maxfev': 100, 'disp': True})
        
    if solver == 'DualAnnealing':
        ## DUAL ANNEALING
        solution = dual_annealing(func=objective, bounds=bnds, args = args, x0 = segment_area_param_fit_0)
    print('I was being called')
    return solution

# %% EVALUTATION (Statistics and Plot)
'''
This function evaluates the quality of the model gammas in comparison with the
experimental gammas
Input:
    - gamma: experimental and model activity coefficient
    - fit_extrapolation: choose wheter the test or trainings set should be considered
    - Experiments: List of experiments that are to be considered
Return:
    - Experiments_with_Error: Experiments with model activity coeefficients and their
                              error added
    - RSS: residual sum of squares between experimental and model activity coefficients
    - Stats: Mean and standard deviation
'''
def evaluation(gamma, fit_extrapolation,Experiments):
    
    # unpack vector of gammas
    gamma_opt = gamma[0]
    gamma_opt = gamma_opt.flatten()
    gamma_exp = gamma[1]
    gamma_UNISAC = gamma[2]
    gamma_UNISAC = gamma_UNISAC.flatten()
    
    # PARITY PLOT
    plt.plot(gamma_opt,gamma_exp,'r*',label='gamma fit') # x vs y
    plt.plot(gamma_UNISAC,gamma_exp,'b*', label = 'gamma UNISAC')
    plt.plot(gamma_exp,gamma_exp+0.2*gamma_exp,'k-')  
    plt.plot(gamma_exp,gamma_exp-0.2*gamma_exp,'k-')  
    axes = plt.gca()
    gamma_min = np.amin(gamma_exp)
    gamma_max = np.amax(gamma_exp)
    axes.set_xlim([gamma_min, gamma_max])
    axes.set_ylim([gamma_min, gamma_max+0.5*gamma_max])
    plt.xlabel('model gamma')
    plt.ylabel('experimental gamma')
    plt.legend(loc='lower right')
    title = "Parityplot - "  + fit_extrapolation
    plt.title(title)
    plt.show()

    # calculate absolute Error between model and experimental activity coefficients
    Err_opt = abs(gamma_opt-gamma_exp)
    Experiments['Squared Error Fit'] = Err_opt
    Err_UNISAC = abs(gamma_UNISAC-gamma_exp)
    Experiments['Squared Error UNISAC'] = Err_UNISAC
    RSS = [sum(Err_UNISAC), sum(Err_opt)]
    
    # Output File
    Experiments['exp activity coef'] = gamma_exp
    Experiments['optimized activity coef'] = gamma_opt
    Experiments['UNISAC activity coef'] = gamma_UNISAC
    Experiments_with_Error = Experiments
    
    # Comparison to UNISAC
    comparison = np.sort(np.log10(Err_opt/Err_UNISAC)+eps)
    percentage = np.arange(len(comparison))/len(comparison)*100
    plt.plot(percentage, 0*comparison, 'k-')
    plt.plot(percentage, comparison, 'g-')
    plt.xlabel('Percentage of Data Sets')
    plt.ylabel('log10(Squared Error fit / Squared Error UNISAC)')
    title = "Comparison fitted results vs UNISAC groups results - " + fit_extrapolation
    plt.title(title)
    plt.show()
    
    # Statistics (Mean, Standarddeviation)
    Stats_UNISAC = stats.describe(Err_UNISAC)
    Stats_fit = stats.describe(Err_opt)
    Stats = [Stats_UNISAC, Stats_fit]
    
    # Boxplot
    plt.boxplot((Err_opt,Err_UNISAC))
    plt.title('Deviation of the squared error in a Box Plot for '+fit_extrapolation)
    plt.xticks([1, 2], ['gamma fit - ', 'gamma UNISAC'])
    axes = plt.gca()
    axes.set_ylim([0, 2])
    plt.show()
    
    return Experiments_with_Error, RSS, Stats

# %% RUN ALL
'''
This function runs all subfunctions in the right order
inputs:
    - seed: check reproducibility at different set of starting values
    - name group: name to load the right excel tables
    - statistics: boolean input whether statistic evaluation is wanted or not
outputs:
    - segment_area_param: updated matrix of segment area params with optimized values
    - Experiments_with_Error: Table of Experiments with activity coefficient
    - RSS: residual sum of squares for original UNIFAC and UNIFAC with own groups
    - Stats: Table of statistical evaluation (mean value, standard deviation)
'''
def run_all(seed, name_group, statistics):
    
    # Load Data andExperiments
    Exp, Order_Fit = load_Experiments(name_group,'fit')
    Exp_fit = Exp.values
    Data = load_Data(name_group)
    Database = load_Database(name_group)
    
    # disassemble segment area parameters that are unknown from those who are known
    Properties = Database[0]
    segment_area_param_where = disassemble(Properties)
    
    # Calculate activity coefficients from solubilites
    gamma_fit_exp = activity_coefficient_from_solubility(Data, Exp)
    
    # Run constant UNISAC
    # own groups
    combinatorial = list(map(UNISAC_const, repeat(Data), Exp_fit, repeat(Database), repeat(False)))
    
    # Set up optimization problem
    # bounds on interaction parameters
    bnds_min = 0
    bnds_max = 2.5
    bnds_new = [(bnds_min, bnds_max)]
    bnds_new = pd.Series(bnds_new)
    bnds = []
    for i in range (int(len(segment_area_param_where[0])*7)):
        bnds.append(bnds_new)
    bnds = pd.concat(bnds)
    bnds = tuple(bnds)
    
    ## number of sets of intial values
    dim = len(segment_area_param_where[0])*7
    num = 2
    
    # create semi-random initial values via LHS
    np.random.seed(seed)
    LHS_design = lhs(dim, samples=num)
    slope = bnds_max - bnds_min
    offset = bnds_min
    SLOPE = np.ones((dim,num))
    OFFSET = np.ones((dim,num))
    for i in range(num):
        SLOPE[:,i] = np.ones(dim)*slope
        OFFSET[:,i] = np.ones(dim)*offset
    segment_area_param_sampling = SLOPE.T*LHS_design+OFFSET.T
    
    ## Optimization
    # initial value of obj funciton
    obj = 1e1000
    # vector of summarized inputs
    inputs = [Data, Exp, Database, combinatorial,
              segment_area_param_where, gamma_fit_exp]
    # call BFGS optimizer
    solutions_BFGS = list(map(optimization,repeat(inputs), segment_area_param_sampling, 
                         repeat(bnds), repeat('L-BFGS-B')))

    # go trough all solution of the optimizer
    for solution in solutions_BFGS:
        # accept solution if better than the best yet existing
        if solution.fun < obj:
            obj = solution.fun
            solution_save = solution
    
    # sort optimization results for best solutions
    segment_area_param_fit_BFGS = []
    for solution in solutions_BFGS:
        segment_area_param_fit_BFGS.append([solution.x, solution.fun])
    segment_area_param_fit_BFGS = np.concatenate(segment_area_param_fit_BFGS)
    segment_area_param_fit_BFGS = np.reshape(segment_area_param_fit_BFGS, (-1, 2))
    segment_area_param_fit_BFGS = segment_area_param_fit_BFGS[segment_area_param_fit_BFGS[:,1].argsort()]
    segment_area_param_fit_BFGS = np.delete(segment_area_param_fit_BFGS,1,1)
   # select best 10% of the solutions obtained by the BFGS solver (ceiled to the next integer)
    num_simplex = np.ceil(0.1*num)
    segment_area_param_fit_BFGS = segment_area_param_fit_BFGS[0:int(num_simplex),:]

    # get parameter as inputs for starmap
    params = []
    for segment_area_param_fit in segment_area_param_fit_BFGS:
        segment_area_param_fit = np.vstack(segment_area_param_fit).astype(np.float)
        params.append(tuple((inputs, segment_area_param_fit, bnds, 'Simplex')))
        
    # Run optimization with Simplex for the results of the BFGS optimizer
    # set num_workers equal to the number of CPU cores
    num_workers = mp.cpu_count() 
    # create a pool of "num workers" subprocesses
    with Pool(num_workers) as p:    
        solutions_Simplex_starmap = p.starmap(optimization, params)
    solutions_Simplex = list(solutions_Simplex_starmap)
    
    # go through all solutions obtained by the simplex solver
    for solution in solutions_Simplex:
        # accept solution if better than the best yet existing
        if solution.fun < obj:
            obj = solution.fun
            solution_save = solution
            print("new Simplex solution accepted")
        
    #save best fitting interaction parameters
    solution = solution_save
    solution_x = solution.x.reshape(len(solution.x),1)
    
    # save best solution and put it into the database
    Database_groups = Database[0].copy()
    segment_area_param = Database_groups.iloc[:,4:11]
    segment_area_param_missing = segment_area_param_where[0].flatten()
    segment_area_param.iloc[segment_area_param_missing,:] = solution_x.T
    Database_groups.iloc[:,4:11] = segment_area_param
    Database[0] = Database_groups

    ### ALL EXPERIMENTS FROM FIT GORUP
    # Run results for all Experiments in the Fit group
    # Run with own groups
    # run constant part of UNIFAC
    gamma_fit_opt = list(map(UNISAC_fit,repeat('None'),repeat(Data),
                         Exp_fit, repeat(Database), combinatorial,
                         repeat(segment_area_param_where),repeat(False)))
    # run fit part of UNISAC with optimized parameter
    gamma_fit_opt = np.array(gamma_fit_opt)
    gamma_fit_opt = gamma_fit_opt[:,0] 

    # Run with UNISAC groups
    # run constant part of UNISAC
    combinatorial_UNISAC = list(map(UNISAC_const, repeat(Data), Exp_fit, repeat(Database), repeat(True)))
    # run fit part of UNISAC with UNISAC parameter
    gamma_fit_UNISAC = list(map(UNISAC_fit,repeat('None'),repeat(Data),
                         Exp_fit, repeat(Database), combinatorial_UNISAC,
                         repeat(segment_area_param_where),repeat(True)))
    gamma_fit_UNISAC = np.array(gamma_fit_UNISAC)
    gamma_fit_UNISAC = np.array(gamma_fit_UNISAC)
    gamma_fit_UNISAC = gamma_fit_UNISAC[:,0] 
    
    # evaluate the results - if wanted
    if statistics:
        gamma_fit = [gamma_fit_opt, gamma_fit_exp, gamma_fit_UNISAC]
        Experiments_Error_fit, RSS_fit, Stats_fit = evaluation(
        gamma_fit, 'fit', Exp)
    else:
        Experiments_Error_fit = {}
        RSS_fit = {}
        Stats_fit = {}
           
    ### EXPERIMENTS FROM EXTRAPOLATION GROUP
    # Run UNISAC fit with similar molecules and optimized interaction parameters
    # Load extrapolation Experiments
    Exp, _ = load_Experiments(name_group,'fit')
    Exp_ext = Exp.values
    
    ## experimentally obtained activity coefficients
    gamma_ext_exp = activity_coefficient_from_solubility(Data, Exp)
    gamma_ext_exp = np.array(gamma_ext_exp)
    
    ## Run with own groups
    # run constant part of UNISAC
    combinatorial = list(map(UNISAC_const, repeat(Data), Exp_ext, repeat(Database), repeat(False)))
    # fit UNISAC with optimized interaction paramters
    gamma_ext_opt = list(map(UNISAC_fit,repeat('None'),repeat(Data),
                         Exp_ext, repeat(Database), combinatorial,
                         repeat(segment_area_param_where),repeat(False)))
    gamma_ext_opt = np.array(gamma_ext_opt)
    gamma_ext_opt = gamma_ext_opt[:,0] 

    # Run with own group split to UNIFAC groups
    # run combinatorial part of UNISAC
    combinatorial_UNISAC = list(map(UNISAC_const, repeat(Data), Exp_ext, repeat(Database), repeat(True)))
    # fit UNISAC, all interaction parameters given by UNISAC
    gamma_ext_UNISAC = list(map(UNISAC_fit,repeat('None'),repeat(Data),
                         Exp_ext, repeat(Database), combinatorial_UNISAC,
                         repeat(segment_area_param_where),repeat(True)))
    gamma_ext_UNISAC = np.array(gamma_ext_UNISAC)
    gamma_ext_UNISAC = gamma_ext_UNISAC[:,0] 
    
    # evaluate the results - if wanted
    if statistics:
        gamma_ext = [gamma_ext_opt, gamma_ext_exp, gamma_ext_UNISAC]
        Experiments_Error_ext, RSS_ext, Stats_ext = evaluation(
        gamma_ext, 'fit', Exp)
    else:
        Experiments_Error_ext = {}
        RSS_ext = {}
        Stats_ext = {}
    
    # revise output data
    Experiments_with_Error = [Experiments_Error_fit,Experiments_Error_ext]
    RSS = [RSS_fit,RSS_ext]
    RSS = pd.DataFrame(RSS)
    gamma = [gamma_fit, gamma_ext]
    gamma = pd.DataFrame(gamma)
    RSS.columns = ['UNISAC params', 'optimized params']
    RSS.index = ['Fit', 'Extrapolation']
    Stats = [Stats_fit,Stats_ext]
    return segment_area_param, Experiments_with_Error, RSS, Stats

# %% MAIN
# If called as main programm
if __name__ == "__main__":
    # different groups
    name_group = 'LC2Kerner'
    # name_groups.append('OLEDs')
    # small number
    eps = np.finfo(float).eps
    # how often should the handler run?
    rep = 10
    # seed to get a different set of starting values each time the handler runs
    seed = np.zeros(rep)
    for i in range(rep):
        seed[i] = np.random.randint(1e9)
    seed = seed.astype(int)
        
    # Run Handler
    ListOfResults = list(map(run_all,seed,repeat(name_group),repeat(True)))
    
    print('Hannes, kannst du die Bilder abspeichern?')