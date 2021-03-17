"""
Created: 2020/11/11
@author: Luisa Peterson

Validiert, wenn CH2-CH2CO (1-9) = 476.4 und CH2CO-CH2 (9-1) = 26.76

Handler for UNIFAC

The goal of the handler is to obtain interaction parameters between molecule 
and solvent for certain molecules. The activity is given by experimental results. 
With the activity, experimental activity coefficients can be calculated. 
Model activity coefficients are obtained via UNIFAC. To use UNIFAC, interaction 
parameters are needed. Some interaction parameters are given
while others are not.

Thus, the interaction parameters can be found by minimizing the weighted squared 
error between the model and the experimental interaction parameter. Based on these
interaction parameters, molecules similar to the original molecule should be fitted.
For these molecules the activity coefficient is determined with the previously 
determined interaction parameters. The structure of the handler is based 
on the following theory: Similarity hypothesis: Similar molecules have similar solubilities
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
from smt.sampling_methods import LHS
# Standard Scaler
from sklearn.preprocessing import StandardScaler

# %% LOAD DATA
'''
The data is imported from Excel tables. The data is imported from Excel tables.
The tables are saved as Pandas.frames
'''    
def load_Data(name_group):
    # load Structural Information
    # subgroup of each component and their counts for all species
    path_empty = 'Data_Test/05_StructuralInformation_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    structural_information = pd.read_excel(path, delimiter = ' ', header = 0, index_col ="Name")
    # only consider the subgroups that are part of used solutes and solvents
    # get considered subgroups
    FilterSubgroups = structural_information.loc[:, (structural_information != 0).any(axis=0)]
    FilterSubgroups = np.array(list(FilterSubgroups))
    FilterSubgroups = np.subtract(FilterSubgroups,1.0)
    # apply subset of subgroups to structural information
    structural_information = structural_information.iloc[:,FilterSubgroups]
    # combine output
    structural_information
    return structural_information

def load_Experiments(name_group, fit_extrapolation):
    # load experimental data for a given solute and different solvents
    # Columns: Temperature in °C, solvent name, solute name, solubility
    path_empty = 'Data_Test/05_Experiments_'
    path_to_ctc2 = '.xlsx'
    path = path_empty + str(name_group) + path_to_ctc2
    Experiments = pd.read_excel(path, sheet_name = fit_extrapolation, delimiter = ' ', header = 0)
    # Conversion from °C to Kelvin
    Experiments['Temperature'] = Experiments['Temperature']
    return Experiments

def load_Database(structural_information, name_group):
    # source: http://www.aim.env.uea.ac.uk/aim/info/UNIFACgroups.html
    # Filter again Subgroups by used subgroups in structural information
    FilterSubgroups = structural_information.loc[:, (structural_information != 0).any(axis=0)]
    FilterSubgroups = np.array(list(FilterSubgroups))
    FilterSubgroups = np.subtract(FilterSubgroups,1.0)

    # Access data from database
    database = {}
    excel_file = pd.ExcelFile("Data_Test/05_database.xlsx")
    # load all sheets in excel file
    for sheet_name in excel_file.sheet_names:
        database[sheet_name] = excel_file.parse(sheet_name)
        
    # Subset of van der Waals Properties
    # vdW Properties are numberered again
    vdW_Properties = database.get("UNIFAC_vanDerWaalsProperties").copy()
    vdW_Properties = vdW_Properties.iloc[FilterSubgroups,:]
    Group_ID = vdW_Properties.Group_ID.unique()
    Group_ID_new = np.arange(len(Group_ID))
    vdW_Properties_Group_ID = vdW_Properties.loc[:,'Group_ID'].replace(to_replace = Group_ID, value = Group_ID_new)
    vdW_Properties.loc[:,'Group_ID'] = np.array(vdW_Properties_Group_ID)
    vdW_Properties.set_index(["Group_ID","Subgroup"], inplace=True)
    
    # Translates Subgroups to Maingroups
    # Subset ist used for groups and main groups and numeration is changed
    group_to_mainGroup = database.get("group_to_mainGroup").copy()
    group_to_mainGroup.set_index("Substructure", inplace=True)
    group_to_mainGroup.sort_values(by=['ID'], inplace=True) # sort by unique subgroup ID
    group_to_mainGroup = group_to_mainGroup.iloc[FilterSubgroups,:]
    group_to_mainGroup['ID'] = np.arange(len(group_to_mainGroup))
    MainGroups = group_to_mainGroup.Main_Group_ID.unique()
    MainGroups_new = np.arange(len(MainGroups))
    group_to_mainGroup['Main_Group_ID'] = group_to_mainGroup['Main_Group_ID'].replace(MainGroups, MainGroups_new)
    
    # Subset of interaction parameter
    # This values are fitted but used for the comparison of the fitted values
    interaction_param = database.get("UNIFAC_GroupInteractionParam").copy()
    interaction_param.set_index("a_nm", inplace=True)
    interaction_param = interaction_param.iloc[Group_ID-1,Group_ID-1]
    interaction_param.index = np.arange(len(interaction_param))
    interaction_param.columns = interaction_param.index
    Database = [vdW_Properties, group_to_mainGroup, interaction_param]
    return Database

# %% DISASSEMBLE
'''
This function dissembles the interaction parameter for the main groups into 
two parts: Firstly, interaction parameters given by UNIFAC and secondly, interaction
parameters not given by UNIFAC. Those can be interaction parameter that are not
available for the considered main groups or interaction parameters for additional
groups.
Input: interaction_param: Interaction params from UNIFAC consortium including Nones
Output: interaction_param_where: interaction parameters disassmbled to whether they
are given or not
'''
def disassemble(interaction_param):
    # get the position of interaction parameter that are not given
    interaction_param_None_where = np.array(np.where(interaction_param == 'None')).T
    # get the position of interaction parameter that are given
    interaction_param_UNIFAC_where = np.array(np.where(interaction_param != 'None')).T
    interaction_param_where = [interaction_param_None_where, interaction_param_UNIFAC_where]
    return interaction_param_where

# %% UNIFAC CONSTANT
'''
This function calculates the part of UNIFAC that remains constant regardless of
the interaction parameters. The function does not need to be called by the opti-
mizer.

Input:
- Data:
    - structural information (UNIFAC groups) for all solvents and solute
- Experiments: Solute and Solventnames
- Database:
    - van der Waals properties for each Subgroup
    - group to main Group translation for activation coefficients
- molecule: name of solute
Output:
- UNIFAC_const_out
    - combinatorial UNIFAC
    - Theta
    - Theta_i
'''
def UNIFAC_const(structural_information, Experiments, Database):
    # Unpack Input
    structInfo = structural_information.copy()
    Experiments = Experiments.copy()
    vdW_Properties = Database[0].copy()
    # Mole fraction of solute and solvent
    xi = np.array((Experiments.loc[:,'x_1'],Experiments.loc[:,'x_2']))
    
    # length for loops
    i_length = 2
    k_length = len(structInfo.columns)
    m_length = len(vdW_Properties)
    e_length = len(Experiments)
    
    # Loop over every experiment
    for e in range(e_length):
        # Set up parameters
        # mole fraction for experiment e
        mole_fraction = xi[:,e]
        # Solute name for experiment e
        solute_name = Experiments.loc[e,'Solute (1)']
        # Solvent name for experiment e
        solvent_name = Experiments.loc[e,'Solvent (2)']
        structural_information = structInfo.loc[[solute_name,solvent_name],:]
    
        ## COMBINATORIAL PART
        # The van der Waals properties of the two compounds can be calculated with the
        # help of the van der Waals properties of the groups
        r = np.zeros(i_length) 
        q = np.zeros(i_length) 
        # Loop over every component of every group
        for k in range(k_length): # for every group
            # van der Waals volume
            si = np.array(structural_information.iloc[:,k])
            r += si*vdW_Properties.iloc[k,0]
            # van der Waals surface
            q += si*vdW_Properties.iloc[k,1]
            
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
        combinatorial_e = 1-V+np.log(V+eps)-5*q[0]*(1-(V/(F +eps))+np.log(V/(F+eps)+eps))
        ## RESIDUAL PART
        # MIX
        # The following group mole fractions and sruface area fractions are obtained
        # for the considered binary system at given mole fractions
        X_num = np.zeros(m_length)
        X = np.zeros(m_length)
        # Loop over every component and every group
        for m in range(m_length):
            # numerator of X
            si = np.array(structural_information.iloc[:,m])
            X_num[m] += np.dot(si,mole_fraction)
            # mole fraction of group m in mix
            X = X_num/(sum(X_num)+eps)
    
        Theta_num = np.zeros(m_length)
        Theta_e = np.zeros(m_length)
        # Loop over every group
        for m in range(m_length):
            # numerator
            Theta_num[m] += vdW_Properties.iloc[m,1]*X[m]
            # surface fraction of group m in mix
        Theta_e = Theta_num/(sum(Theta_num)+eps)
        
        # PURE COMPOUNDS 
        # group mole fractions
        X_i = np.zeros((m_length))
        for k in range(k_length): # for every group
            # mole fraction of group m and pure component i
            X_i[k] = structural_information.iloc[0,k]/(sum(structural_information.iloc[0,:])+eps)
    
        # surface area fractions
        Theta_i_num = np.zeros(m_length)
        Theta_i_e = np.zeros(m_length)
        Theta_i_num += vdW_Properties.iloc[:,1]*X_i
        Theta_i_sum = sum(Theta_i_num.T)
        Theta_i_e = Theta_i_num/Theta_i_sum

        combinatorial_e = pd.Series(combinatorial_e)
        Theta_e = pd.Series(Theta_e)
        Theta_i_e = pd.DataFrame(Theta_i_e)
    
        if e == 0:
            combinatorial = combinatorial_e
            Theta = Theta_e
            Theta_i = Theta_i_e
        else:
            combinatorial = pd.concat([combinatorial,combinatorial_e])
            Theta = pd.concat([Theta,Theta_e],axis=1)
            Theta_i = pd.concat([Theta_i,Theta_i_e],axis=1)
        
    UNIFAC_const_out = [combinatorial]
    UNIFAC_const_out.append(Theta)
    UNIFAC_const_out.append(Theta_i)
    return UNIFAC_const_out

# %% UNIFAC FIT
'''
This function calculates the second part of UNIFAC. This part is dependent from
the Interaction Parameter and, thus, needs to be considered by the objective func-
tion.
Input:
- interaction_param: either given from UNIFAC or values that are to be fitted
- Data:
    - structural information (UNIFAC groups) for all solvents and solute
- Experiments: Solute and Solventnames
- Database:
    - van der Waals properties for each Subgroup
    - group to main Group translation for activation coefficients

Output:
- gamma: Activity coefficients of all species in the mixture, [-]
'''
def UNIFAC_fit(interaction_param_fit, structural_information, Experiments, Database, 
               interaction_param, interaction_param_where, UNIFAC_const_out):
    
    eps = np.finfo(float).eps
    # Allocate interaction param
    interaction_param[interaction_param_where[0]] = interaction_param_fit.reshape(-1,1)
    # Reshape Interaction Param
    dim = int(np.sqrt(interaction_param.size))
    interaction_param = np.reshape(interaction_param,(dim, dim))
    interaction_param = pd.DataFrame(interaction_param)   
    
    # Unpack Input
    structInfo = structural_information.copy()
    Experiments = Experiments.copy()
    vdW_Properties = Database[0].copy()
    group_to_mainGroup = Database[1].copy() 
    
    # get Temperature    
    T = Experiments.loc[:,'Temperature']
    
    # length for loops
    k_length = len((structInfo.columns))
    m_length = len(group_to_mainGroup)
    n_length = len(group_to_mainGroup)
    e_length = len(Experiments)
        
    # Loop over all experiments
    for e in range(e_length):
        # Temperature for experiment e
        temperature = T[e]
        # Solute and Solvent for experiment e
        solute_name = Experiments.loc[e,'Solute (1)']
        solvent_name = Experiments.loc[e,'Solvent (2)']
        # Structural information for experiment e
        structural_information = structInfo.loc[[solute_name,solvent_name],:]
        # UNIFAC const output for experiment e
        combinatorial = UNIFAC_const_out[0].iloc[e]
        Theta = UNIFAC_const_out[1].iloc[:,e]
        Theta_i = UNIFAC_const_out[2].iloc[:,e]
        
        ## MIX
        # temperature depsendence of the activity coefficients
        interaction_coeff_temp = np.zeros((m_length, n_length))
        # Loop over all components
        for m in range(m_length):
            for n in range (n_length):
                # translate unique group ID to main group ID
                # find interaction parameter for row m and column n
                group1 = group_to_mainGroup.iloc[m,1] # adjust for list indexing starting at 0
                group2 = group_to_mainGroup.iloc[n,1]
                if group1 == group2:
                    # no interactions within the same group
                    # find interaction parameter for row m and column n
                    interaction_param.iloc[group1,group2] = 0.0
                a_mn = interaction_param.iloc[group1,group2]
                interaction_coeff_temp[m,n] = np.exp(-1*a_mn/(temperature+eps))
                
                      
        # Now all data are avalable to calculate the group activity coefficients in
        # the binary system
        sum1 = np.zeros(m_length)
        # Loop over all components
        for m in range(m_length):
            Theta = Theta.replace(np.nan,0)
            interaction_coeff_temp = np.nan_to_num(interaction_coeff_temp)
            sum1[m] += np.dot(Theta,interaction_coeff_temp[:,m])
            
        sum2 = np.zeros(k_length)
        sum3 = np.zeros(k_length)
        # Loop over all components
        for k in range(k_length):
                sum2[k] += np.dot(Theta,interaction_coeff_temp[:,k])    
                # Loop over all components
                for m in range(m_length):
                    sum3[k] += [(Theta[m]*interaction_coeff_temp[k,m])
                                /(sum1[m]+eps)]       
                    
        group_activity_coeff = np.zeros((len(structural_information.columns)))
        # Group activity coefficients
        vdW_Prop = np.array(vdW_Properties.iloc[:,1])
        group_activity_coeff = vdW_Prop*(1-np.log(sum2+eps)-sum3) 
                
        # group activity coefficient
        sum1_i = np.zeros(len(group_to_mainGroup))
        sum2_i = np.zeros(len(structural_information.columns))
        sum3_i = np.zeros(len(structural_information.columns))

        # Loop over all components
        for m in range(m_length):
            sum1_i += Theta_i.iloc[m]*interaction_coeff_temp[m,:]
        # Loop over all components            
        for n in range(n_length):
            sum2_i += Theta_i.iloc[n]*interaction_coeff_temp[n,:]
        # Loop over all components       
        for k in range(k_length):
            for m in range (m_length):
                sum3_i[k] += [(Theta_i.iloc[m]*interaction_coeff_temp[k,m])/(sum2_i[m]+eps)]
               
        group_activity_coeff_i = np.zeros(len(structural_information.columns))
        # Loop over all components
        for k in range(k_length):
            # group activity coeff
            group_activity_coeff_i[k] = vdW_Properties.iloc[k,1]*(1-np.log(sum1_i[k])-sum3_i[k]+eps)
               
        group_activity_coeff_i = 0
        group_activity_coeff_i = vdW_Properties.iloc[:,1]*(1-np.log(sum1_i+eps)-sum3_i)
                
        # Herewith all data are available to calculate the residual part of the
        # activity coefficients following the solutions of group concepst and finally
        # to calculate the required activity coefficient
        # residual part
        residual = 0
        residual += np.dot(structural_information.iloc[0,:],(group_activity_coeff-group_activity_coeff_i))
        ### ACTIVITY COEFFICIENT
        # activity coefficient
        gamma_e = np.exp(combinatorial+residual)
        gamma_e = pd.Series(gamma_e)
        
        if e == 0:
            gamma = gamma_e
        else:
            gamma = pd.concat([gamma,gamma_e],axis=1)
    return gamma

#%% OBJECTIVE FUNCTION
'''
The target function calls the part of UNIFAC that depends on the interaction
parameters. With the given interaction parameters the activity coecient is cal-
culated. The result is compared with an experimentally determined activity coef-
cient. For this purpose, the error square sum is determined. If this function is
called by the optimizer, it tries to 
nd the interaction parameters so that the sum
of squares of errors is minimal.
Input:
    - Interaction Parameter (Parameter to be fitted)
    - Data, Experiments, Database, Output UNIFAC const, name of the group to pass
    to UNIFAC fit
    - experimental gamma
Output:
    - Residual Sum of Square (Sum over the squared error of all experiments)
'''
def objective(interaction_param_fit, structural_information, Experiments, Database, 
              interaction_param, interaction_param_where, UNIFAC_const_out, gamma_exp):
    # Get gammas for interaction params that are to be fitted
    gamma_mod = UNIFAC_fit(interaction_param_fit, structural_information, Experiments, Database,
                           interaction_param, interaction_param_where, UNIFAC_const_out)
    gamma_mod = np.array(gamma_mod.iloc[0,:])
    No_exp = len(gamma_mod)    # get Error for every experiment
    Err = ((gamma_mod-gamma_exp)/gamma_exp)**2
    # Sum of Squares
    RSS = np.sum(Err)/No_exp
    return RSS

# %% OPTIMIZATION
'''
Set Up of the Optimization Problem:
Initial values and boundary conditions that are given to the optimizer.
Input:
    - Interaction Parameter (Parameter to be fitted)
    - Data, Experiments, Database, Output UNIFAC const, name of the group to pass
    to UNIFAC fit
    - experimental gamma
Output:
    - initial Residual Sum of Square (Sum over the squared error of all experiments)
    - Solution: Output of the optimizer including fitted interaction parameters,
    final Residual Sum of Squares and stats from the optimizer
'''
def optimization(structural_information, Experiments, Database, UNIFAC_const_out, gamma_exp, 
                 interaction_param, interaction_param_where, interaction_param_fit_0, bnds,
                 solver):    
    # RUN OPTIMIZATION
    args = (structural_information, Experiments, Database, interaction_param, interaction_param_where,
            UNIFAC_const_out, gamma_exp)
    func = objective
    if solver == 'L-BFGS-B':
        ## MINIMIZE PARALLEL
        solution = minimize_parallel(func, interaction_param_fit_0, args = args, bounds=bnds,
                                     options={'maxfun': 15000, 'maxiter': 15000},
                                     parallel={'loginfo': True, 'time': True})
    if solver == 'Simplex':
        ## SIMPLEX
        solution =  minimize(func, interaction_param_fit_0, method='nelder-mead', args = args,
                             options={'maxiter': 10000, 'maxfev': 10000, 'disp': True})
    if solver == 'DualAnnealing':
        ## DUAL ANNEALING
        solution = dual_annealing(func, bounds=bnds, args = args, x0 = interaction_param_fit_0)
    return solution

# %% EVALUTATION (Statistics and Plot)
'''
This function evaluates the quality of the model gammas in comparison with the
experimental gammas
Input:
    - gamma_mod: gammas obtaines by UNIFAC
    - gamma_exp: gammas obtained by experiments
Return:
    - Stats: Mean and Std derivation
'''
def evaluation(gamma, fit_extrapolation,Experiments):
    
    gamma_opt = gamma[0]
    gamma_exp = gamma[1]
    gamma_ini = gamma[2]
    
    gamma_opt = gamma_opt.flatten()
    gamma_ini = gamma_ini.flatten()
    # PARITY PLOT
    plt.plot(gamma_opt,gamma_exp,'r*',label='gamma fit') # x vs y
    plt.plot(gamma_ini,gamma_exp,'b*', label = 'gamma initial')
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

    # Squared Error
    Err_ini = ((gamma_ini-gamma_exp)/gamma_exp)**2
    Experiments['Squared Error initial'] = Err_ini
    Err_opt = ((gamma_opt-gamma_exp)/gamma_exp)**2
    Experiments['Squared Error Fit'] = Err_opt
    RSS = [sum(Err_ini), sum(Err_opt)]
    
    # Output File
    Experiments['exp activity coef'] = gamma_exp
    Experiments['initial activity coef'] = gamma_ini
    Experiments['optimized activity coef'] = gamma_opt
    Experiments_with_Error = Experiments
    
    
    # Squared Error over Temperature
    labels = Experiments['Temperature']
    x = np.arange(len(labels))  # the label locations
    width = 0.2
    plt.bar(x - width/2,Err_opt, width, label='gamma opt')
    plt.bar(x + width/2,Err_ini, width, label='gamma initial')
    plt.ylabel('weighted squared Error')
    title = "Squared Error over Temperature - "  + fit_extrapolation
    plt.title(title)
    plt.xlabel('Temperature')
    axes = plt.gca()
    axes.set_ylim([0, 2])
    plt.legend()
    plt.show()
    
    # Regression Results
    comparison = np.sort(np.log10(Err_opt/Err_ini))
    percentage = np.arange(len(comparison))/len(comparison)*100
    plt.plot(percentage, 0*comparison, 'k-')
    plt.plot(percentage, comparison, 'g-')
    plt.xlabel('Percentage of Data Sets')
    plt.ylabel('log10(Squared Error fit / Squared Error initial)')
    title = "Comparison fitted results vs initial results - " + fit_extrapolation
    plt.title(title)
    plt.show()
    
    # Statistics (Mean, Standarddeviation)
    Stats_ini = stats.describe(Err_ini)
    Stats_fit = stats.describe(Err_opt)
    Stats = [Stats_ini, Stats_fit]
    
    # Boxplot
    plt.boxplot((Err_opt,Err_ini))
    plt.title('Deviation of the squared error in a Box Plot for '+fit_extrapolation)
    plt.xticks([1, 2], ['gamma fit - ', 'gamma initial'])
    axes = plt.gca()
    axes.set_ylim([0, 2])
    plt.show()
    
    return Experiments_with_Error, RSS, Stats

# %% EVALUTATION OF ALTERNATIVES (Statistics and Plot)
'''
This function evaluates the quality of the model gammas from interaction parameters
that did not proved to fit superior but where still tried for extrapolation
Input:
    - gamma_mod_alternative: gammas obtaines with alternative interaction coefficients
    - gamma_exp: gammas obtained by experiments
    - RSS_extrapolation: RSS for extrapolation with best fitting interaction parameters
    - sol_LHS: interaction parameters obtained by fitting
Return:
    - Stats: Mean and Std derivation
'''
def evaluation_of_alternatives(gamma_extrapolation, RSS_extrapolation, obj_fit_LHS, sol_fit_LHS):
    gamma_exp = gamma_extrapolation[1]
    gamma_ext_alternative = gamma_extrapolation[3]
    
    RSS_extrapolation_alt = []
    for i in range (len(gamma_ext_alternative)):
        # Squared Error
        gamma_ext_alt = gamma_ext_alternative[i]
        Err_extrapolation_alt = ((gamma_ext_alt-gamma_exp)/gamma_exp)**2
        RSS_ext_alt = Err_extrapolation_alt.sum()
        RSS_extrapolation_alt.append(RSS_ext_alt)
    obj_extrapolation_alt = np.divide(RSS_extrapolation_alt,len(gamma_ext_alternative))
        
    # How well can different optimization results extrapolate the data?
    best_fit_where = np.where(np.round(RSS_extrapolation_alt,2) == round(RSS_extrapolation[1],2))
    plt.plot(np.log10(RSS_extrapolation_alt),'g.',label='alternative fits')
    plt.plot(best_fit_where, np.log10(RSS_extrapolation[1]),'r.',label='best fit')
    plt.xlabel('Runs for different initial values')
    plt.ylabel('Weighted RSS')
    plt.title("weighted RSS over runs for different initial values - extrapolation" )
    axes = plt.gca()
    axes.set_ylim([0, 2])
    plt.legend()
    plt.show()

    # Variation within interaction parameters
    # Sensititvity
    # Normalizing Features
    sol_fit_LHS_norm = StandardScaler().fit_transform(sol_fit_LHS)
    plt.boxplot(sol_fit_LHS_norm)
    plt.xlabel('Interaction Parameter')
    plt.ylabel('normalized Variance')
    title = "Boxplot to show the variation within interaction parameters" 
    plt.title(title)
    plt.show()
    
    # Histogram
    plt.hist(RSS_extrapolation_alt,10)
    plt.xlabel('Sum of weighted squared errors (RSS)')
    plt.ylabel('Count')
    plt.title("Histogram - distribution of the RSS" )
    plt.show()
    
    # Relative Fit
    plt.plot(obj_fit_LHS, obj_extrapolation_alt, 'r.') # x vs y
    plt.xlabel('Fit')
    plt.ylabel('Extrapolation')
    plt.title('RSS-fit vs. RSS-extrapolation')
    axes = plt.gca()
    axes.set_ylim([0, 2])
    plt.show()
    return
# %% RUN ALL
'''
This function runs all functions in the right order. It is run for every group
Input:  - name group, statistics: bool (either yes or no)
Output: - tbd
'''

def run_all(name_group, statistics):
    ## Load Data
    # Data = [Solutes,Solvents,Experiments,structural_information]
    structural_information = load_Data(name_group)
    
    ## Load Experiments
    Experiments = load_Experiments(name_group,'fit')
    gamma_fit_exp = Experiments.loc[:,'gamma_1']
    gamma_fit_exp = np.array(gamma_fit_exp)
    
    ## Load Database
    # Database = [vdW_Properties, group_to_mainGroup, interaction_param]
    Database = load_Database(structural_information, name_group)
    
    # Get UNIFAC interaction parameters (for comparison or initial value)
    interaction_param = Database[2]
    interaction_param = interaction_param.values.flatten()
    # disassemble interaction parameters into UNIFAC and unknowns
    interaction_param_where = disassemble(interaction_param)
    
    ## Run constant UNIFAC
    UNIFAC_const_out = UNIFAC_const(structural_information, Experiments, Database)
    
    ## Set Up optimiation Problem
    # Bounds on variables
    bnds_new = [(-2500, 10000)]
    #bnds_new = [(-25000, 100000)]
    bnds_new = pd.Series(bnds_new)
    for i in range (int(len(interaction_param_where[0]))):
        if i == 0:
            bnds = bnds_new
        else:
            bnds = pd.concat([bnds,bnds_new])
    bnds = tuple(bnds)
    # Initial interaction parameters
    xlimits = np.array(bnds)
    # Latinhypercube Sampling; near random values within boundaries
    sampling = LHS(xlimits=xlimits)
    # Number of runs
    num = 10
    interaction_params =  sampling(num)
    
    ## Optimization
    # initial value of obj funciton
    obj = 1e10
    obj_fit_LHS = []
    sol_fit_LHS = []
    cut_value = len(interaction_param_where[0])*0.1
    
    for i in range (len(interaction_params)):
        interaction_param_fit_0 = interaction_params[i,:]
        # for different starting values the L-BFGS-B is used
        solution = optimization(structural_information, Experiments, Database, 
                                UNIFAC_const_out, gamma_fit_exp,interaction_param, 
                                interaction_param_where,
                                interaction_param_fit_0, bnds, 'L-BFGS-B')
        print('Optimization', i+1, 'out of', num, 'done')
        # save results
        obj_fit_LHS.append(solution.fun)
        sol_fit_LHS.append(solution.x)
        # accept solution if better than the best yet existing
        if solution.fun < cut_value:
            interaction_param_fit_Simplex = solution.x
            # the result which is able to fit the data best is used as input for optimization
            # via simplex
            solution = optimization(structural_information, Experiments, Database, 
                            UNIFAC_const_out, gamma_fit_exp,interaction_param, 
                            interaction_param_where,
                            interaction_param_fit_Simplex, bnds, 'Simplex')
            # save results
            obj_fit_LHS.append(solution.fun)
            sol_fit_LHS.append(solution.x)
            # accept solution if better than the best yet existing
            
        if solution.fun < obj:
            obj = solution.fun
            solution_save = solution
            print("new solution accepted")
    #save best fitting interaction parameters
    solution = solution_save
    # Save results of all LHSations
    obj_fit_LHS = np.array(obj_fit_LHS)
    sol_fit_LHS = np.array(sol_fit_LHS)
    # kick results where the optimization failed
    Filter_obj_LHS = np.array(np.where(obj_fit_LHS>cut_value))
    obj_fit_LHS = np.delete(obj_fit_LHS,Filter_obj_LHS,0)
    sol_fit_LHS = np.delete(sol_fit_LHS,Filter_obj_LHS,0)
    results_fit_LHS = [obj_fit_LHS, sol_fit_LHS]
    # Fitted interaction parameters
    interaction_param_opt = solution.x
    
    ## Run variable part of UNIFAC to get gammas for the fit
    # optimized interaction paramters
    gamma_fit_opt = UNIFAC_fit(interaction_param_opt, structural_information, Experiments, Database, 
                           interaction_param, interaction_param_where, UNIFAC_const_out)
    gamma_fit_opt = np.array(gamma_fit_opt)
    interaction_param_ini = np.zeros(len(interaction_param_where[0]))
    #interaction_param_ini = np.array([476.4, 26.76]) 
    # initial interaction parameters
    gamma_fit_ini = UNIFAC_fit(interaction_param_ini, structural_information, Experiments, Database, 
                              interaction_param, interaction_param_where, UNIFAC_const_out)
    gamma_fit_ini = np.array(gamma_fit_ini)
    
    ## evaluate the results - if wanted
    if statistics:
        gamma_fit = [gamma_fit_opt, gamma_fit_exp, gamma_fit_ini]
        Experiments_Error_fit, RSS_fit, Stats_fit = evaluation(
        gamma_fit, 'fit', Experiments)
    else:
        Experiments_Error_fit = {}
        RSS_fit = {}
        Stats_fit = {}
        
        
    ## EXTRAPOLATION
    # Run UNIFAC fit with other molecules and optimized interaction parameters
    Experiments = load_Experiments(name_group,'extrapolation')
    UNIFAC_const_out = UNIFAC_const(structural_information, Experiments, Database)
    
    # get activity coefficients
    # activity coefficients with optimized
    gamma_ext_opt = UNIFAC_fit(interaction_param_opt, structural_information, Experiments, Database, 
                   interaction_param, interaction_param_where, UNIFAC_const_out)
    gamma_ext_opt = np.array(gamma_ext_opt)
    # Run alternative optimization Results to compare them for extrapolation
    gamma_ext_alternative = []
    for i in range (len(sol_fit_LHS)):
        interaction_param_alternative = sol_fit_LHS[i,:]
        gamma_ext_alt = UNIFAC_fit(interaction_param_alternative, structural_information, Experiments, 
        Database, interaction_param, interaction_param_where, UNIFAC_const_out)
        gamma_ext_alt = np.array(gamma_ext_alt)
        gamma_ext_alternative.append(gamma_ext_alt)
    # activity coefficient for when unknown interaction params = 0
    gamma_ext_ini = UNIFAC_fit(interaction_param_ini, structural_information, Experiments, Database, 
                      interaction_param, interaction_param_where, UNIFAC_const_out)
    gamma_ext_ini = np.array(gamma_ext_ini)
    # experimentally obtained activity coefficients
    gamma_ext_exp = Experiments.loc[:,'gamma_1']
    gamma_ext_exp = np.array(gamma_ext_exp)
    
    # Statistic
    if statistics:
        gamma_extrapolation = [gamma_ext_opt, gamma_ext_exp, gamma_ext_ini, gamma_ext_alternative]
        Experiments_Error_extrapolation, RSS_extrapolation, Stats_extrapolation = evaluation(
        gamma_extrapolation, 'extrapolation', Experiments)
        evaluation_of_alternatives(gamma_extrapolation, RSS_extrapolation, obj_fit_LHS, sol_fit_LHS)
    else:
        Experiments_Error_extrapolation = {}
        RSS_extrapolation = {}
        Stats_extrapolation = {}
                
    dim = int(np.sqrt(interaction_param.size))
    interaction_param[interaction_param_where[0]] = interaction_param_opt.reshape(-1,1)
    interaction_param = np.reshape(interaction_param,(dim, dim))
    interaction_param = pd.DataFrame(interaction_param)  
    Experiments_with_Error = [Experiments_Error_fit,Experiments_Error_extrapolation]
    RSS = [RSS_fit,RSS_extrapolation]
    RSS = pd.DataFrame(RSS)
    RSS.columns = ['initial params', 'optimized params']
    RSS.index = ['Fit', 'Extrapolation']
    Stats = [Stats_fit,Stats_extrapolation]
    return interaction_param, Experiments_with_Error, RSS, Stats, results_fit_LHS
# %% MAIN
# If called as main programm
if __name__ == "__main__":
    # different groups
    Groups = 'Validierung'
    # Groups.append('OLEDs')
    # small number
    eps = np.finfo(float).eps
    #for name_group in Groups:
    name_group = Groups
    # Run program
    interaction_param, Experiments_with_Error, RSS, Stats, results_fit_LHS  = run_all(name_group,True)
    print('Hannes, kannst du die Bilder abspeichern?')