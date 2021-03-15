import pandas as pd

def load_prereq2 (name_material, path_empty):
    #load simulated data from DDBST from disc.
    #now load simulated values (first contains VAP data)
    #extended version loads big matrix as well
    #path_empty ='C:\\Users\\X202722\\Desktop\\Prozesssynthese\\13 DSC DATEN\\Strukturen_ohneN\\'
    path_to_ctc2 = '.ctc2.txt'
    path = path_empty + name_material + path_to_ctc2
    # contains general modelled data
    path_to_help_data = '.ctchelpDat2.txt'
    # put all together change hard coded name to looped index
    path_to_BPT = path_empty + name_material + path_to_help_data 
    # big matrix which contains permutations of different properties
    # ['T', 'AZ', 'TC', 'PC', 'BPT', 'MPT', 'HFUS']
    path_to_rest_val=".ctcmatrixSize.txt"
    sum_rest_val = path_empty + name_material + path_to_rest_val 
    #load data where path has been specified
    Sim_vap_0 = pd.read_csv(path, delimiter = ' ', header = 0)
    basic_data = pd.read_csv(path_to_BPT, delimiter = ' ', header = None)
    rest_val = pd.read_csv(sum_rest_val, delimiter = ' ', header = None)
    rest_val.columns  = ['T', 'AZ', 'TC', 'PC', 'BPT', 'MPT', 'HFUS']
    return basic_data, Sim_vap_0, rest_val