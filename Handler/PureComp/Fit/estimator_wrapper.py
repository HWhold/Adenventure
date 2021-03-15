# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 18:33:47 2021
This function contains the wrapper for the estimators
1) Vapor pressure
2) Discrete properties

TODO test needed
@author: X202722
"""

import pandas as pd
from Adenventure.PhysModels.PureComp.testBPT import useOffset
from Adenventure.PhysModels.PureComp.testMPT2 import useOffset_MPT 
from Adenventure.PhysModels.PureComp.testHFUS import useOffset_HFUS
from Adenventure.Handler.PureComp.Fit.load_prereq2 import load_prereq2
from Adenventure.PhysModels.PureComp.testVAP import ParaMollerRarey, ParaNannoolalRarey, PresMollerRarey, PresNannoolalRarey

def estimateValue(fileName, properti, method, pathToFile, offset, na, M):
    # print(fileName)
    # this function calculates values for fitted estimators
    # filename: name of file or molecule in ctc file for vap1m module [String]
    # properti: estimated porperti [integer]
    # method: used method for fitting [integer]
    # path to file: path to the ddbst-result files (comes out of vap1m module) [String]
    # offset:calculated offset in parameter search room  [float]
    value = 0
    if properti == 'BPT':
        _, _, rest_val =  load_prereq2 (fileName, pathToFile)
        T_sim2 = rest_val.iloc[method*3,4]
        value = useOffset (offset, T_sim2, method, na, M)
    elif properti == 'MPT':
        _, _, rest_val =  load_prereq2 (fileName, pathToFile)
        T_sim2 = rest_val.iloc[method,5]
        value = useOffset_MPT (offset, T_sim2, method, na, M)
    elif properti == 'HFUS':
        _, _, rest_val =  load_prereq2 (fileName, pathToFile)
        T_sim2 = rest_val.iloc[method,6]
        value = useOffset_HFUS (offset, T_sim2, method, na, M)
    else:
        value = 0        
    return value

def estimateValue_VAP(fileName, method, pathToFile, offset, Tb_opt, T_calc):
    # print(fileName)
    # this function calculates values for fitted estimators
    # filename filename of evaluated compound
    # method method as int
    # pathtofile path to file where DDBST data lies
    # offset calculated offset
    # Tb_opt optimized Boiling point
    # T_calc temperature of evaluation
    VAPopt = 0
    _, vap_val, rest_val =  load_prereq2 (fileName, pathToFile)
    p_sim = vap_val.iloc[:,method+2]    
    Tb_sim = rest_val.iloc[0,4]
    # 345 because 330 is taken as header
    T_sim = [345+i*15 for i  in range(14) ]
    if method == 0:
        #MollerMethod
        Bs, Ds = ParaMollerRarey(p_sim[1], p_sim[3], Tb_sim, T_sim[1], T_sim[3])
        VAPopt = PresMollerRarey(Bs+offset[0],Ds+offset[1],Tb_opt,T_calc)
        
    if method == 1:
        #Nannoolal Method
        Ds_sim= ParaNannoolalRarey(p_sim[1], Tb_sim, T_sim[1])
        VAPopt = PresNannoolalRarey(Ds_sim+offset,Tb_opt,T_calc)
        
    return VAPopt

def estimator_wrapper_main(line_ends,exp_data2, process_data, pathToFile,VAP_T):
    # line ends defines the structure of the process data dataframe where which property ends: [list]
    # exp_data2 is the offset data containing best fit method as well as all other fits [pd.DataFrame]
    # process data is the data which is supposed to be fit. contains IDs, fit molecule, structure, filename and occurence [pd.DataFrame]
    # path to file gives path to where all files of vap1m module lay
    # process_data contains all fitted data
       ####### catch fractions of area% equaling zero

    process_data['FlächenAnteil in Prozent'][process_data['FlächenAnteil in Prozent']<1e-10] = 1e-2

    properti = 'BPT'
    [BPT_offset_start, BPT_offset_end, MPT_offset_start, MPT_offset_end, HFUS_offset_start, HFUS_offset_end, VAP_offset_start, VAP_offset_end] = line_ends
    # reads best fit methods for attached molecules __BPT__!
    process_data['BPT_FIT_method'] = process_data['Zuordnung Fit'].apply(lambda x:exp_data2['BPT Best estimator'].iloc[x-1] ).copy()
    # this makes the fit : estimateValue(fileName, properti, method, pathToFile, offset, na, M)
    process_data['BPT_FIT'] = process_data.apply(lambda x: estimateValue(str(x.iloc[4]), properti, 0, pathToFile, exp_data2.iloc[x.iloc[11]-1,int(x.iloc[14])+BPT_offset_start], x.iloc[13], x.iloc[12]), axis = 1).copy()
    
    # reads best fit methods for attached molecules __MPT__!
    properti = 'MPT'
    process_data['MPT_FIT_method'] = process_data['Zuordnung Fit'].apply(lambda x:exp_data2['MPT Best estimator'].iloc[x-1] ).copy()
    # this makes the fit : estimateValue(fileName, properti, method, pathToFile, offset, na, M)
    process_data['MPT_FIT'] = process_data.apply(lambda x: estimateValue(str(x.iloc[4]), properti, int(x.iloc[16]), pathToFile, exp_data2.iloc[x.iloc[11]-1,int(x.iloc[16])+MPT_offset_start], x.iloc[13], x.iloc[12]), axis = 1).copy()
    
    # reads best fit methods for attached molecules __HFUS__!
    properti = 'HFUS'
    process_data['HFUS_FIT_method'] = process_data['Zuordnung Fit'].apply(lambda x:exp_data2['HFUS Best estimator'].iloc[x-1]).copy()
    # this makes the fit : estimateValue(fileName, properti, method, pathToFile, offset, na, M)
    process_data['HFUS_FIT'] = process_data.apply(lambda x: estimateValue(str(x.iloc[4]), properti, int(x['HFUS_FIT_method']), pathToFile, exp_data2.iloc[x.iloc[11]-1,int(x['HFUS_FIT_method'])+HFUS_offset_start], x.iloc[13], x.iloc[12]), axis = 1).copy()
    
    
    # reads best fit methods for attached molecules __VAP__!
    properti = 'VAP'
    process_data['VAP_FIT_method'] = process_data['Zuordnung Fit'].apply(lambda x:exp_data2['VAP Best estimator'].iloc[x-1] ).copy()
    # this makes the fit : estimateValue_VAP(fileName, method, pathToFile, offset, Tb_opt, T_calc)
    process_data['VAP_FIT'] = process_data.apply(lambda x: estimateValue_VAP(str(x.iloc[4]), int(x['VAP_FIT_method']), pathToFile, exp_data2.iloc[x.iloc[11]-1,int(x['VAP_FIT_method'])+VAP_offset_start], x['BPT_FIT'], VAP_T), axis = 1).copy()
    
    #tada we have fitted all data
    return process_data