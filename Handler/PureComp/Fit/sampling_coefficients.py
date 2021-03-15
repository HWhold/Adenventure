import pandas as pd
import itertools
from Adenventure.PhysModels.PureComp.testVAP import *

def samplingCoefficient(T2,p4, Tb, para,T,p, model):
# this function creates reliable fitted GCT data for the GCT methots nannolal and mollerby maximum likelihood estimation with assumed gaussiaun distribution
#loop through all data points(p4,T2) from experiment and creates fits for BPT according to analytical solutions
# takes average of all fitted information and spits it out as result
    
    res_fits = []
    # p2 = p[p<1000]
    # T2 = T[p<1000]
    if model  == 'Moller':
        combi = pd.DataFrame(itertools.combinations(range(len(p4)), 2))
        for a in combi.index:
                # print(a)
                a,b = combi.iloc[a,:]
                # print(ParaMollerRarey_Tb(p4[a], p4[b],p4[c], Tb, T2[a], T2[b], T2[c],para.Moller_Bs, para.Moller_Ds))
                res_fits.append(ParaMollerRarey(p4[a], p4[b], Tb, T2[a], T2[b]))
                # res_fits.append(ParaMollerRarey_Tb(p2[a], p2[b],p2[c], Tb, T2[a], T2[b], T2[c],para.Moller_Bs, para.Moller_Ds))
    elif model == 'Nannolal':
        combi = pd.DataFrame(itertools.combinations(range(len(T2)), 2))
        for a in combi.index:
            # print(a)
            a,b = combi.iloc[a,:]
            # print(ParaMollerRarey_Tb(p4[a], p4[b],p4[c], Tb, T2[a], T2[b], T2[c],para.Moller_Bs, para.Moller_Ds))
            Tb = ParaNannoolalRarey_Tb(p4[a], T2[a], p4[b], T2[b])
            Bs = ParaNannoolalRarey(p4[a], Tb, T2[a])
            res_fits.append([Tb, Bs])
            
    res_fits = pd.DataFrame(res_fits)
    if model == 'Moller':
        # res_fits['diffTb'] = np.abs(res_fits.iloc[:,2]-Tb)
        # res_fits['diffpT'] = res_fits.apply(lambda x:np.sum(np.abs(PresMollerRarey(x.iloc[0],x.iloc[1],x.iloc[2],T2)-p4)/p4),axis = 1)
        # res_fits_fit = res_fits[[0,1]][np.logical_and(res_fits[2]>Tb-200,res_fits[2]<2000)]    
        res_fits_fit = res_fits
    elif model == 'Nannolal':
        res_fits['diffTb'] = np.abs(res_fits.iloc[:,0]-Tb)
        res_fits['diffpT'] = res_fits.apply(lambda x:np.sum(np.abs(PresNannoolalRarey(x.iloc[1],x.iloc[0],T2)-p4)/p4),axis = 1)      
        res_fits_fit = res_fits[[0,1]][np.logical_and(res_fits[0]>Tb-200,res_fits[0]<2000)] 
    # res_fits = res_fits.sort_values(by = 'diffTb', axis = 0).reset_index()
    
      
    return np.mean(res_fits_fit )     
