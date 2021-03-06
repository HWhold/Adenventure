# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 10:35:23 2020

@author: X202722
"""

import itertools
import functools
import pandas as pd
import numpy as np
from runVAPS_rev5 import parameters
from fitVapToExpValuesTest import clapeyron, clapeyronFit
# from samplingCoefficients_fitall import samplingCoefficient

def make_parameter_VAP_fit (Tb_sim, Tb_exp,T_sim, p_sim, p_exp,T_exp, method):
    import numpy as np
    # this function returns offsets for moller and nannoolal
    # TODO Schwerpunkt finden
    # fit if possible BPT guess to 
    OffsetP = 0
    # creates parameters based on fit withexperimental data
    para = parameters()    
    
    # use fit from antoine equation to average experimental errors
    # p3 = antoine(antoine_result.x,T)
    # using clapeyron to avoid falling into antoine optimization fuck up problems
    para2 = clapeyronFit(T_exp, p_exp)
    p3 = clapeyron (para2,T_sim)    
    

    if method == 0:
        #Moller
        #get parameters from GCT via explicit solution, important to use simulation values
        # so gct values can be reverse engineered
        para.Moller_Bs, para.Moller_Ds = ParaMollerRarey(p_sim[1], p_sim[3], Tb_sim, T_sim[1], T_sim[3])
        # iteration throws out bullshit when p_sample > p_atmosphere   
        T2=T_exp[p_exp<1000]
        p4=p_exp[p_exp<1000]
        # print(p3, para2)
        #calculate parameters for models with fitted experimental values and BPT from DDBST
        # para.Bs3,para.Ds3 = ParaMollerRarey(p4[0], p4[-1],Tb_exp, T2[0], T2[-1])
        para.Bs3,para.Ds3 = samplingCoefficient(T2,p4, Tb_exp, para,T_exp,p_exp, 'Moller')
        OffsetP = [para.Bs3-para.Moller_Bs,para.Ds3-para.Moller_Ds]
        p_sim_scal = PresMollerRarey(para.Moller_Bs,para.Moller_Ds,Tb_sim,T_exp)
    elif method == 1:
        #Nannoolal
        # get GCT values from DDBST
        Ds_sim= ParaNannoolalRarey(p_sim[1], Tb_sim, T_sim[1])
        # create true values from experimental data
        # Ds_exp= ParaNannoolalRarey(p_exp[1], Tb_exp, T_sim[1]) 
        _, Ds_exp = samplingCoefficient(T_exp, p_exp, Tb_exp, para,T_exp,p_exp, 'Nannolal')
        #calculate offset
        OffsetP= Ds_exp-Ds_sim    
        p_sim_scal = PresNannoolalRarey(Ds_sim, Tb_sim, T_exp)  
        
    return OffsetP, p_sim_scal


def samplingCoefficient(T2,p4, Tb, para,T,p, model):
    #loop through all data points from experiment and check whether the fit works or not
    # take closest fit to boiling point
    
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


def useOffset_VAP(OffsetP, p_sim, Tb_sim, T_sim,Tb_opt, T_calc, method):
    # this function calculates optimized vapor pressure
    # First get GCT values then calculate needed P
    if method == 0:
        #MollerMethod
        Bs, Ds = ParaMollerRarey(p_sim[1], p_sim[3], Tb_sim, T_sim[1], T_sim[3])
        VAPopt = PresMollerRarey(Bs+OffsetP[0],Ds+OffsetP[1],Tb_opt,T_calc)
        
    if method == 1:
        #Nannoolal Method
        Ds_sim= ParaNannoolalRarey(p_sim[1], Tb_sim, T_sim[1])
        VAPopt = PresNannoolalRarey(Ds_sim+OffsetP,Tb_opt,T_calc)
        
    return VAPopt

def ParaMollerRarey(p1,p2,Tb,T1,T2):
    import numpy as np
    # calculate GI parameters from RareyMoller pressure is taken in mbar
    C =  -2.65+np.power(Tb,1.485)/135
    # print(C)
    S1 = (T1-Tb)/(T1-C)
    S2 = (T2-Tb)/(T2-C)
    
    F1 = np.log(T1/Tb)
    F2 = np.log(T2/Tb)
    
    lp1 = np.log(p1/1013.25)
    lp2 = np.log(p2/1013.25)
   
    Ds = (lp1/S1-lp2/S2)/(F1/S1-F2/S2)
    Bs = (lp1-Ds*F1)/S1
    
    # print('check',PresMollerRarey(Bs, Ds, Tb, T2)/p2)
    # print('check',PresMollerRarey(Bs, Ds, Tb, T1)/p1)
    return Bs,Ds

def equationsMoller(k,p,T):
    import numpy as np
    T_b  =  k[0] 
    B =   k[1] 
    D  = k[2]
    # print(np.log(T[0]/T_b)-p[0])
    return [(B*(T[0]-T_b))/(T[0]-(-2.65+np.power(T_b,1.435)/135))+D*np.log(T[0]/T_b)-p[0],
            (B*(T[1]-T_b))/(T[1]-(-2.65+np.power(T_b,1.435)/135))+D*np.log(T[1]/T_b )-p[1],
            (B*(T[2]-T_b))/(T[2]-(-2.65+np.power(T_b,1.435)/135))+D*np.log(T[2]/T_b )-p[2]]

def ParaMollerRarey_Tb(p1_in,p2_in,p3_in,Tb1,T1,T2,T3,Bs, Ds):
    # calculate GI parameters from RareyMoller without boiling point 
    # pressure is taken in mbar
    import numpy as np
    Tb =Tb1
    p1 = np.log(p1_in/1013.25)
    p2 = np.log(p2_in/1013.25)
    p3 = np.log(p3_in/1013.25)

    from scipy.optimize import fsolve, broyden1
    import math
    import numpy as np
    from functools import partial
    
    p1 = np.log(p1_in/1013.25)
    p2 = np.log(p2_in/1013.25)
    p3 = np.log(p3_in/1013.25)
    T_iter = 0
    dTb = 1
    i = 1
    while abs(dTb) > 1e-7:
        #reduce equation size through substitution
        s1 = (T1-Tb)/(T1+2.65-np.power(Tb,1.485)/135)
        s2 = (T2-Tb)/(T2+2.65-np.power(Tb,1.485)/135)
        s3 = (T3-Tb)/(T3+2.65-np.power(Tb,1.485)/135)
       
        f1 = np.log(T1/Tb)
        f2 = np.log(T2/Tb)
        f3 = np.log(T3/Tb)
        
        a =( (p1/f1-p2/f2)/(s1/f1-s2/f2)- (p1/f1)/(s1/f1-s3/f3))*(s1/f1-s3/f3)
        
        T_iter = T3/np.exp(-p3/a)

        dTb = T_iter-Tb
        Tb=T_iter
        i = i+1
        if i> 10000:
            break

    Tb2 = Tb
    Bs2,Ds2 = ParaMollerRarey(p1_in, p2_in, Tb2, T1, T2)     
    
    return Bs2, Ds2, Tb2


def PresMollerRarey(Bs,Ds,Tb,T):
    # calculate vapor pressure with GI parameters from moller and rarey
    #return value is in mbar
    # p in mbar T in K
    import numpy as np
    C =  -2.65+(np.power(Tb,1.485))/135
    
    S = (T-Tb)/(T-C)
    
    F =  np.log(T/Tb)
    
    rightSide = Bs*S+Ds*F
    P = np.exp(rightSide)*1013.25
    return P

def ParaNannoolalRarey(p,Tb,T):
 # calculate vapor pressure with GI parameters from nannoolal and rarey
    # p in mbar T in K
    import numpy as np
    Trb = T/Tb
    p_rel = np.log10(p/1013.25)
    K = (Trb-1)/(Trb-1/8)
    dB=p_rel/K-4.1012
    dB=np.log10(p/1013.25)/((Trb-1)/(Trb-1/8))-4.1012
    return dB


def ParaNannoolalRarey_Tb(p1_in,T1,p2_in,T2):
 # calculate vapor pressure with GI parameters from nannoolal and rarey
    # p in mbar T in K
    import numpy as np
    
    p1 = np.log10(p1_in/1013.25)
    p2 = np.log10(p2_in/1013.25)
    
    a = T1*T2*(p1-p2)
    
    b =  -1/8*T2 * p1 + 1/8 * T1* p2  + T2 * p2 - T1* p1
    
    c = 1/8*(p1-p2)  
    
    # print('a b c', a,b,c,(b*b)-4*a*c, (b*b)-4*a*c)
    Tb1 = 1
    Tb2 = 1
    if (b*b)-4*a*c > 0:
        x_1=(-b+np.sqrt((b*b)-4*a*c))/2/a
        Tb1 = 1/x_1
        # print('asdf',Tb)
        # we should check which one is closer to BPT_sim for reality reasons but in general X-2 is better
    # else:
        x_2=(-b-np.sqrt((b*b)-4*a*c))/2/a
        Tb2 = 1/x_2
        # print('asdf',Tb)
    if Tb1<Tb2:
        Tb = Tb1
    else:
        Tb=Tb2    
    return Tb

def PresNannoolalRarey(dB,Tb,T):
    # return value in mbar
    # p in mbar T in K
    import numpy as np
    Trb = T/Tb
    K = (Trb-1)/(Trb-1/8)
    P = 1013.25*np.power(10,(4.1012+dB)*(T/Tb-1)/(T/Tb-1/8))
    return P


if __name__ == "__main__":
    from helfFvalidation import import_Vap_small, vapor_eq_calc
    import matplotlib.pyplot as plt
    import numpy as np
    # everything in °C
    Tb = 159.192
    T1 = 158
    T2 = 100
    
    var1 = -8018.64
    var2 = 25.4741
    
    p1 =  vapor_eq_calc(273+T1,var1,var2)
    p2 =  vapor_eq_calc(273+T2,var1,var2)
    print(p1, p2)
    Bs,Ds = ParaMollerRarey(p1,p2, Tb+273, T1+273,T2+273)
    print(Ds,Bs)
    dB = ParaNannoolalRarey(p1 ,Tb, T1+273)
    dB= 0
    save = np.ones(100)
    save2 = np.ones(100)
    Temps = np.ones(100)
    for i in range(100):
        T=i*1+350
        Temps[i]= T
        #save[i]=PresMollerRarey(Bs, Ds, Tb+273, T)
        save[i]=PresNannoolalRarey( -dB,Tb+273, T)
        save2[i] =  vapor_eq_calc(T,var1,var2)
    print(PresMollerRarey(Bs, Ds, Tb+273, 190+273))
    plt.plot(Temps,save)
    plt.plot(Temps,save2)
    #plt.plot(Temps,save/save2)
    plt.show()