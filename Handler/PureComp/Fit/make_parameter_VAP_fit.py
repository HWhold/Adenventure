
from Adenventure.PhysModels.PureComp.fitParameterClass import parameters
from Adenventure.PhysModels.Thermodynamics.clapeyronEQ import clapeyronFit, clapeyron
from Adenventure.PhysModels.PureComp.testVAP import (
    ParaMollerRarey, 
    equationsMoller, 
    ParaMollerRarey_Tb, 
    PresMollerRarey, 
    ParaNannoolalRarey, 
    ParaNannoolalRarey_Tb, 
    PresNannoolalRarey
)
from Adenventure.Handler.PureComp.Fit.sampling_coefficients import samplingCoefficient  

def make_parameter_VAP_fit (Tb_sim, Tb_exp,T_sim, p_sim, p_exp,T_exp, method):
    import numpy as np
    # this function returns offsets for moller and nannoolal
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
        #calculate parameters for  with fitted experimental values and BPT from DDBST
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

