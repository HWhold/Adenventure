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