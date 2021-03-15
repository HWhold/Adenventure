# -*- coding: utf-8 -*-
"""
Created on Tue May 26 07:49:48 2020

@author: X202722
"""


def make_parameter_BPT_fit (T_sim, T_exp, method, na, M):
    import numpy as np
    # fit if possible BPT guess to 
    #nannoolal
    if method == 0:
        a = 0.6583
        b= 1.6868
        c= 84.3395
        
        value = (T_sim-c)*(np.power(na,a)+b)
        
        OffsetP = (T_exp-c)*(np.power(na,a)+b)-value
        # OffsetP = T_exp-T_sim
    elif method == 1:
        #champion
        OffsetP = T_exp-T_sim
    elif method == 2:
        #stein and brown
        A=-94.84
        B=1.5577
        C=-0.0007705
        
        A2=282.7
        B2= 0.5209
                
        if  (T_sim-A2)/(1-B2) < 700:          
            a = C
            b = B
            c = (A-T_sim)
            # value =  (-b-np.sqrt(np.power(b,2)-4*a*c))/(2*a)
            value = (-b+np.sqrt(np.power(b,2)-4*a*c))/(2*a)
            
            a = A+B*value+C*np.power(value,2)-T_exp
            b = B + 2 * C * value
            c = C
                        
            OffsetP =-1/2*(2*C*value + B + np.sqrt(B**2 - 4*A*C + 4*C*T_exp))/C
            # -1/2*(2*C*value + B + sqrt(B^2 - 4*A*C + 4*C*T_exp))/C, O == -1/2*(2*C*value + B - sqrt(B^2 - 4*A*C + 4*C*T_exp))/C]
            # OffsetP = -(A - (A**2 + 2*A + 4*C*T_exp + 1)**(1/2) + 2*C*value + 1)/(2*C)
            
        else:
            # value = (T_sim-A2)/(1-B2)
            OffsetP = (T_exp-A2)/(1-B2)-T_sim
            
    elif method == 3:
        #devottta
        OffsetP = T_exp-T_sim
    elif method == 4:
        #joback
        OffsetP = T_exp-T_sim
    elif method == 5:
        #GVS
        OffsetP = T_exp-T_sim
    elif method == 6:
        #gani
        k = 204.359
        OffsetP = np.exp(T_exp/k) -np.exp(T_sim/k)
        # OffsetP = T_exp-T_sim
    elif method == 7:
        #marrero
        k = -0.366
        if M == []:
            OffsetP = T_exp-T_sim
        else:
            OffsetP = (T_exp-T_sim)/np.power(M,k)
            
    elif method == 8:
        #marrero simple groups
        OffsetP = T_exp-T_sim
    elif method == 9:
        
        #marrero simple groups, simple approach
        OffsetP = T_exp-T_sim
        
    return OffsetP

def useOffset (OffsetP, T_sim, method, na, M):
    import numpy as np
    # use fit to create new value
    #nannoolal
    if method == 0:
        a = 0.6583
        b= 1.6868
        c= 84.3395
        value = (T_sim-c)*(na**a+b)
        Topt = (value + OffsetP)/(na**a+b)+c
        # Topt = T_sim + OffsetP
    elif method == 1:
        #champion
        Topt = T_sim + OffsetP
    elif method == 2:
        # stein and brown
        
        A=-94.84
        B=0.5577
        C=-0.0007705
        A2=282.7
        B2= 0.5209
        
        bb = (B+1)
        
        if  (T_sim-A2)/(1-B2) < 700:          
            a = C
            b = bb
            c = (A-T_sim)
            # value =  (-b-np.sqrt(np.power(b,2)-4*a*c))/(2*a)
            value2 = (-b+np.sqrt(np.power(b,2)-4*a*c))/(2*a)
            
            Topt = A + bb*(value2+OffsetP) + C * (value2+OffsetP)**2
            
        else:
            value = (T_sim-A2)/(1-B2)
            Topt = (T_sim+OffsetP)*(1-B2)+A2
            
    elif method == 3:
        #devottta
        Topt = T_sim + OffsetP
    elif method == 4:
        #joback
        Topt = T_sim + OffsetP
    elif method == 5:
        #GVS
        Topt = T_sim + OffsetP
    elif method == 6:
        #gani
        k = 204.359
        value = np.exp(T_sim/k)
        if value + OffsetP<0:
            print(value)
            print(OffsetP)
            
        Topt = k * np.log(value + OffsetP)
        # Topt = T_sim + OffsetP
    elif method == 7:
        #marrero
        k = -0.366
        b = 149.84
        
        if M == []:    
            Topt = T_sim + OffsetP
        else:
            value = (T_sim-b)/(M**(k))
            Topt = (value+OffsetP)*M**k+b
        # Topt = T_sim + OffsetP
    elif method == 8:
        #marrero simple groups
        Topt = T_sim + OffsetP
    elif method == 9:
        #marrero simple groups, simple approach
        Topt = T_sim + OffsetP
        
    return Topt

def boiling_point_wrapper (res, meta_real):
    #['T', 'AZ', 'TC', 'PC', 'BPT', 'MPT', 'HFUS']
    #loop through columns and compare available data with simulated ones
    # choose best fit and calculate best fit
    import pandas as pd
    import numpy as np

    #create variables
    
    T_fit = np.zeros((10,len(res)))
    T_fit_list = []
    names = []
    M = []
    
    #run through all interesting molecules
    for p in range(len(res)):
        # 10 methods to look at
        if not np.isnan(meta_real.iloc[p,0]):

            # get experimental value
            T_exp = res[p].iloc[0,0]+273
            
            #grab auxiliary values                
            na = res[p].iloc[0,10]
            M = res[p].iloc[0,11]
            
            #loop through methods
            for i in range(10):
                #get simulated values
                T_sim = res[p].iloc[i*3+1,0]
                method = i
                # check if there is a simulated value
                if np.isnan(T_sim):
                    offset = np.nan
                    T_fit[i,:] =  np.nan
                else:  
                    #if yes calculate offset
                    offset = make_parameter_BPT_fit(T_sim, T_exp, method, na, M)
                    #offset = 0
                    #then use offset to calculate new values
                    for k in range(len(res)):
                        
                        T_sim_fit = res[k].iloc[i*3+1,0]
                        T_exp_fit = res[k].iloc[0,0]+273
                        na_fit = res[k].iloc[0,10]
                        M_fit = res[k].iloc[0,11]
                        
                        #check if offset is nan
                        if np.isnan(T_sim_fit):
                            T_fit[i,k] = np.nan
                        else:    
                            # fit value with calculated offset
                            T_fit[i,k] =  useOffset(offset,T_sim_fit, method, na_fit,M_fit, T_exp_fit)
                        
            T_fit_list.append(pd.DataFrame(T_fit))
            names.append(meta_real.index[p])
            T_fit = np.zeros((10,len(res)))
    
    return T_fit_list, names

def BPT_sim_data(res_clean, meta_real):
    # calculate b
    import numpy as np
    import pandas as pd
    T_fit_list = []
    T_sim = np.zeros(10)
    for p in range(len(res_clean)):
        # 10 methods to look at
        if not np.isnan(meta_real.iloc[p,0]):
            for i in range(10):
                T_sim[i] = res_clean[p].iloc[i*3+1,0]
        T_fit_list.append(T_sim)
        T_sim = np.zeros(10)
    
    return T_fit_list
                
def prepare_Violin_plot(data, title, ax):
    # print violinplot
    import numpy as np
    import matplotlib.pyplot as plt
    data_plot_T_pure = []
    #remove nan for boxplot
    for T in data.T.columns:
        T_sim_flat = data.T[T].to_numpy().flatten()
        T_sim_flat = T_sim_flat[~np.isnan(T_sim_flat)]
        data_plot_T_pure.append(T_sim_flat)
        
    ax.set_title(title)
    ax.boxplot(data_plot_T_pure)
    
    return ax
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        