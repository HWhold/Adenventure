def CalcBinary(prods, table, db, T):
    # careful with COSMO1 and COSMO2
    import pandas as pd
    import cCOSMO
    import numpy as np
    gamma = []
    numbers = []
    for i in range(len(table)):
        #here names are defined
        name1 = prods
        name2 = table.iloc[i]['CAS#']
        names = [name1, name2]
        # names = ['BENZENE','ETHANOL']
        # print(names)
        #loaded from harddrive and into virtual db
        for iden in names:
            # print(iden)
            n = db.normalize_identifier(iden)
            db.add_profile(n)
        # COSMO = cCOSMO.COSMO3(names, db)
        COSMO = cCOSMO.COSMO3(names, db)
        cConsts = COSMO.get_mutable_COSMO_constants()
        cConsts.fast_Gamma = False
        # here I should implement solvent dependent temperature
        T = 293
        # define inf dilution
        dist = 1e-3
        z = np.array([dist, 1-dist])
        # get activity coefficient
        a = np.exp(COSMO.get_lngamma(T, z))
        gamma.append(a)
        # print(i)
        #names.append(prof.name)
        numbers.append(names)
        # print(COSMO.get_lngamma(T, z))

    df = pd.DataFrame(list(gamma), columns=["gamma1", "gamma2"])

    dg = pd.DataFrame(list(numbers), columns=['COMP', 'Solvent'])
    dh = pd.concat([df, dg], axis=1, join='inner')
    return dh

def CalcBinary_single(names , db, T):
    # careful with COSMO1 and COSMO2
    # call me with table['comb'] =  list([hexene,name1] for name1 in table['CAS#'])
    # list( CalcBinary_single(names,db,T) for names in table['comb'])
    import pandas as pd
    import cCOSMO
    import numpy as np
    # gamma = []
    # numbers = []
    # names = ['BENZENE','ETHANOL']
    # print(names)
    #loaded from harddrive and into virtual db
    # for iden in names:
    #     # print(iden)
    #     n = db.normalize_identifier(iden)
    #     db.add_profile(n)
    COSMO = cCOSMO.COSMO3(names, db)
    # here I should implement solvent dependent temperature
    T = 293
    # define inf dilution
    dist = 1e-3
    z = np.array([dist, 1-dist])
    # get activity coefficient
    a = np.exp(COSMO.get_lngamma(T, z))
    # dh2 =   pd.DataFrame([[a[0]],[a[1]],[names[0]],[names[1]]]).T.rename(columns={0:"gamma1", 1:"gamma2",2:'COMP', 3:'Solvent'})
    return a

