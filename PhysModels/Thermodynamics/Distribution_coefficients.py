def CreateDist(df):
    import numpy as np
    # calculate distribution coefficients for inf dilution
    values = df['gamma1'].values
    # print(len(values))
    res = np.zeros((len(values), len(values)))
    i = 0
    j = 0
    length = len(values)
    for k in range(length):
        for m in range(length):
            res[k, m] = values[m]/values[k]

    return res

def CreateDist_acc(df):
    import numpy as np
    import itertools
    # calculate distribution coefficients for inf dilution
    # accelerated version
    # this is calculated by K = x1/x2 = gamma_(inf,2)/gamma_(inf,2) (see Gmehling Paper Solvent choice)
    # input is dataframe with column gamma1
    values = df['gamma1'].values
    res = np.reshape(np.array(list((b/a for a,b in itertools.product(values,repeat = 2)))),(len(values), len(values)))

    return res