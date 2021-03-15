
from calc_KPI_gmehling import calc_KPI_gmehling



def test_calc_KPI_gmehling():
    #Four parts first capacaity
    # solvent loss
    # this is the validation with data from gmeling
    x_SR = 0.01
    alpha_dis = 1/0.46
    ws = 1
    wc = 1
    wsl = 3
    wsd = 0.75
    # distribution coefficient of components in extract =gamma_LM1/gamma_LM2
    K_E_i = np.array([0.4843,	0.3625])
    # distribution coefficient of components in raffinat - here product
    K_R_j = 0.0106
    
    #number of components
    num_components = len(K_E_i)+1
    #capacity = first part of formula
    capacity =  sum(K_E_i/K_R_j)/(num_components-1)
    
    # selectivity
    selectivity = sum(K_E_i)/len(K_E_i)
    
    # solvent loss of raffinat solvet
    solvent_loss = 1-x_SR
    
    #  distillation difficulty
    distillation = np.log10(alpha_dis)
    
    res = capacity**ws * selectivity**wc * solvent_loss**wsl * distillation**wsd
    vec = [ws, wc, wsl, wsd]
    compare = res - calc_KPI_gmehling(K_E_i, K_R_j, x_SR, alpha_dis, vec )
    
    assert compare==0 , 'test_calc_KPI_gmehling_failed'