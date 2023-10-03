import headers as h
import scipy.optimize
import matplotlib.pyplot as plt

def objective_function(_xAexp, _xBexp, _steps, _texp):
    _xAcalc = h.np.zeros(_steps)
    _xAdiff = h.np.zeros(_steps)
    _xBcalc = h.np.zeros(_steps)
    _xBdiff = h.np.zeros(_steps)

    for x in range(_steps):
        _xcalc1 = h.spo.fsolve
        _xcalc2 = h.spo.fsolve
        _xAcalc[x] = _xcalc1[0]
        _xBcalc[x] = _xcalc2[0]
        _xAdiff[x] = _xAexp[x] - _xAcalc[x]
        _xBdiff[x] = _xBexp[x] - _xBcalc[x]


    OF = (h.np.abs(h.np.divide(_xAdiff, _xAexp))) ** 2 + (h.np.abs(h.np.divide(_xBdiff, _xBexp))) ** 2
    return OF


def phase_balance(x_1_alpha, x_1_beta, p):
    _x_2_alpha = 1 - x_1_alpha
    _x_2_beta = 1 - x_1_beta

    con1 = (h.np.log(x_1_alpha/x_1_beta))
    con2 = ((1-x_1_beta)**2)/((1-x_1_beta*p)**2)
    con3 = ((1-x_1_alpha)**2)/((1-x_1_alpha*p)**2)
    eq1 = con1/(con2-con3)

    con4 =
    con5 =
    con6 =

    return 0