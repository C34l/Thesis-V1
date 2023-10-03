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


def phase_balance(p, x_1_alpha, x_1_beta, _texp):
    _x_2_alpha = 1 - x_1_alpha
    _x_2_beta = 1 - x_1_beta

    con1 = (h.np.log(x_1_alpha/x_1_beta))
    con2 = ((1-x_1_beta)**2)/((1-x_1_beta*p)**2)
    con3 = ((1-x_1_alpha)**2)/((1-x_1_alpha*p)**2)
    eq1 = con1/(con2-con3)

    con4 = (h.np.log(1-x_1_alpha/1-x_1_beta))
    con5 = x_1_beta**2 * (1-p) / (1-x_1_beta*p)**2
    con6 = x_1_alpha**2 * (1-p) / (1-x_1_alpha*p)**2
    eq2 = con4/(con5-con6)

    a= 0.9499 + 12.587/_texp

    #_func = 1 - (eq2/eq1)
    _func = eq1 - a

    return _func

def x_crit(x,p):
    a = 2*x-1/(x*(x-2))
    func = 1 - a / p
    return func


guess_p = 0.9
guess_x = 0.2
t_exp_01593 = h.np.array([273.15, 282.65, 292.75, 302.65, 312.45, 322.35, 332.35])
x_1_alpha = h.np.array([0.00155311, 0.00141498, 0.00126175, 0.00120053, 0.00112407, 0.00116994, 0.00120053])
x_1_beta = h.np.array([0.997587, 0.996519, 0.995637, 0.994455, 0.992705, 0.991104, 0.989026])
steps = len(x_1_alpha)



p_output = h.spo.least_squares(phase_balance, guess_p, args=(x_1_alpha,x_1_beta,t_exp_01593),bounds=((-20,1)))
x_c = h.spo.fsolve(x_crit, guess_x, args=(p_output.x))
print('value for p:', p_output.x)
print('value for xc:', x_c)



