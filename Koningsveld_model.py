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

    con4 = (h.np.log((1-x_1_alpha)/(1-x_1_beta)))
    con5 = x_1_beta**2 * (1-p) / (1-x_1_beta*p)**2
    con6 = x_1_alpha**2 * (1-p) / (1-x_1_alpha*p)**2
    eq2 = con4/(con5-con6)

    a = 0.9499 + 12.587/_texp

    #_func = 1 - (eq2/eq1)
    _func = eq1 - a

    return _func

def x_crit(x,p):
    a = (2*x-1)/(x*(x-2))
    func = 1 - p / a
    return func

def phase_balance_for_x( x_1_beta, x_1_alpha, _texp):
    _x_2_alpha = 1 - x_1_alpha
    _x_2_beta = 1 - x_1_beta
    p = 0.935

    con1 = (h.np.log(x_1_alpha/x_1_beta))
    con2 = ((1-x_1_beta)**2)/((1-x_1_beta*p)**2)
    con3 = ((1-x_1_alpha)**2)/((1-x_1_alpha*p)**2)
    eq1 = con1/(con2-con3)

    con4 = (h.np.log(1-x_1_alpha/1-x_1_beta))
    con5 = x_1_beta**2 * (1-p) / (1-x_1_beta*p)**2
    con6 = x_1_alpha**2 * (1-p) / (1-x_1_alpha*p)**2
    eq2 = con4/(con5-con6)

    a = 0.9499 + 12.587/_texp

    #_func = 1 - (eq2/eq1)
    _func = eq1 - a

    return _func

guess_p = -6
guess_x = 0.2
t_exp_01593 = h.np.array([273.15, 282.65, 292.75, 302.65, 312.45, 322.35, 332.35])
x_1_alpha = h.np.array([0.00155311, 0.00141498, 0.00126175, 0.00120053, 0.00112407, 0.00116994, 0.00120053])
x_1_beta = h.np.array([0.997587, 0.996519, 0.995637, 0.994455, 0.992705, 0.991104, 0.989026])

x_1_alpha_test = h.np.array([.0001,0.00011,0.00012,0.00013,0.00014,0.00015,0.00016,0.00017,0.0002,0.001,0.002,0.003,0.004])
x_1_beta_test = h.np.array([.997,0.995,0.993,0.991,0.988,0.986,0.984,0.982])
t_test = h.np.array([273.15,280,290,300,310,320,330,340,360,380,400,420,450])
steps = len(x_1_alpha_test)
x_output_beta = h.np.zeros(steps)
x_output_alpha = h.np.zeros(steps)
x_diff_beta = h.np.zeros(steps)
x_diff_alpha = h.np.zeros(steps)
p=h.np.array([1])

p_output = h.spo.least_squares(phase_balance, guess_p, args=(x_1_alpha,x_1_beta,t_exp_01593),bounds=((-20,1)))
x_c = h.spo.fsolve(x_crit, guess_x, args=(p_output.x))
print('value for p:', p_output.x)
print('value for xc:', x_c)

# for x in range(steps):
#     x_output_beta[x] = h.spo.fsolve(phase_balance_for_x, 1-2*x_1_alpha_test[x], args=(x_1_alpha_test[x], t_test[x]))
#     print('point_beta:',x_output_beta[x])
#
# for x in range(steps):
#     x_diff_beta[x] = x_output_beta[x] - x_1_beta[x]
#     xpercent = ((x_output_beta[x] / x_1_beta[x]) - 1) * 100
#     print('error_pointwise_percentile_beta', xpercent)
#
# fqs_norm_beta = (h.np.abs(h.np.divide(x_diff_beta, x_1_beta))) ** 2
# fqs_summe_beta = h.np.sum(fqs_norm_beta)
# print('error_sum_beta:', fqs_summe_beta)

# for x in range(steps):
#     x_output_alpha[x] = h.spo.fsolve(phase_balance_for_x, x_1_alpha_test[x], args=(x_1_beta[x], t_test[x]))
#     print('point_alpha:', x_output_alpha[x])
#
# for x in range(steps):
#     x_diff_alpha[x] = x_output_alpha[x] - x_1_alpha[x]
#     xpercent = ((x_output_alpha[x] / x_1_alpha[x]) - 1) * 100
#     print('error_pointwise_percentile_alpha', xpercent)
#
# fqs_norm_alpha = (h.np.abs(h.np.divide(x_diff_alpha, x_1_alpha))) ** 2
# fqs_summe_alpha = h.np.sum(fqs_norm_alpha)
# print('error_sum_alpha:', fqs_summe_alpha)