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

    a = 0.99 + 2.022 / _texp

    # con1 = (h.np.log(x_1_alpha/x_1_beta))
    # con2 = ((1-x_1_beta)**2)/((1-x_1_beta*p)**2)
    # con3 = ((1-x_1_alpha)**2)/((1-x_1_alpha*p)**2)
    # eq1 = con1/(con2-con3)
    #
    # con4 = (h.np.log((1-x_1_alpha)/(1-x_1_beta)))
    # con5 = x_1_beta**2 * (1-p) / (1-x_1_beta*p)**2
    # con6 = x_1_alpha**2 * (1-p) / (1-x_1_alpha*p)**2
    # eq2 = con4/(con5-con6)

    con1 = (h.np.log(x_1_alpha / x_1_beta))
    con2 = ((1 - x_1_beta) ** 2) / ((1 - x_1_beta * p) ** 2)
    con3 = ((1 - x_1_alpha) ** 2) / ((1 - x_1_alpha * p) ** 2)
    eq1 = 1 - (con1 / ((con2 - con3) * a))

    con4 = (h.np.log(_x_2_alpha / _x_2_beta))
    con5 = x_1_beta ** 2 * (1 - p) / (1 - x_1_beta * p) ** 2
    con6 = x_1_alpha ** 2 * (1 - p) / (1 - x_1_alpha * p) ** 2
    eq2 = 1 - con4 / ((con5 - con6) * a)

    _func = 1 - (eq1/eq2)
    #_func = eq1 - a

    return _func

def x_crit(x,p):
    a = - (2*(x-1))/(x*(x-2))
    func = 1 - p / a
    return func

def phase_balance_for_x( x_1_beta, x_1_alpha, _texp):
    _x_2_alpha = 1 - x_1_alpha
    _x_2_beta = 1 - x_1_beta
    p = -0.501
    a = 0.90233 + 25.1624/ _texp

    con1 = (h.np.log(x_1_alpha/x_1_beta))
    con2 = ((1-x_1_beta)**2)/((1-x_1_beta*p)**2)
    con3 = ((1-x_1_alpha)**2)/((1-x_1_alpha*p)**2)
    eq1 = 1 - (con1/((con2-con3)*a))

    con4 = (h.np.log(_x_2_alpha/_x_2_beta))
    con5 = x_1_beta**2 * (1-p) / (1-x_1_beta*p)**2
    con6 = x_1_alpha**2 * (1-p) / (1-x_1_alpha*p)**2
    eq2 = 1 - con4/((con5-con6)*a)

    _func = 1 - (eq1/eq2)
    #_func = eq1 - a

    return _func

#datablock
guess_p = -1
guess_x = 0.4

t_exp_01593 = h.np.array([273.15, 282.65, 292.75, 302.65, 312.45, 322.35, 332.35])
t_exp_019198 = h.np.array([283.15, 293.15, 298.15, 303.15, 313.15, 323.15])
x_1_alpha_019198 = h.np.array([0.0008,0.0009,0.0011,0.0013,0.0016,0.0018])
x_1_beta_019198 = h.np.array([0.992,0.991,0.989,0.987,0.984,0.982])
x_1_alpha = h.np.array([0.00155311, 0.00141498, 0.00126175, 0.00120053, 0.00112407, 0.00116994, 0.00120053])
#01593 values
x_1_beta = h.np.array([0.997587, 0.996519, 0.995637, 0.994455, 0.992705, 0.991104, 0.989026])
#005740 values
#x_1_beta = h.np.array([])

# x_1_alpha_test = h.np.array([0.0017326, 0.001451498, 0.00116175, 0.00130053, 0.001002407, 0.001545845, 0.0011878455])
# x_1_beta_test = h.np.array([.9971212,0.9912235,0.99312312,0.991122,0.9883123,0.981231236,0.981231234])
x_1_alpha_test = h.np.array([0.0017326, 0.001451498, 0.00116175, 0.00130053, 0.001002407, 0.001545845])
x_1_beta_test = h.np.array([.9971212,0.9912235,0.99312312,0.991122,0.9883123,0.981231236])
t_test = h.np.array([273.15,280,290,300,310,320,330,340,360,380,400,420,450])
steps = len(x_1_alpha)
x_output_beta = h.np.zeros(steps)
x_output_alpha = h.np.zeros(steps)
x_diff_beta = h.np.zeros(steps)
x_diff_alpha = h.np.zeros(steps)
xpercent = h.np.zeros(steps)
p = h.np.array([1])

# #calculation of p and x crit
p_output = h.spo.least_squares(phase_balance, guess_p, args=(x_1_alpha_019198, x_1_beta_019198, t_exp_019198),bounds=((-20,1)),  method='trf')
x_c = h.spo.fsolve(x_crit, guess_x, args=(p_output.x))
print('value for p:', p_output.x)
print('value for xc:', x_c)

#gives values for x1 beta
# for x in range(steps):
#     x_output_beta[x] = h.spo.root(phase_balance_for_x, x_1_beta[x], args=(x_1_alpha[x], t_exp_01593[x]))
#     print('point_beta:',x_output_beta[x],'temp:',t_exp_01593[x])
x_output_alpha= h.spo.root(phase_balance_for_x, x_1_alpha_test, args=(x_1_beta_019198, t_exp_019198), method='lm')
print('point_alpha:',x_output_alpha.x,'temp:',t_exp_019198)

args_alpha = x_output_alpha.x

print('array alpha:',args_alpha)

x_output_beta = h.spo.root(phase_balance_for_x, x_1_beta_test, args=(x_1_alpha_019198, t_exp_019198), method='lm')
print('point_beta:',x_output_beta.x,'temp:',t_exp_019198)

# for x in range(steps):
#     x_diff_beta[x] = x_output_beta[x] - x_1_beta[x]
#     xpercent[x] = ((x_output_beta[x] / x_1_beta[x]) - 1) * 100
#     print('error_pointwise_percentile_beta', xpercent[x])

# fqs_test = h.np.sum(xpercent)/steps
# print('median error:', fqs_test)
#
# fqs_norm_beta = (h.np.abs(h.np.divide(x_diff_beta, x_1_beta))) ** 2
# fqs_summe_beta = h.np.sum(fqs_norm_beta)
# print('error_sum_beta:', fqs_summe_beta)

#gives values for x1 alpha

# for x in range(steps):
#     x_output_alpha[x] = h.spo.fsolve(phase_balance_for_x, x_1_alpha[x]+0.0001, args=(x_1_beta[x], t_exp_01593[x]))
#     print('point_alpha:', x_output_alpha[x])
#
# for x in range(steps):
#     x_diff_alpha[x] = x_output_alpha[x] - x_1_alpha[x]
#     xpercent = ((x_output_alpha[x] / x_1_alpha[x]) - 1) * 100
#     print('error_pointwise_percentile_alpha', xpercent)
#
# fqs_norm_alpha = (h.np.abs(h.np.divide(x_diff_alpha, x_1_alpha)))**2
# fqs_summe_alpha = h.np.sum(fqs_norm_alpha)
# print('error_sum_alpha:', fqs_summe_alpha)