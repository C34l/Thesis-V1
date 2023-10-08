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

    a = -.8 + 320 / _texp

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

    con4 = (h.np.log(1 - x_1_alpha / 1 - x_1_beta))
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
    p = -.3125
    a = -.8 + 320 / _texp

    con1 = (h.np.log(x_1_alpha/x_1_beta))
    con2 = ((1-x_1_beta)**2)/((1-x_1_beta*p)**2)
    con3 = ((1-x_1_alpha)**2)/((1-x_1_alpha*p)**2)
    eq1 = 1 - (con1/((con2-con3)*a))

    con4 = (h.np.log(1-x_1_alpha/1-x_1_beta))
    con5 = x_1_beta**2 * (1-p) / (1-x_1_beta*p)**2
    con6 = x_1_alpha**2 * (1-p) / (1-x_1_alpha*p)**2
    eq2 = 1 - con4/((con5-con6)*a)

    _func = 1 - (eq1/eq2)
    #_func = eq1 - a

    return _func

guess_p = -.3
guess_x = 0.2
t_exp_01593 = h.np.array([200, 282.65, 292.75, 302.65, 332.45, 350.35, 400])
x_1_alpha = h.np.array([0.1, 0.1513, 0.2273, 0.2551, 0.2997, 0.3488, 0.4])
x_1_beta = h.np.array([0.899, 0.832, 0.799, 0.734, 0.6247, .4999, 0.3998])

steps = len(x_1_alpha)
x_output_beta = h.np.zeros(steps)
x_output_alpha = h.np.zeros(steps)
x_diff_beta = h.np.zeros(steps)
x_diff_alpha = h.np.zeros(steps)
p = h.np.array([1])

# p_output = h.spo.least_squares(phase_balance, guess_p, args=(x_1_alpha, x_1_beta, t_exp_01593),bounds=((-20,1)),  method='trf')
# x_c = h.spo.fsolve(x_crit, guess_x, args=(p_output.x))
# print('value for p:', p_output.x)
# print('value for xc:', x_c)

for x in range(steps):
    x_output_beta[x] = h.spo.fsolve(phase_balance_for_x, x_1_beta[x]+0.007, args=(x_1_alpha[x], t_exp_01593[x]))
    print('point_beta:',x_output_beta[x])

for x in range(steps):
    x_diff_beta[x] = x_output_beta[x] - x_1_beta[x]
    xpercent = ((x_output_beta[x] / x_1_beta[x]) - 1) * 100
    print('error_pointwise_percentile_beta', xpercent)

fqs_norm_beta = (h.np.abs(h.np.divide(x_diff_beta, x_1_beta))) ** 2
fqs_summe_beta = h.np.sum(fqs_norm_beta)
print('error_sum_beta:', fqs_summe_beta)

for x in range(steps):
    x_output_alpha[x] = h.spo.fsolve(phase_balance_for_x, x_1_alpha[x]+0.007, args=(x_1_beta[x], t_exp_01593[x]))
    print('point_alpha:', x_output_alpha[x])

for x in range(steps):
    x_diff_alpha[x] = x_output_alpha[x] - x_1_alpha[x]
    xpercent = ((x_output_alpha[x] / x_1_alpha[x]) - 1) * 100
    print('error_pointwise_percentile_alpha', xpercent)

fqs_norm_alpha = (h.np.abs(h.np.divide(x_diff_alpha, x_1_alpha)))**2
fqs_summe_alpha = h.np.sum(fqs_norm_alpha)
print('error_sum_alpha:', fqs_summe_alpha)

x_Butanol_left = [0.01650, 0.01900, 0.01600, 0.01800, 0.01750, 0.01750, 0.01800, 0.01700, 0.01700, 0.01750, 0.01750, 0.01600, 0.01700, 0.01600, 0.01850, 0.01250, 0.01850, 0.01850, 0.02250, 0.02250, 0.02250, 0.02450, 0.02500, 0.02850, 0.03100, 0.03450, 0.03700, 0.04100, 0.04650, 0.05300, 0.05700, 0.07750, 0.08450, 0.09700, 0.10200]
T_butanol_left = [290.55649, 292.41147, 294.94098, 296.12142, 297.63912, 299.49410, 302.02361, 303.87858, 307.58853, 310.45531, 312.98482, 323.60877, 333.55818, 338.61720, 344.01349, 348.90388, 353.45700, 355.81788, 363.23777, 367.11636, 371.83811, 373.18718, 380.60708, 383.81113, 386.50927, 388.02698, 391.39966, 392.58010, 394.94098, 396.96459, 397.30185, 397.63912, 397.97639, 398.48229, 398.65093]
x = [0.4915, 0.4900, 0.4890, 0.4875, 0.4855, 0.4830, 0.4810, 0.4785, 0.4740, 0.4715, 0.4575, 0.4370, 0.4390, 0.4190, 0.4035, 0.3675, 0.3590, 0.3525, 0.3315, 0.3240, 0.2990, 0.2920, 0.2700, 0.2500, 0.2395, 0.2095, 0.1955, 0.1925, 0.1860, 0.1695, 0.1610, 0.1455, 0.1290, 0.1195, 0.1105, 0.1020, 0.09700]
y = [291.1, 293.4, 296.5, 297.8, 299.8, 302.0, 305.2, 308.4, 311.0, 313.7, 323.1, 331.9, 334.1, 343.5, 353.6, 360.5, 363.4, 365.8, 370.5, 372.3, 379.6, 379.9, 384.1, 387.9, 388.7, 393.4, 394.6, 395.8, 395.4, 397.3, 397.5, 397.5, 398.5, 398.7, 398.8, 398.7, 398.5]