
import headers as h
import scipy.optimize
import matplotlib.pyplot as plt

# needs a class for parameter/ dataset structure
# which params are given/ needed?
# are there different starting values available in literature?
# how many differentiations are needed? how do i implement these?
# --> chemical potential can be derived to gain additional information about the systems behaviour
# is there any need for eos, and or excess quantity models?

# we need to figure out a way to import our experimental data and add it here

class FitFunctionsBinary:
    # calculating activity coefficient for margules
    @staticmethod
    def y_margules(_t, _x, _a):
        return 0

    # calculating activity coefficient for porter with 3 coefficients in
    # gamma explicit form with no respect to component.
    # component information is given in the func arguments during call of func
    @staticmethod
    def y_porter(_a, _T, _x):
        _func = h.np.exp((_a[0]+_a[1]/_T+_a[2]/_T**2)*_x**2)
        return _func

    #solving porter bilance for T to remodel the whole phase diagram
    @staticmethod
    def y_porter_for_T(_T, _a, _x):
        _func = h.np.exp((_a[0] + _a[1] / _T + _a[2] / _T ** 2) * _x ** 2)
        return _func

    # calculating x_i in phase beta from x_i in phase alpha via solution of equation system,
    # both equations set to zero and with normation to 1
    # not suitable for koningsveld due to being in generalized form and not evolved out of porter in explicit form
    # solver works towards zero for both
    # params of a for porter as array
    # we're currently fitting the T way of things, whereas in literature oftentimes x based opti is used(!)
    @staticmethod
    def _general_phase_partition_balance_porter(_T, _x_1_alpha, _x_1_beta, _a,):
        _x_2_alpha = 1 - _x_1_alpha
        _x_2_beta = 1 - _x_1_beta

        _gamma_1_alpha = FitFunctionsBinary.y_porter(_a, _T, _x_2_alpha)
        _gamma_1_beta = FitFunctionsBinary.y_porter(_a, _T, _x_2_beta)
        _gamma_2_alpha = FitFunctionsBinary.y_porter(_a, _T, _x_1_alpha)
        _gamma_2_beta = FitFunctionsBinary.y_porter(_a, _T, _x_1_beta)

        _eq1 = ((_x_1_beta*_gamma_1_beta)/(_x_1_alpha*_gamma_1_beta))-1
        _eq2 = ((_x_2_beta*_gamma_2_beta)/(_x_2_alpha*_gamma_2_alpha))-1

        _func = 1 - (_eq1/_eq2)
        return _func

    # changed partition balance to account for more precise remodeling based on x_1_alpha and T as leading guesses.
    # this should help with the general problem of imprecise guesses for x_1_beta,
    # as x_1_alpha is often measured and T is constant for both comp points
    @staticmethod
    def _general_phase_partition_balance_porter_x(_x_1_beta, _x_1_alpha, _T, _a, ):
        _x_2_alpha = 1 - _x_1_alpha
        _x_2_beta = 1 - _x_1_beta

        _gamma_1_alpha = h.spo.fsolve(FitFunctionsBinary.y_porter_for_T, _T, args=(_a, _x_2_alpha))
        _gamma_1_beta = h.spo.fsolve(FitFunctionsBinary.y_porter_for_T, _T, args=(_a, _x_2_beta))
        _gamma_2_alpha = h.spo.fsolve(FitFunctionsBinary.y_porter_for_T, _T, args=(_a, _x_1_alpha))
        _gamma_2_beta = h.spo.fsolve(FitFunctionsBinary.y_porter_for_T, _T, args=(_a, _x_1_beta))

        _eq1 = ((_x_1_beta * _gamma_1_beta) / (_x_1_alpha * _gamma_1_beta)) - 1
        _eq2 = ((_x_2_beta * _gamma_2_beta) / (_x_2_alpha * _gamma_2_alpha)) - 1

        _func = 1 - (_eq1 / _eq2)
        return _func

    # target temperature difference for optimization based on minimum least square comparision
    # note that we carefully need to align our function parameters to pass the correct values
    # through to our optimization
    @staticmethod
    def _least_square_error_sum(_a, _x_1_alpha, _x_1_beta, _texp, _steps):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(h.pm1.FitFunctionsBinary._general_phase_partition_balance_porter,
                                   _texp[x], args=(_x_1_alpha[x], _x_1_beta[x], _a), full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]

        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    # method to calculate the porter coefficients via optimization
    # / minimization going the route of Tdiff into balance into gamma
    # cleanup of old declarations still needed!
    @staticmethod
    def _porter_parameter_fit(_a, _x_1_alpha, _x_1_beta, _texp, _steps):
        _steps = len(_texp)

        _res = h.spo.minimize(h.pm1.FitFunctionsBinary._least_square_error_sum, _a,
                              args=(_x_1_alpha, _x_1_beta, _texp, _steps), method='Powell',bounds=((None,None),(None,None),(0,0)))
        print('A1 = ' + str(_res.x[0]))
        print('A2 = ' + str(_res.x[1]))
        print('A3 = ' + str(_res.x[2]))


        t_diff_norm = h.np.zeros(_steps)
        t_diff = h.np.zeros(_steps)
        _tcalc = h.np.zeros(_steps)
        _gamma = h.np.zeros(_steps)

        _a_fit = (_res.x[0], _res.x[1], _res.x[2])

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(h.pm1.FitFunctionsBinary._general_phase_partition_balance_porter, _texp[x],
                                   args=(_x_1_alpha[x], _x_1_beta[x], _a_fit), full_output=True)
            _tcalc[x] = _tcalc1[0]
            t_diff[x] = abs((_texp[x] - _tcalc[x]))
            _gamma[x] = (h.pm1.FitFunctionsBinary.y_porter(_a_fit, _tcalc[x], _x_1_alpha[x]))

        t_diff_norm = h.np.abs(h.np.divide(t_diff, _texp))
        ard_neu_norm = (100 / _steps) * sum(t_diff_norm)
        if ard_neu_norm <= 10 and h.np.sum(t_diff) > 0:
            print('Calculation successful and optimization precise')
            print('ARD_normiert [%] =', ard_neu_norm)
            print('T-Calc [K] =', _tcalc)
            print('T-Exp [K] =', _texp)
            print('T-Diff [K] =', t_diff)
        elif ard_neu_norm <= 10 and h.np.sum(t_diff) == 0:
            print('Calculation failed')
        else:
            print('Optimization too imprecise for evaluation')
        #_params = (_res.x[0], _res.x[1], _res.x[2], ard_neu_norm, _tcalc, _gamma,)
        _params = ard_neu_norm
        return _params

    # calculating activity coefficient for koningsveld? or whole different model?
    @staticmethod
    def y_koningsveld(_t, _x, _a1, _a2, _a3, _p):
        return 0

    # calculating activity coefficient for nrtl
    @staticmethod
    def y_nrtl():
        return 0

    # basic function for every possible gamma, maybe need for optimized forms depending on model
    @staticmethod
    def eq_for_x():
        return 0

    # equation to solve EQs for T, again model dependent
    @staticmethod
    def eq_for_T():
        return 0

    # which are our two phases?
    #we're calculating the concentration of chcl3 in the two phases. if phase 1 is close to xchcl3 = 0
    # it's the waterphase, if phase 2 corresponds to xchcl3 = 1 phase 2, which means we're basically looking at
    # chloroform in the water rich phase and chloroform in the chloroform rich phase over the
    # molar concentration of chloroform in the entire system which means: SUMx_chlor_in_1+x_chlor_in_2 =/= 1
    #lets rename the water rich phase alpha and the chloroformrich phase beta
    #
    @staticmethod
    def method_caller():
        a0 = -100
        a1 = 100
        a2 = 0

        a = h.np.array([a0, a1, a2])
        test_steps = 400
        b = h.np.zeros(test_steps)

        #good dataset
        t_exp_01593 = h.np.array([273.15, 282.65, 292.75, 302.65, 312.45, 322.35, 332.35])
        x_chloroform_in_alpha_01593 = h.np.array([0.00155311, 0.00141498, 0.00126175, 0.00120053, 0.00112407, 0.00116994
                                                     , 0.00120053])
        x_chloroform_in_beta_01593 = h.np.array([0.997587, 0.996519, 0.995637, 0.994455, 0.992705, 0.991104, 0.989026])

        #highly imprecise dataset
        t_exp_019198 = h.np.array([283.15, 293.15, 298.15, 303.15, 313.15, 323.15])
        x_chloroform_in_alpha_019198 = h.np.array([0.0010, 0.001, 0.001, 0.002, 0.003, 0.004])
        x_chloroform_in_beta_019198 = h.np.array([0.992, 0.991, 0.989, 0.987, 0.984, 0.982])

        #good dataset but only halve the data needed, useful for referencing tho
        t_exp_015816 = h.np.array([273.15, 276.35, 290.55, 302.55, 314.75, 328.05])
        x_chloroform_in_alpha_015816 = h.np.array([0.00148727, 0.0013413, 0.00107333, 0.00115314, 0.00107333, 0.00116819])

        _x_1_alpha = x_chloroform_in_alpha_01593
        _x_1_beta = x_chloroform_in_beta_01593
        _texp = t_exp_01593
        _steps = len(_texp)
        for x in range(test_steps):
            b = ([a0, a1 + x*5, a2])
            print('Start-Parametersatz')
            print('A1-Start = ' + str(b[0]))
            print('A2-Start = ' + str(b[1]))
            print('A3-Start = ' + str(b[2]))
            h.pm1.FitFunctionsBinary._porter_parameter_fit(b, _x_1_alpha, _x_1_beta, _texp, _steps)


        return 0

    # plotting the wanted data sets and outputting them to excel
    @staticmethod
    def plotter():
        a0 = 6.64
        a1 = 1618.233
        a2 = 0

        a = h.np.array([a0, a1, a2])

        # good dataset
        t_exp_01593 = h.np.array([273.15, 282.65, 292.75, 302.65, 312.45, 322.35, 332.35])
        x_chloroform_in_alpha_01593 = h.np.array([0.00155311, 0.00141498, 0.00126175, 0.00120053, 0.00112407, 0.00116994
                                                     , 0.00120053])
        x_chloroform_in_beta_01593 = h.np.array([0.997587, 0.996519, 0.995637, 0.994455, 0.992705, 0.991104, 0.989026])

        # highly imprecise dataset
        t_exp_019198 = h.np.array([283.15, 293.15, 298.15, 303.15, 313.15, 323.15])
        x_chloroform_in_alpha_019198 = h.np.array([0.0010, 0.001, 0.001, 0.002, 0.003, 0.004])
        x_chloroform_in_beta_019198 = h.np.array([0.992, 0.991, 0.989, 0.987, 0.984, 0.982])

        # good dataset but only halve the data needed, useful for referencing tho
        t_exp_015816 = h.np.array([273.15, 276.35, 290.55, 302.55, 314.75, 328.05])
        x_chloroform_in_alpha_015816 = h.np.array(
            [0.00148727, 0.0013413, 0.00107333, 0.00115314, 0.00107333, 0.00116819])

        calc_steps = 20
        x_to_plot1 = h.np.zeros(calc_steps)
        x_to_plot2 = h.np.zeros(calc_steps)
        x_to_print1 = h.np.zeros(calc_steps)
        x_to_print2 = h.np.zeros(calc_steps)
        t_start1 = h.np.zeros(int(calc_steps / 2))
        t_start2 = h.np.zeros(int(calc_steps / 2))

        for x in range (int(calc_steps / 2)):
            t_start1[x] = 333.15 - 9 * x
            t_start2[x] = 333.15 - 4.5 * x

        t_start = h.np.append(t_start1, t_start2, axis=0)

        for x in range(int(calc_steps / 2)):
            x_to_plot1[x] = 0.001 + x * 0.0001
            x_to_plot2[x] = 0.989 + x * 0.0001

            x_to_print1[x] = h.spo.fsolve(h.pm1.FitFunctionsBinary._general_phase_partition_balance_porter_x, x_to_plot2[x],
                                        args=(x_to_plot1[x], t_start1[x], a))


        n1 = h.spo.fsolve(h.pm1.FitFunctionsBinary._general_phase_partition_balance_porter_x, 0.99,
                         args=(0.001109, 320, a))
        dev1 = ((0.99/n1)-1)*100
        n2 = h.spo.fsolve(h.pm1.FitFunctionsBinary._general_phase_partition_balance_porter_x, 0.99,
                          args=(0.00107, 325, a))
        dev2 = ((0.99 / n2) - 1) * 100
        n3 = h.spo.fsolve(h.pm1.FitFunctionsBinary._general_phase_partition_balance_porter_x, 0.99,
                          args=(0.00103, 330, a))
        dev3 = ((0.99 / n3) - 1) * 100
        print('x1', n1, 'dev1', dev1)
        print('x2', n2, 'dev2', dev2)
        print('x3', n3, 'dev3', dev3)

        #data = {'x2_start': x_to_plot2, 'x': x_to_print1, 't': t_start, 'x1_start': x_to_plot1, 't2': t_start}
        #dataframe = h.pd.DataFrame(data=data)
        #dataframe.to_excel(r'C:\Users\Ulf\Desktop\Promo\Daten\chloroform_water_x_based_5.xlsx')
        #print('export finished')

        return 0