
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

        _func = _eq1 - _eq2
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
                              args=(_x_1_alpha, _x_1_beta, _texp, _steps), method='Nelder-Mead',bounds=((None,None),(None,None),(0,0)))
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
        print('ARD_normiert [%] =', ard_neu_norm)
        print('T-Calc [K] =', _tcalc)
        print('T-Exp [K] =', _texp)
        print('T-Diff [K] =', t_diff)
        _params = (_res.x[0], _res.x[1], _res.x[2], ard_neu_norm, _tcalc, _gamma,)
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
        a0 = 10
        a1 = 100
        a2 = 0

        a = h.np.array([a0, a1, a2])
        test_steps = 100
        b = h.np.zeros(test_steps)

        t_exp_01593 = h.np.array([273.15, 282.65, 292.75, 302.65, 312.45, 322.35, 332.35])
        x_chloroform_in_alpha_01593 = h.np.array([0.00155311, 0.00141498, 0.00126175, 0.00120053, 0.00112407, 0.00116994
                                                     , 0.00120053])
        x_chloroform_in_beta_01593 = h.np.array([0.997587, 0.996519, 0.995637, 0.994455, 0.992705, 0.991104, 0.989026])

        _x_1_alpha = x_chloroform_in_alpha_01593
        _x_1_beta = x_chloroform_in_beta_01593
        _texp = t_exp_01593
        _steps = len(_texp)
        for x in range(test_steps):
            b = ([a0+x, a1 + x, a2])
            print('Start-Parametersatz')
            print('A1-Start = ' + str(b[0]))
            print('A2-Start = ' + str(b[1]))
            print('A3-Start = ' + str(b[2]))
            h.pm1.FitFunctionsBinary._porter_parameter_fit(a, _x_1_alpha, _x_1_beta, _texp, _steps)

        return 0