import headers as h
import scipy.optimize
import matplotlib.pyplot as plt

# needs a class for parameter/ dataset structure
# which params are given/ needed?
# are there different starting values available in literature?
# how many differentiations are needed? how do i implement these? --> chemical potential can be derived to gain additional information about the systems behaviour
# is there any need for eos, and or excess quantity models?


class fit_functions_binary:
    # calculating activity coefficient for margules
    @staticmethod
    def y_margules(_t, _x, _a):
        return 0

    # calculating activity coefficient for porter with 3 coefficients in gamma explicit form with no respect to component.
    # component information is given in the func arguments during call of func
    @staticmethod
    def y_porter(_T, _x, _a1, _a2, _a3):
        _func = h.np.exp((_a1+_a2/_T+_a3/_T**2)*_x**2)
        return _func

    # calculating x_i in phase beta from x_i in phase alpha via solution of equation system, both equations set to zero and with normation to 1
    # not suitable for koningsveld due to being in generalized form and not evolved out of porter in explicit form
    # solver works towards zero for both
    @staticmethod
    def _general_phase_partition_balance(_x_1_alpha, _x_1_beta, _T, _a1, _a2, _a3):
        _x_2_alpha = 1 - _x_1_alpha
        _x_2_beta = 1 - _x_1_beta

        _gamma_1_alpha = fit_functions_binary.y_porter(_T, _x_2_alpha, _a1, _a2, _a3)
        _gamma_1_beta = fit_functions_binary.y_porter(_T, _x_2_beta, _a1, _a2, _a3)
        _gamma_2_alpha = fit_functions_binary.y_porter(_T, _x_1_alpha, _a1, _a2, _a3)
        _gamma_2_beta = fit_functions_binary.y_porter(_T, _x_1_beta, _a1, _a2, _a3)

        _eq1 = ((_x_1_beta*_gamma_1_beta)/(_x_1_alpha*_gamma_1_beta))-1
        _eq2 = ((_x_2_beta*_gamma_2_beta)/(_x_2_alpha*_gamma_2_alpha))-1

        _func = _eq1 - _eq2
        return 0
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

    #criteria for optimization of system, quality criterium
    @staticmethod
    def min_fqs():
        return 0
