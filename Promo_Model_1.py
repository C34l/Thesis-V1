import headers as h
import scipy.optimize
import matplotlib.pyplot as plt

# needs a class for parameter/ dataset structure
# which params are given/ needed?
# are there different starting values available in literature?
# how many differentiations are needed? how do i implement these? --> chemical potential can be derived to gain additional information about the systems behaviour
class fit_functions_binary:
    #calculating activity coefficient for margules
    @staticmethod
    def y_margules(_t, _x, _a):
        return 0

    # calculating activity coefficient for porter with 3 coefficients in gamma explicit form
    @staticmethod
    def y_porter(_t, _x, _a1, _a2, _a3):
        _func = h.np.exp((_a1+_a2/_t+_a3/_t**2)*(1-_x)**2)
        return _func

    #calculating x_i in phase beta from x_i in phase alpha via solution of equationsystem
    @staticmethod
    def _phase_partition_balance():
        _eq1 = 0
        _eq2 = 0

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
