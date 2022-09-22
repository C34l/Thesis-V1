import headers as h
import scipy.optimize
import matplotlib.pyplot as plt
from shapely.geometry import LineString

_r = 8.3141
_t0S = 404.75
_t0RS = 393.35
_t0R = 404.65
_h0S = 24500.0
_h0S_fit = 28100.0
_h0RS = 25600.0
_h0RS_fit = 31200.0
_h0R = 25736.0
#_h0R_fit = 29936.0
_h0R_fit = 28100.0
_Alpha = 0.4
XEXP = h.np.array([0.0, 0.1994, 0.3192, 0.4, 0.5, 0.5514, 0.605, 0.631, 0.6807, 0.6902, 0.7508, 0.7997, 0.8504, 0.9002, 0.9406, 1])
TEXP = h.np.array([404.75, 393.65, 387.55, 391.35, 393.35, 392.85, 391.35, 390.95, 388.35, 387.95, 391.75, 393.65, 397.15, 399.15, 400.65, 404.65])

class Diagrams:
    @staticmethod
    def y_nrtl(_g, _xi, _texp, _alpha):
        _xj = float(1.0 - _xi)
        _tauij = (_g[0] / (_r * _texp))
        _tauji = (_g[1] / (_r * _texp))
        _Gij = (h.np.exp(-_alpha * _tauij))
        _Gji = (h.np.exp(-_alpha * _tauji))
        _func = ((_xj ** 2) * (_tauji * (_Gji / (_xi + _xj * _Gji)) ** 2 + _tauij * (_Gij / (_xi * _Gij + _xj) ** 2)))

        return _func

    @staticmethod
    def y_porter(_x, _t, _a1, _a2):
        _func = h.np.exp((1 - _x) ** 2 * (_a1 + (_a2 / _t)))
        return _func

    @staticmethod
    def t_sle_nrtl(_texp, _xi, _h0, _t0, _g, _alpha):
        _yi = Diagrams.y_nrtl(_g, _xi, _texp, _alpha)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx) - 1
        return _func

    @staticmethod
    def t_sle_porter(_texp, _xi, _h0, _t0, _g, _alpha):
        _yi = Diagrams.y_porter(_g, _xi, _texp, _alpha)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx) - 1
        return _func

    @staticmethod
    def PD_ideal():
        return 0

    @staticmethod
    def PD_real_literatur():
        return 0

    @staticmethod
    def PD_Porter():
        return 0

    @staticmethod
    def PD_NRTL():
        return 0

    @staticmethod
    def Bilanz_A_porter():
        return 0

    @staticmethod
    def Bilanz_A_nrtl():
        return 0

    @staticmethod
    def Bilanz_B_porter():
        return 0

    @staticmethod
    def Bilanz_B_nrtl():
        return 0

    @staticmethod
    def Bilanz_C_porter():
        return 0

    @staticmethod
    def Bilanz_C_nrtl():
        return 0