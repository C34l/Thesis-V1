import headers as h
import scipy.optimize
import matplotlib.pyplot as plt
from shapely.geometry import LineString

_r = 8.3141
_t0S = 404.75
_t0RS = 393.35
_t0R = 404.65

_h0S = 24500.0
_h0RS = 25600.0
_h0R = 25736.0

_Alpha = 0.4
XEXP = h.np.array([0.0, 0.1994, 0.3192, 0.4, 0.5, 0.5514, 0.605, 0.631, 0.6807, 0.6902, 0.7508, 0.7997, 0.8504, 0.9002, 0.9406, 1])
TEXP = h.np.array([404.75, 393.65, 387.55, 391.35, 393.35, 392.85, 391.35, 390.95, 388.35, 387.95, 391.75, 393.65, 397.15, 399.15, 400.65, 404.65])

_a1a = -5.597
_a2a = 1952.672
_aSa = h.np.array([_a1a, _a2a, ])

_a1b = -32.196
_a2b = 12794.746
_aRa = h.np.array([_a1b, _a2b, ])

_gabS = -1403.640
_gbaS = 2972.166
_gSa = h.np.array([_gabS, _gbaS, ])

_gabRS = -3143.16
_gbaRS = 9443.10
_gRSa = h.np.array([_gabRS, _gbaRS, ])

_gabR = -23.289
_gbaR = 7526.385
_gRa = h.np.array([_gabR, _gbaR, ])


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

    #ana
    @staticmethod
    def PD_ideal(x, t):
        _func = h.np.log(1/((4*x)-4*x**2))-((_h0RS/_r)*((1/t)-(1/_t0RS)))
        return _func

    #num
    @staticmethod
    def PD_real_literatur(c, x, t):
        _func = -(t/_h0RS)*((x/(1-x))-1)*(((_r*t)/x)-2*c*(1-x))
        return _func

    #num
    @staticmethod
    def PD_Porter_pia(a, x, t):
        _func = ((_r*t**2)/_h0RS)*((1/(1-x))-(1/x)+2*a[0]*(-2*x+1)+2*a[1]*(((-2*x+1))/t))
        return _func

    #num
    @staticmethod
    def PD_NRTL_pia(g, x, t):
        _func =
        return 0

    @staticmethod
    def PD_Porter_fit():
        return 0

    @staticmethod
    def PD_NRTL_fit():
        return 0

    @staticmethod
    def Bilanz_A_porter_pia():
        return 0

    @staticmethod
    def Bilanz_A_nrtl_pia():
        return 0

    @staticmethod
    def Bilanz_A_porter_fit():
        return 0

    @staticmethod
    def Bilanz_A_nrtl_fit():
        return 0

    @staticmethod
    def Bilanz_B_porter_pia():
        return 0

    @staticmethod
    def Bilanz_B_nrtl_pia():
        return 0

    @staticmethod
    def Bilanz_B_porter_fit():
        return 0

    @staticmethod
    def Bilanz_B_nrtl_fit():
        return 0

    @staticmethod
    def Bilanz_C_porter_pia():
        return 0

    @staticmethod
    def Bilanz_C_nrtl_pia():
        return 0

    @staticmethod
    def Bilanz_C_porter_fit():
        return 0

    @staticmethod
    def Bilanz_C_nrtl_fit():
        return 0