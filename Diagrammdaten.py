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
XEXP_korr = h.np.array([0.0, 0.1994, 0.3192, 0.4, 0.5, 0.486, 0.395, 0.369, 0.3193, 0.3098, 0.2492, 0.2003, 0.1496, 0.0998, 0.0594, 0])
XEXP_korr_links = h.np.array([0.0, 0.1994, 0.3192, 0.4, 0.5, ])
XEXP_korr_rechts = h.np.array([0.5, 0.486, 0.395, 0.369, 0.3193, 0.3098, 0.2492, 0.2003, 0.1496, 0.0998, 0.0594, 0])
TEXP = h.np.array([404.75, 393.65, 387.55, 391.35, 393.35, 392.85, 391.35, 390.95, 388.35, 387.95, 391.75, 393.65, 397.15, 399.15, 400.65, 404.65])
TEXP_links = h.np.array([404.75, 393.65, 387.55, 391.35, 393.35, ])
TEXP_rechts = h.np.array([393.35, 392.85, 391.35, 390.95, 388.35, 387.95, 391.75, 393.65, 397.15, 399.15, 400.65, 404.65])

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
    def PD_ideal(t, x):
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

    # num
    @staticmethod
    def PD_Porter_pia_num(x, t, a):
        _tinitial = 393.35

        _initial = [_tinitial]
        _func = h.spi.solve_ivp(Diagrams.PD_Porter_pia, [0.5, 1.0], _initial, method='RK45',
                                    args=(_Alpha, a[0], a[1], _h0RS), t_eval=x, dense_output=True)
        return _func

    #num
    @staticmethod
    def PD_NRTL_pia(g, x, t):
        a = _Alpha
        u = g[0]
        U = h.np.exp((a*u)/(_r*t))
        v = g[1]
        V = h.np.exp((a*v)/(_r*t))

        _func = ((x/(1-x))-1)*(t/_h0RS)*(((_r*t)/x)+(2*(1-x)**2)*(((u*(U-1))/(x*(U-1)+1)**3)-((v*V*(V-1))/((x-1)*V-x)**3))-(2*(1-x))*((u/(x*(U-1)+1)**2)+((v*V)/(x-(x-1)*V)**2)))
        return _func

    #num
    @staticmethod
    def PD_NRTL_pia_num(x, t, g):
        _tinitial = 393.35

        _initial = [_tinitial]
        _func = h.spi.solve_ivp(Diagrams.PD_NRTL_pia, [0.5, 1.0], _initial, method='RK45',
                                args=(_Alpha, g[0], g[1], _h0RS), t_eval=x, dense_output=True)
        return _func

    @staticmethod
    def PD_Porter_fit_minfqs(_a, _xi, _texp, _h0, _t0, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.PD_Porter_pia_num, _texp[x], args=(_xi[x], _h0, _t0, _a,), full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    @staticmethod
    def PD_NRTL_fit_minfqs(_g, _xi, _texp, _h0, _t0, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.PD_NRTL_pia_num, _texp[x], args=(_xi[x], _h0, _t0, _g,), full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #ana
    @staticmethod
    def Bilanz_A_porter_pia(t, x):
        ya = Diagrams.y_porter(x, t, _aSa[0], _aSa[1])
        yb = Diagrams.y_porter(x, t, _aRa[0], _aRa[1])
        _func = (((1/2)*(((-_h0S/(_r*t))*(1-(t/_t0S)))+(((-_h0R/(_r*t))*(1-(t/_t0R))))))/((h.np.log(0.25)-h.np.log(x*ya*(1-x)*yb))))-1
        return _func

    @staticmethod
    def Bilanz_A_nrtl_pia(t, x):
        ya = Diagrams.y_nrtl(_gSa, x, t, _Alpha)
        yb = Diagrams.y_nrtl(_gRa, x, t, _Alpha)
        _func = (((1 / 2) * (((-_h0S / (_r * t)) * (1 - (t / _t0S))) + (((-_h0R / (_r * t)) * (1 - (t / _t0R)))))) / (
        (h.np.log(0.25) - h.np.log(x * ya * (1 - x) * yb)))) - 1
        return _func

    @staticmethod
    def Bilanz_A_porter_fit_minfqs(_a, _xi, _texp, _h0, _t0, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_A_porter_pia, _texp[x], args=(_xi[x], _h0, _t0, _a,), full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    @staticmethod
    def Bilanz_A_nrtl_fit_minfqs(_g, _xi, _texp, _h0, _t0, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_A_NRTL_pia, _texp[x], args=(_xi[x], _h0, _t0, _g,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

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