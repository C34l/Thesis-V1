import headers as h
import scipy.optimize as spo
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

XEXP_rechts = h.np.array([0.5, 0.5514, 0.605, 0.631, 0.6807, 0.6902, 0.7508, 0.7997, 0.8504, 0.9002, 0.9406, 1])
XEXP_korr = h.np.array([0.0, 0.1994, 0.3192, 0.4, 0.5, 0.486, 0.395, 0.369, 0.3193, 0.3098, 0.2492, 0.2003, 0.1496, 0.0998, 0.0594, 0])
XEXP_korr_links_Porter = h.np.array([0.000000000000001, 0.1994, 0.3192, 0.4, 0.5, ])
XEXP_korr_links= h.np.array([0.000000001, 0.1994, 0.3192, 0.4, 0.5, ])
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

_aGes = h.np.array([_a1a, _a2a, _a1b, _a2b])


_gabS = -1403.640
_gbaS = 2972.166
_gSa = h.np.array([_gabS, _gbaS, ])

_gabRS = -3143.16
_gbaRS = 9443.10
_gRSa = h.np.array([_gabRS, _gbaRS, ])

_gabR = -23.289
_gbaR = 7526.385
_gRa = h.np.array([_gabR, _gbaR, ])

_gGes = h.np.array([_gabS, _gbaS, _gabR, _gbaR])

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
    def y_porter(_x, _t, _a):
        _func = h.np.exp((1 - _x) ** 2 * (_a[0] + (_a[1] / _t)))
        return _func

    @staticmethod
    def t_sle_nrtl(_texp, _xi, _h0, _t0, _g, _alpha):
        _yi = Diagrams.y_nrtl(_g, _xi, _texp, _alpha)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx) - 1
        return _func

    @staticmethod
    def t_sle_porter(_texp, _xi, _h0, _t0, _g, _alpha):
        _yi = Diagrams.y_porter(_g, _xi, _texp,)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx) - 1
        return _func

    #ana, löst PD ideal
    @staticmethod
    def PD_ideal(t, x):
        _func = h.np.log(1/((4*x)-4*x**2))-((_h0RS/_r)*((1/t)-(1/_t0RS)))
        return _func

    #num
    @staticmethod
    def PD_real_literatur(c, x, t):
        _func = -(t/_h0RS)*((x/(1-x))-1)*(((_r*t)/x)-2*c*(1-x))
        return _func

    #num, berechne t aus PD-real
    @staticmethod
    def PD_real_literatur_num(t, x, a,):
        _tinitial = 393.35

        _initial = [_tinitial]
        _func = h.spi.solve_ivp(Diagrams.PD_real_literatur, [0.0, 1.0], _initial, method='RK45', args=(a,),
                                dense_output=True)
        a = h.np.ndarray.flatten(_func.y)
        return a

    # soll PD-real optimieren und c bestimmen
    @staticmethod
    def PD_real_literatur_minfqs(c, _xi, _texp, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
                _tcalc1 = h.spo.fsolve(Diagrams.PD_real_literatur_num, _texp[x], args=(_xi[x], c,), full_output=True)
                _tcalc[x] = _tcalc1[0]
                _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #num, Gleichung PD Porter Pia_
    @staticmethod
    def PD_Porter_pia(x, t, a1, a2):
        _func = ((_r*t**2)/_h0RS)*((1/(1-x))-(1/x)+2*a1*(-2*x+1)+2*a2*(((-2*x+1))/t))
        return _func

    # num, löst PD porter Pia
    @staticmethod
    def PD_Porter_pia_num(t, x, a):
        _tinitial = 393.35

        _initial = [_tinitial]
        _func = h.spi.solve_ivp(Diagrams.PD_Porter_pia, [0.5, 1.0], _initial, method='RK45', args=(a, ), t_eval=x, dense_output=True)
        return _func

    # soll PD-Porter optimieren und a bestimmen
    @staticmethod
    def PD_Porter_fit_minfqs(_a, _xi, _texp, _h0, _t0, _steps, ):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.PD_Porter_pia_num, _texp[x], args=(_xi, _a,), full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #num, Gleichung PD-Nrtl-Pia
    @staticmethod
    def PD_NRTL_pia(x, t, g1, g2):
        a = _Alpha
        u = g1
        U = h.np.exp((a*u)/(_r*t))
        v = g2
        V = h.np.exp((a*v)/(_r*t))

        _func = ((x/(1-x))-1)*(t/_h0RS)*(((_r*t)/x)+(2*(1-x)**2)*(((u*(U-1))/(x*(U-1)+1)**3)-((v*V*(V-1))/((x-1)*V-x)**3))-(2*(1-x))*((u/(x*(U-1)+1)**2)+((v*V)/(x-(x-1)*V)**2)))
        return _func

    #num, löst PD-Nrtl-Pia
    @staticmethod
    def PD_NRTL_pia_num(t, x, g):
        _tinitial = 393.35

        _initial = [_tinitial]
        _func = h.spi.solve_ivp(Diagrams.PD_NRTL_pia, [0.5, 1.0], _initial, method='RK45',
                                args=(g), t_eval=x, dense_output=True)
        return _func

    #num, soll PD-NRTL optimieren und g bestimmen
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

    #ana, lässt t berechnen
    @staticmethod
    def Bilanz_A_porter_pia(t, x, a):
        A = [a[0], a[1]]
        B = [a[2], a[3]]
        ya = Diagrams.y_porter(x, t, A)
        yb = Diagrams.y_porter(x, t, B)
        _func = (((1/2)*(((-_h0S/(_r*t))*(1-(t/_t0S)))+(((-_h0R/(_r*t))*(1-(t/_t0R))))))/((h.np.log(0.25)-h.np.log(x*ya*(1-x)*yb))))-1
        return _func

    #soll A Porter optimieren und as bestimmen
    @staticmethod
    def Bilanz_A_porter_fit_minfqs(_a, _xi, _texp, _steps, ):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_A_porter_pia, _texp[x], args=(_xi[x], _a,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #ana, lässt t berechnen
    @staticmethod
    def Bilanz_A_nrtl_pia(t, x, g):
        A = h.np.array([g[0], g[1]])
        B = h.np.array([g[2], g[3]])
        ya = Diagrams.y_nrtl(A, x, t, _Alpha)
        yb = Diagrams.y_nrtl(B, x, t, _Alpha)
        _func = (((1 / 2) * (((-_h0S / (_r * t)) * (1 - (t / _t0S))) + (((-_h0R / (_r * t)) * (1 - (t / _t0R)))))) / (
         (h.np.log(0.25) - h.np.log(x * ya * (1 - x) * yb)))) - 1
        return _func

    #soll nrtl koeffizienten bestimmen
    @staticmethod
    def Bilanz_A_nrtl_fit_minfqs(_g, _xi, _texp, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_A_nrtl_pia, _texp[x], args=(_xi[x], _g,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #berechnet t für Ansatz B
    @staticmethod
    def Bilanz_B_porter_pia(t, x, a):
        A = h.np.array([a[0], a[1]])
        B = h.np.array([a[2], a[3]])
        ya = Diagrams.y_porter(x, t, A)
        yb = Diagrams.y_porter(x, t, B)

        _func = (((1 / 2) * (((-_h0S / (_r * t)) * (1 - (t / _t0S))) + (((-_h0R / (_r * t)) * (1 - (t / _t0R)))))) / (
         (h.np.log(4 * x * ya * (1 - x) * yb)))) - 1
        return _func

    #soll Porterkoeff für Ansatz B optimieren
    @staticmethod
    def Bilanz_B_porter_fit_minfqs(_a, _xi, _texp, _steps, ):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_B_porter_pia, _texp[x], args=(_xi[x], _a,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #soll t für Ansatz B nrtl berechnen
    @staticmethod
    def Bilanz_B_nrtl_pia(t, x, g):
        A = h.np.array([g[0], g[1]])
        B = h.np.array([g[2], g[3]])
        ya = Diagrams.y_nrtl(A, x, t, _Alpha)
        yb = Diagrams.y_nrtl(B, x, t, _Alpha)
        _func = (((1 / 2) * (((-_h0S / (_r * t)) * (1 - (t / _t0S))) + (((-_h0R / (_r * t)) * (1 - (t / _t0R)))))) / ((h.np.log(4 * x * ya * (1 - x) * yb)))) - 1
        return _func

    #soll für ansatz B g's optimieren
    @staticmethod
    def Bilanz_B_nrtl_fit_minfqs(_g, _xi, _texp, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_B_nrtl_pia, _texp[x], args=(_xi[x], _g,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    #soll t für Bilanz C berechnen
    @staticmethod
    def Bilanz_C_porter_pia(t, x, a):
        _func = (_r*t)*((x/(1-x))-1)*((1/x)-2*(1-x)*(a[0]+(a[1]/t)))
        return _func

    @staticmethod
    def Bilanz_C_porter_fit_minfqs(_a, _xi, _texp, _steps, ):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_C_porter_pia, _texp[x], args=(_xi[x], _a,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe


    @staticmethod
    def Bilanz_C_nrtl_pia(x, t, g,):
        a = _Alpha
        u = g[0]
        U = h.np.exp((a * u) / (_r * t))
        v = g[1]
        V = h.np.exp((a * v) / (_r * t))

        _func = ((x / (1 - x)) - 1) * (t / _h0RS) * (((_r * t) / x) + (2 * (1 - x) ** 2) * (
                    ((u * (U - 1)) / (x * (U - 1) + 1) ** 3) - ((v * V * (V - 1)) / ((x - 1) * V - x) ** 3)) - (
                                                                 2 * (1 - x)) * ((u / (x * (U - 1) + 1) ** 2) + (
                    (v * V) / (x - (x - 1) * V) ** 2)))
        return _func

    @staticmethod
    def Bilanz_C_nrtl_fit_minfqs(_g, _xi, _texp, _steps,):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(Diagrams.Bilanz_C_nrtl_pia, _texp[x], args=(_xi[x], _g,),
                                   full_output=True)
            _tcalc[x] = _tcalc1[0]
            _tdiff[x] = _texp[x] - _tcalc[x]
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    @staticmethod
    def PD_Gleichungen_literatur():
        t_PD_ideal_links = h.np.zeros(len(TEXP_links))
        t_PD_ideal_links_diff = h.np.zeros(len(TEXP_links))

        t_PD_ideal_rechts = h.np.zeros(len(TEXP_rechts))
        t_PD_ideal_rechts_diff = h.np.zeros(len(TEXP_rechts))
        #print(TEXP_links)

        for x in range(len(TEXP_links)):
          t_PD_ideal_links[x] = spo.fsolve(Diagrams.PD_ideal, TEXP_links[x], args=XEXP_korr_links[x])
          t_PD_ideal_links_diff[x] = abs(TEXP_links[x]-t_PD_ideal_links[x])

        t_diff_norm_links = h.np.abs(h.np.divide(t_PD_ideal_links_diff, TEXP_links))
        ard_neu_norm_links = (100 / len(TEXP_links)) * sum(t_diff_norm_links)
        print('ARD_normiert für PD_ideal_links [%] =', ard_neu_norm_links)

        for x in range(len(TEXP_rechts)):
          t_PD_ideal_rechts[x] = spo.fsolve(Diagrams.PD_ideal, TEXP_rechts[x], args=XEXP_korr_rechts[x])
          t_PD_ideal_rechts_diff[x] = abs(TEXP_rechts[x]-t_PD_ideal_rechts[x])

        t_diff_norm_rechts = h.np.abs(h.np.divide(t_PD_ideal_rechts_diff, TEXP_rechts))
        ard_neu_norm_rechts = (100 / len(TEXP_rechts)) * sum(t_diff_norm_rechts)
        print('ARD_normiert für PD_ideal_rechts [%] =', ard_neu_norm_rechts)
        #print(t_PD_ideal_links)

        #steps_t_links = len(TEXP_links)
        #c = h.np.array([1])
        #res_1 = spo.minimize(Diagrams.PD_real_literatur_minfqs, c, args=(XEXP_korr_links, TEXP_links, steps_t_links), method='Nelder-Mead',)
        #print('deltag_AB s = ' + str(res_1.x[0]))

        return 0

    @staticmethod
    def PD_Gleichungen():

        t_Porter_links_Pia = h.np.zeros(len(TEXP_links))
        t_Porter_links_diff_Pia = h.np.zeros(len(TEXP_links))
        t_NRTL_links_Pia = h.np.zeros(len(TEXP_links))
        t_NRTL_links_diff_Pia = h.np.zeros(len(TEXP_links))

        t_Porter_rechts_Pia = h.np.zeros(len(TEXP_rechts))
        t_Porter_rechts_diff_Pia = h.np.zeros(len(TEXP_rechts))
        t_NRTL_rechts_Pia = h.np.zeros(len(TEXP_rechts))
        t_NRTL_rechts_diff_Pia = h.np.zeros(len(TEXP_rechts))

        _tinitial = 393.35

        _initial = [_tinitial]
        _tcalcRSS_Porter = h.spi.solve_ivp(Diagrams.PD_Porter_pia, [0.000000001, 0.5], _initial, method='RK45',args=(_aSa), t_eval=XEXP_korr_links, dense_output=True)
        _tcalcRSR_Porter = h.spi.solve_ivp(Diagrams.PD_Porter_pia, [0.5, 1.0], _initial, method='RK45',
                                    args=(_aRa), t_eval=XEXP_rechts, dense_output=True)

        _tcalcRSS_NRTL = h.spi.solve_ivp(Diagrams.PD_NRTL_pia, [0.000000001, 0.5], _initial, method='RK45',
                                           args=(_gSa), t_eval=XEXP_korr_links, dense_output=True)
        _tcalcRSR_NRTL = h.spi.solve_ivp(Diagrams.PD_NRTL_pia, [0.5, 1.0], _initial, method='RK45',
                                           args=(_gRa), t_eval=XEXP_rechts, dense_output=True)

        t_Porter_links = h.np.reshape(_tcalcRSS_Porter.y, len(TEXP_links))
        t_Porter_rechts = h.np.reshape(_tcalcRSR_Porter.y, len(TEXP_rechts)-1)
        t_NRTL_links = h.np.reshape(_tcalcRSS_NRTL.y, len(TEXP_links))
        t_NRTL_rechts = h.np.reshape(_tcalcRSR_NRTL.y, len(TEXP_rechts) - 1)

        print(t_Porter_links, t_NRTL_links)

        for x in range(len(TEXP_links)):
          t_Porter_links_diff_Pia[x] = abs(TEXP_links[x]-t_Porter_links[x])
          t_NRTL_links_diff_Pia[x] = abs(TEXP_links[x] - t_NRTL_links[x])

        #print(t_Porter_links_diff_Pia)
        t_diff_norm_links_P = h.np.abs(h.np.divide(t_Porter_links_diff_Pia, TEXP_links))
        ard_neu_norm_links_P = (100 / len(TEXP_links)) * sum(t_diff_norm_links_P)
        t_diff_norm_links_N = h.np.abs(h.np.divide(t_NRTL_links_diff_Pia, TEXP_links))
        ard_neu_norm_links_N = (100 / len(TEXP_links)) * sum(t_diff_norm_links_N)


        for x in range(len(TEXP_rechts)-1):
          t_Porter_rechts_diff_Pia[x] = abs(TEXP_rechts[x]-t_Porter_rechts[x])
          t_NRTL_rechts_diff_Pia[x] = abs(TEXP_rechts[x] - t_NRTL_rechts[x])

        t_diff_norm_rechts_P = h.np.abs(h.np.divide(t_Porter_rechts_diff_Pia, TEXP_rechts))
        ard_neu_norm_rechts_P = (100 / len(TEXP_rechts)) * sum(t_diff_norm_rechts_P)
        t_diff_norm_rechts_N = h.np.abs(h.np.divide(t_NRTL_rechts_diff_Pia, TEXP_rechts))
        ard_neu_norm_rechts_N = (100 / len(TEXP_rechts)) * sum(t_diff_norm_rechts_N)

        print('ARD_normiert für PD_Porter_links [%] =', ard_neu_norm_links_P)
        print('ARD_normiert für PD_Porter_rechts [%] =', ard_neu_norm_rechts_P)
        print('ARD_normiert für PD_NRTL_links [%] =', ard_neu_norm_links_N)
        print('ARD_normiert für PD_NRTL_rechts [%] =', ard_neu_norm_rechts_N)

        return 0

    @staticmethod
    def Ansatz_A():
        t_Porter_A_Pia = h.np.zeros(len(TEXP))
        t_Porter_A_diff_Pia = h.np.zeros(len(TEXP))
        t_NRTL_A_Pia = h.np.zeros(len(TEXP))
        t_NRTL_A_diff_Pia = h.np.zeros(len(TEXP))

        t_Porter_A = h.np.zeros(len(TEXP))
        t_Porter_A_diff = h.np.zeros(len(TEXP))
        t_NRTL_A = h.np.zeros(len(TEXP))
        t_NRTL_A_diff = h.np.zeros(len(TEXP))

        for x in range(len(TEXP)):
            t_Porter_A_Pia[x] = spo.fsolve(Diagrams.Bilanz_A_porter_pia, TEXP[x], args=(XEXP[x], _aGes))
            t_Porter_A_diff_Pia[x] = abs(TEXP[x] - t_Porter_A_Pia[x])
            t_NRTL_A_Pia[x] = spo.fsolve(Diagrams.Bilanz_A_nrtl_pia, TEXP[x], args=(XEXP[x], _gGes))
            t_NRTL_A_diff_Pia[x] = abs(TEXP[x] - t_NRTL_A_Pia[x])
            #print(t_Porter_A_Pia[x])

        t_diff_norm_P_Pia = h.np.abs(h.np.divide(t_Porter_A_diff_Pia, TEXP))
        ard_neu_norm_P_Pia = (100 / len(TEXP)) * sum(t_diff_norm_P_Pia)
        t_diff_norm_N_Pia = h.np.abs(h.np.divide(t_NRTL_A_diff_Pia, TEXP))
        ard_neu_norm_N_Pia = (100 / len(TEXP)) * sum(t_diff_norm_N_Pia)
        print('ARD_normiert für A_Porter_Pia [%] =', ard_neu_norm_P_Pia)
        print('ARD_normiert für A_NRTL_Pia [%] =', ard_neu_norm_N_Pia)

        steps_t = len(XEXP)
        res_Porter = spo.minimize(Diagrams.Bilanz_A_porter_fit_minfqs, _aGes, args=(XEXP, TEXP, steps_t,),
                                        method='Nelder-Mead',)
        print('A_A1 = ' + str(res_Porter.x[0]))
        print('A_A2 = ' + str(res_Porter.x[1]))
        print('A_B1 = ' + str(res_Porter.x[2]))
        print('A_B1 = ' + str(res_Porter.x[3]))
        _aNeu = (res_Porter.x[0], res_Porter.x[1], res_Porter.x[2], res_Porter.x[3])

        res_NRTL = spo.minimize(Diagrams.Bilanz_A_nrtl_fit_minfqs, _gGes, args=(XEXP, TEXP, steps_t,),
                                  method='Nelder-Mead', )
        print('A_gab = ' + str(res_NRTL.x[0]))
        print('A_gba = ' + str(res_NRTL.x[1]))
        print('A_gAB = ' + str(res_NRTL.x[2]))
        print('A_gBA = ' + str(res_NRTL.x[3]))
        _gNeu = (res_NRTL.x[0], res_NRTL.x[1], res_NRTL.x[2], res_NRTL.x[3])

        t_diff_P = h.np.zeros(steps_t)
        _tcalc_P = h.np.zeros(steps_t)
        _gamma_P = h.np.zeros(steps_t)

        t_diff_N = h.np.zeros(steps_t)
        _tcalc_N = h.np.zeros(steps_t)
        _gamma_N = h.np.zeros(steps_t)

        for x in range(len(TEXP)):
            t_Porter_A[x] = spo.fsolve(Diagrams.Bilanz_A_porter_pia, TEXP[x], args=(XEXP[x], _aNeu))
            t_Porter_A_diff[x] = abs(TEXP[x] - t_Porter_A[x])
            t_NRTL_A[x] = spo.fsolve(Diagrams.Bilanz_A_nrtl_pia, TEXP[x], args=(XEXP[x], _gNeu))
            t_NRTL_A_diff[x] = abs(TEXP[x] - t_NRTL_A[x])

        t_diff_norm_P = h.np.abs(h.np.divide(t_Porter_A_diff, TEXP))
        ard_neu_norm_P = (100 / len(TEXP_links)) * sum(t_diff_norm_P)
        t_diff_norm_N = h.np.abs(h.np.divide(t_NRTL_A_diff, TEXP))
        ard_neu_norm_N = (100 / len(TEXP_links)) * sum(t_diff_norm_N)
        print('ARD_normiert für A_Porter [%] =', ard_neu_norm_P)
        print('ARD_normiert für A_NRTL [%] =', ard_neu_norm_N)
        return 0

    @staticmethod
    def Ansatz_B():
        t_Porter_B_Pia = h.np.zeros(len(TEXP))
        t_Porter_B_diff_Pia = h.np.zeros(len(TEXP))
        t_NRTL_B_Pia = h.np.zeros(len(TEXP))
        t_NRTL_B_diff_Pia = h.np.zeros(len(TEXP))

        t_Porter_B = h.np.zeros(len(TEXP))
        t_Porter_B_diff = h.np.zeros(len(TEXP))
        t_NRTL_B = h.np.zeros(len(TEXP))
        t_NRTL_B_diff = h.np.zeros(len(TEXP))

        for x in range(len(TEXP)):
            t_Porter_B_Pia[x] = spo.fsolve(Diagrams.Bilanz_B_porter_pia, TEXP[x], args=(XEXP[x], _aGes))
            t_Porter_B_diff_Pia[x] = abs(TEXP[x] - t_Porter_B_Pia[x])
            t_NRTL_B_Pia[x] = spo.fsolve(Diagrams.Bilanz_B_nrtl_pia, TEXP[x], args=(XEXP[x], _gGes))
            t_NRTL_B_diff_Pia[x] = abs(TEXP[x] - t_NRTL_B_Pia[x])

        t_diff_norm_P_Pia = h.np.abs(h.np.divide(t_Porter_B_diff_Pia, TEXP))
        ard_neu_norm_P_Pia = (100 / len(TEXP)) * sum(t_diff_norm_P_Pia)
        t_diff_norm_N_Pia = h.np.abs(h.np.divide(t_NRTL_B_diff_Pia, TEXP))
        ard_neu_norm_N_Pia = (100 / len(TEXP)) * sum(t_diff_norm_N_Pia)
        print('ARD_normiert für B_Porter_Pia [%] =', ard_neu_norm_P_Pia)
        print('ARD_normiert für B_NRTL_Pia [%] =', ard_neu_norm_N_Pia)


        steps_t = len(XEXP)
        res_Porter = spo.minimize(Diagrams.Bilanz_B_porter_fit_minfqs, _aGes, args=(XEXP, TEXP, steps_t,),
                                  method='Powell', )
        print('B_A1 = ' + str(res_Porter.x[0]))
        print('B_A2 = ' + str(res_Porter.x[1]))
        print('B_B1 = ' + str(res_Porter.x[2]))
        print('B_B1 = ' + str(res_Porter.x[3]))
        _aNeu = (res_Porter.x[0], res_Porter.x[1], res_Porter.x[2], res_Porter.x[3])

        res_NRTL = spo.minimize(Diagrams.Bilanz_B_nrtl_fit_minfqs, _gGes, args=(XEXP, TEXP, steps_t,),
                                method='Powell', )
        print('B_gab = ' + str(res_NRTL.x[0]))
        print('B_gba = ' + str(res_NRTL.x[1]))
        print('B_gAB = ' + str(res_NRTL.x[2]))
        print('B_gBA = ' + str(res_NRTL.x[3]))
        _gNeu = (res_NRTL.x[0], res_NRTL.x[1], res_NRTL.x[2], res_NRTL.x[3])

        for i in range(len(TEXP)):
            t_Porter_B[i] = spo.fsolve(Diagrams.Bilanz_B_porter_pia, TEXP[i], args=(XEXP[i], _aNeu))
            t_Porter_B_diff[i] = abs(TEXP[i] - t_Porter_B[i])
            t_NRTL_B[i] = spo.fsolve(Diagrams.Bilanz_B_nrtl_pia, TEXP[i], args=(XEXP[i], _gNeu))
            t_NRTL_B_diff[i] = abs(TEXP[i] - t_NRTL_B[i])

        t_diff_norm_P = h.np.abs(h.np.divide(t_Porter_B_diff, TEXP))
        ard_neu_norm_P = (100 / len(TEXP)) * sum(t_diff_norm_P)
        t_diff_norm_N = h.np.abs(h.np.divide(t_NRTL_B_diff, TEXP))
        ard_neu_norm_N = (100 / len(TEXP)) * sum(t_diff_norm_N)
        print('ARD_normiert für B_Porter [%] =', ard_neu_norm_P)
        print('ARD_normiert für B_NRTL [%] =', ard_neu_norm_N)
        return 0

    @staticmethod
    def Ansatz_C():
        t_Porter_C_Pia_links = h.np.zeros(len(TEXP_links))
        t_Porter_C_diff_Pia_links = h.np.zeros(len(TEXP_links))
        t_Porter_C_Pia_rechts = h.np.zeros(len(TEXP_rechts))
        t_Porter_C_diff_Pia_rechts = h.np.zeros(len(TEXP_rechts))

        t_NRTL_C_Pia_links = h.np.zeros(len(TEXP_links))
        t_NRTL_C_diff_Pia_links = h.np.zeros(len(TEXP_links))
        t_NRTL_C_Pia_rechts= h.np.zeros(len(TEXP_rechts))
        t_NRTL_C_diff_Pia_rechts = h.np.zeros(len(TEXP_rechts))

        t_Porter_C_links = h.np.zeros(len(TEXP_links))
        t_Porter_C_diff_links = h.np.zeros(len(TEXP_links))
        t_Porter_C_rechts = h.np.zeros(len(TEXP_rechts))
        t_Porter_C_diff_rechts = h.np.zeros(len(TEXP_rechts))

        t_NRTL_C_links = h.np.zeros(len(TEXP_links))
        t_NRTL_C_diff_links = h.np.zeros(len(TEXP_links))
        t_NRTL_C_rechts = h.np.zeros(len(TEXP_rechts))
        t_NRTL_C_diff_rechts = h.np.zeros(len(TEXP_rechts))

        for x in range(len(TEXP_links)):
            t_Porter_C_Pia_links[x] = spo.fsolve(Diagrams.Bilanz_C_porter_pia, TEXP_links[x], args=(XEXP_korr_links[x], _aSa))
            t_Porter_C_diff_Pia_links[x] = abs(TEXP_links[x] - t_Porter_C_Pia_links[x])
            t_NRTL_C_Pia_links[x] = spo.fsolve(Diagrams.Bilanz_C_nrtl_pia, TEXP_links[x], args=(XEXP_korr_links[x], _gSa))
            t_NRTL_C_diff_Pia_links[x] = abs(TEXP_links[x] - t_NRTL_C_Pia_links[x])

        for x in range(len(TEXP_rechts)):
            t_Porter_C_Pia_rechts[x] = spo.fsolve(Diagrams.Bilanz_C_porter_pia, TEXP_rechts[x],
                                                 args=(XEXP_korr_rechts[x], _aSa))
            t_Porter_C_diff_Pia_rechts[x] = abs(TEXP_rechts[x] - t_Porter_C_Pia_rechts[x])
            t_NRTL_C_Pia_rechts[x] = spo.fsolve(Diagrams.Bilanz_C_nrtl_pia, TEXP_rechts[x],
                                               args=(XEXP_korr_rechts[x], _gSa))
            t_NRTL_C_diff_Pia_rechts[x] = abs(TEXP_rechts[x] - t_NRTL_C_Pia_rechts[x])

        t_diff_norm_P_Pia_links = h.np.abs(h.np.divide(t_Porter_C_diff_Pia_links, TEXP_links))
        ard_neu_norm_P_Pia_links = (100 / len(TEXP_links)) * sum(t_diff_norm_P_Pia_links)
        t_diff_norm_P_Pia_rechts = h.np.abs(h.np.divide(t_Porter_C_diff_Pia_rechts, TEXP_rechts))
        ard_neu_norm_P_Pia_rechts = (100 / len(TEXP_rechts)) * sum(t_diff_norm_P_Pia_rechts)

        t_diff_norm_N_Pia_links = h.np.abs(h.np.divide(t_NRTL_C_diff_Pia_links, TEXP_links))
        ard_neu_norm_N_Pia_links = (100 / len(TEXP_links)) * sum(t_diff_norm_N_Pia_links)
        t_diff_norm_N_Pia_rechts = h.np.abs(h.np.divide(t_NRTL_C_diff_Pia_rechts, TEXP_rechts))
        ard_neu_norm_N_Pia_rechts = (100 / len(TEXP_links)) * sum(t_diff_norm_N_Pia_rechts)
        print('ARD_normiert für C_Porter_Pia_links [%] =', ard_neu_norm_P_Pia_links)
        print('ARD_normiert für C_Porter_Pia_rechts [%] =', ard_neu_norm_P_Pia_links)
        print('ARD_normiert für C_NRTL_Pia_links [%] =', ard_neu_norm_N_Pia_links)
        print('ARD_normiert für C_NRTL_Pia_rechts [%] =', ard_neu_norm_N_Pia_rechts)

        steps_t = len(TEXP_links)
        _aSa2 = h.np.array([-10, 4600])
        res_Porter = spo.minimize(Diagrams.Bilanz_C_porter_fit_minfqs, _aSa2, args=(XEXP_korr_links_Porter, TEXP_links, steps_t,),
                                  method='Powell', )
        print('C_A1 = ' + str(res_Porter.x[0]))
        print('C_A2 = ' + str(res_Porter.x[1]))

        _aNeu = (res_Porter.x[0], res_Porter.x[1])
        res_NRTL = spo.minimize(Diagrams.Bilanz_C_nrtl_fit_minfqs, _gSa, args=(XEXP_korr_links, TEXP_links, steps_t,),
                                method='Powell', )
        print('C_gab = ' + str(res_NRTL.x[0]))
        print('C_gba = ' + str(res_NRTL.x[1]))
        #print('C_gAB = ' + str(res_NRTL.x[2]))
        #print('C_gBA = ' + str(res_NRTL.x[3]))
        _gNeu = (res_NRTL.x[0], res_NRTL.x[1],)

        for i in range(len(TEXP_links)):
            t_Porter_C[i] = spo.fsolve(Diagrams.Bilanz_C_porter_pia, TEXP_links[i], args=(XEXP_korr_links[i], _aNeu))
            t_Porter_C_diff[i] = abs(TEXP_links[i] - t_Porter_C[i])
            t_NRTL_C[i] = spo.fsolve(Diagrams.Bilanz_C_nrtl_pia, TEXP_links[i], args=(XEXP_korr_links[i], _gNeu))
            t_NRTL_C_diff[i] = abs(TEXP_links[i] - t_NRTL_C[i])

        t_diff_norm_P = h.np.abs(h.np.divide(t_Porter_C_diff, TEXP_links))
        ard_neu_norm_P = (100 / len(TEXP_links)) * sum(t_diff_norm_P)
        t_diff_norm_N = h.np.abs(h.np.divide(t_NRTL_C_diff, TEXP_links))
        ard_neu_norm_N = (100 / len(TEXP_links)) * sum(t_diff_norm_N)
        print('ARD_normiert für C_Porter [%] =', ard_neu_norm_P)
        print('ARD_normiert für C_NRTL [%] =', ard_neu_norm_N)
        print(t_Porter_C)
        return 0