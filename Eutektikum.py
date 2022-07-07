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


class EutFind:
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
    def t_sle(_texp, _xi, _h0, _t0, _g,  _alpha):

        _yi = EutFind.y_nrtl(_g, _xi, _texp, _alpha)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx)-1
        return _func


    @staticmethod
    def t_sle_ideal(_texp, _xi, _h0, _t0,):
        _yi = 1
        _gx = _yi * (_xi)
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx) - 1
        return _func

    @staticmethod
    def t_sle_rac_ideal(_t, _x, _h0, _t0):
        _xkorr = 4*_x*(1-_x)
        _func = ((_h0/_r)*((1/_t0)-(1/_t)))-h.np.log(_xkorr)

        return _func

    @staticmethod
    def t_sle_rac_nrtl(_t, _x, _h0rac, _t0rac, _gi, _gj, _alpha ):

        _yi = EutFind.y_nrtl(_gi, _x, _t, _alpha)
        _yj = EutFind.y_nrtl(_gj, _x, _t, _alpha)

        _xkorr2 = 4*_x*(1-_x)
        _gx = h.np.exp(_yi) * h.np.exp(_yj) * _xkorr2

        _func = (h.np.exp((-_h0rac / (_r * _t)) * (1 - (_t / _t0rac))) / _gx) - 1

        return _func

    @staticmethod
    def x_sle_ideal(_texp, _h0, _t0, _g, _alpha):
        _func = _yi = 1
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))))
        return _func

    @staticmethod
    def calc_eut_ideal(_p, _h0i, _t0i, _h0j, _t0j):
        _t, _x = _p

        _pure = EutFind.t_sle_ideal(_t, _x, _h0i, _t0i)
        _rac = EutFind.t_sle_rac_ideal(_t, _x, _h0j, _t0j)

        _func = [_pure, _rac]
        return _func

    @staticmethod
    def calc_eut_nrtl2(_p, _h0i, _t0i, _gi, _gj, _h0rac, _t0rac, _alpha):
        _t, _x = _p

        _pure = EutFind.t_sle(_t, _x, _h0i, _t0i, _gi, _alpha)
        _rac = EutFind.t_sle_rac_nrtl(_t, _x, _h0rac, _t0rac, _gi, _gj, _alpha)

        _func = [_pure, _rac]
        return _func

    @staticmethod
    def calc_eut_nrtl1(_p, _h0i, _t0i, _gi, _h0j, _t0j, _alpha):
        _t, _x = _p

        _pure = EutFind.t_sle(_t, _x, _h0i, _t0i, _gi, _alpha)
        _rac = EutFind.t_sle_rac_ideal(_t, _x, _h0j, _t0j)

        _func = [_pure, _rac]
        return _func

    @staticmethod
    def num_equation123(_x, _t, _alpha, _gSR, _gRS, _h0):

        _T = _t[0]
        #_func = (_T / _h0) * (_x / (1 - _x)) * (((_r * _T) / _x) + ((2 * h.np.exp(-_alpha * _gSR / _r * _T) * h.np.exp(
        #    -_alpha * _gSR / _r * _T) * _gSR * (h.np.exp(-_alpha * _gSR / _r * _T) * _x - h.np.exp(
        #    -_alpha * _gSR / _r * _T) + 1)) / ((h.np.exp(-_alpha * _gSR / _r * _T) - 1) * _x - h.np.exp(
        #    -_alpha * _gSR / _r * _T)) ** 3)((2 * _gRS * (_x + h.np.exp(-_alpha * _gRS / _r*_T) - 1)) / ((h.np.exp(-_alpha * _gRS / _r*_T) - 1) * _x + 1) ** 3))


        _func1 = (_t/_h0)*(-(_x/(1-_x))+1)*(((_r*_t)/_x)+((2*(h.np.exp((-_alpha*_gSR)/(_r*_t)))**2)*(h.np.exp((-_alpha*_gSR)/(_r*_t))*_x-h.np.exp((-_alpha*_gSR)/(_r*_t))+1)/(((h.np.exp((-_alpha*_gSR)/(_r*_t)))-1)*_x+1)**3)-((2*_gRS*(_x+h.np.exp((-_alpha*_gRS)/(_r*_t)))-1)/((h.np.exp((-_alpha*_gRS)/(_r*_t)))*_x+1)**3))
        return _func1

    @staticmethod
    def num_equation1234(_x, _t, _alpha, _gSR, _gRS, _h0):

        _T = _t[0]
        # _func = (_T / _h0) * (_x / (1 - _x)) * (((_r * _T) / _x) + ((2 * h.np.exp(-_alpha * _gSR / _r * _T) * h.np.exp(
        #    -_alpha * _gSR / _r * _T) * _gSR * (h.np.exp(-_alpha * _gSR / _r * _T) * _x - h.np.exp(
        #    -_alpha * _gSR / _r * _T) + 1)) / ((h.np.exp(-_alpha * _gSR / _r * _T) - 1) * _x - h.np.exp(
        #    -_alpha * _gSR / _r * _T)) ** 3)((2 * _gRS * (_x + h.np.exp(-_alpha * _gRS / _r*_T) - 1)) / ((h.np.exp(-_alpha * _gRS / _r*_T) - 1) * _x + 1) ** 3))

        _func1 = (_t / _h0) * ((_x / (1 - _x)) - 1) * (((_r * _t) / _x) + (
                    (2 * (h.np.exp((-_alpha * _gSR) / (_r * _t))) ** 2) * (
                        h.np.exp((-_alpha * _gSR) / (_r * _t)) * _x - h.np.exp((-_alpha * _gSR) / (_r * _t)) + 1) / (
                                ((h.np.exp((-_alpha * _gSR) / (_r * _t))) - 1) * _x + 1) ** 3) - ((2 * _gRS * (
                    _x + h.np.exp((-_alpha * _gRS) / (_r * _t))) - 1) / ((h.np.exp(
            (-_alpha * _gRS) / (_r * _t))) * _x + 1) ** 3))
        return _func1

    @staticmethod
    def solve_num_eq(_alpha, _gSR, _gRS, _h0, _teval):

        _tinitial = int(393.35)
        _y = h.np.empty(_tinitial)

        _y0 = _y[0]
        _func = h.spi.solve_ivp(h.Eut.EutFind.num_equation, [0.5, 0], [393.35, 0], args=(_alpha, _gSR, _gRS, _h0), t_eval=_teval,)

        return _func

    @staticmethod
    def pure_ma():

            _gabS = -1403.640
            _gbaS = 2972.166
            _gSa = h.np.array([_gabS, _gbaS, ])

            _gabRS = -3143.16
            _gbaRS = 9443.10
            _gRSa = h.np.array([_gabRS, _gbaRS, ])

            _gabR = -23.289
            _gbaR = 7526.385
            _gRa = h.np.array([_gabR, _gbaR, ])

            def _ideal_sle():
                _xin = h.np.zeros(100)
                _xinRac = h.np.zeros(100)
                _xinRac[0] = 0.5
                _tin = h.np.zeros(100)
                _tin[0] = 273.15

                _tcalcS = h.np.zeros(100)
                _tcalcRSS = h.np.zeros(100)
                _tcalcRSR = h.np.zeros(100)
                _tcalcR = h.np.zeros(100)

                _tTest = h.np.zeros(100)

                # loading up values
                for x in range(len(_xin)):
                    _xin[x] = _xin[x-1]+0.01
                    _xinRac[x] = _xinRac[0]+.005*x

                for x in range(len(_tin)):
                    _tin[x] = _tin[0]+(2.0*x)

                # loading up calc points
                for x in range(len(_xin)):
                    _tcalcS[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0S, _t0S,))
                    _tcalcRSS[x] = h.spo.fsolve(EutFind.t_sle_rac_ideal, _tin[x], args=(_xinRac[x], _h0RS, _t0RS,))
                    _tcalcRSR[x] = h.spo.fsolve(EutFind.t_sle_rac_ideal, _tin[x], args=(_xinRac[x], _h0RS, _t0RS,))
                    _tcalcR[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0R, _t0R,))

                # plotting
                figure, axis = plt.subplots(2, constrained_layout=True)

                axis[0].plot(_xin, _tcalcS, '-g', label='S-Ma-Ideal-SLE')
                axis[0].plot(_xinRac, _tcalcRSS, '--g', label='Rac-Ma-Prigogine')
                axis[0].set_title("S-Rac-ideal")
                axis[0].set_ylabel('Temperatur / [K]')
                axis[0].set_xlabel('x-S-Ma / [-]')
                axis[0].legend()

                line_1 = LineString(h.np.column_stack((_xin, _tcalcS)))
                line_2 = LineString(h.np.column_stack((_xinRac, _tcalcRSS)))
                intersection = line_1.intersection(line_2)

                axis[0].plot(*intersection.xy, 'ro')
                x, y = intersection.xy
                xout0 = x[0]
                xoutstr0 = "{:10.4f}".format(xout0)
                yout0 = y[0]
                youtstr0 = "{:10.4f}".format(yout0)
                stringout0 = str(f"SP({xoutstr0}"f"\\{youtstr0})")
                axis[0].text(x[0], y[0], 'SP')
                axis[0].text(x[0], 250, stringout0)
                print(x, y,)

                axis[1].plot(_xin, _tcalcR, '-b', label='R-Ma-Ideal-SLE')
                axis[1].plot(_xinRac, _tcalcRSR, '--b', label='Rac-Ma-Prigogine')
                axis[1].set_title("R-Rac-ideal")
                axis[1].set_ylabel('Temperatur / [K]')
                axis[1].set_xlabel('x-R-Ma / [-]')
                axis[1].legend()

                line_3 = LineString(h.np.column_stack((_xin, _tcalcR)))
                line_4 = LineString(h.np.column_stack((_xinRac, _tcalcRSR)))
                intersection = line_3.intersection(line_4)

                axis[1].plot(*intersection.xy, 'ro')
                x2, y2 = intersection.xy
                xout1 = x2[0]
                xoutstr1 = "{:10.4f}".format(xout1)
                yout1 = y2[0]
                youtstr1 = "{:10.4f}".format(yout1)
                stringout1 = str(f"SP({xoutstr1}"f"\\{youtstr1})")
                axis[1].text(x2[0], y2[0], 'SP')
                axis[1].text(x2[0], 250, stringout1)
                print(x2, y2)

                plt.show()

                return 0

            def _nrtl_pure_comp_sle():
                _xin = h.np.zeros(100)
                _xinRac = h.np.zeros(100)
                _xinRac[0] = 0.5
                _tin = h.np.zeros(100)
                _tin[0] = 273.15

                _tcalcS = h.np.zeros(100)
                _tcalcRSS = h.np.zeros(100)
                _tcalcRSR = h.np.zeros(100)
                _tcalcR = h.np.zeros(100)

                _tTest = h.np.zeros(100)

                _xStest = h.np.array([1-.0, 1-.0594, 1-.0998, 1-.1496, 1-.2003, 1-.2492, 1-.3098])
                _tStest = h.np.array([404.65, 400.65, 399.15, 397.15, 393.65, 391.75, 387.75])
                _tStestCalc = h.np.zeros(len(_tStest))
                # loading up values
                for x in range(len(_xin)):
                    _xin[x] = _xin[x - 1] + 0.01
                    _xinRac[x] = _xinRac[0]+.005*x

                for x in range(len(_tin)):
                    _tin[x] = _tin[0] + (2.0 * x)

                # loading up calc points
                for x in range(len(_xin)):
                    _tcalcS[x] = h.spo.fsolve(EutFind.t_sle, _tin[x], args=(_xin[x], _h0S, _t0S, _gSa, _Alpha))
                    _tcalcRSS[x] = h.spo.fsolve(EutFind.t_sle_rac_ideal, _tin[x], args=(_xinRac[x], _h0RS, _t0RS,))
                    _tcalcRSR[x] = h.spo.fsolve(EutFind.t_sle_rac_ideal, _tin[x], args=(_xinRac[x], _h0RS, _t0RS,))
                    _tcalcR[x] = h.spo.fsolve(EutFind.t_sle, _tin[x], args=(_xin[x], _h0R, _t0R, _gRa, _Alpha))

                for x in range(len(_xStest)):
                    _tStestCalc[x] = h.spo.fsolve(EutFind.t_sle, _tStest[x], args=(_xStest[x], _h0S, _t0S, _gSa, _Alpha))


                # plotting
                figure, axis = plt.subplots(2, constrained_layout=True)


                axis[0].plot(_xin, _tcalcS, '-g', label='S-Ma-NRTL')
                axis[0].plot(_xinRac, _tcalcRSS, '--g', label='Rac-Ma-Prigogine')
                axis[0].set_title("S-Rac-NRTL")
                axis[0].set_ylabel('Temperatur / [K]')
                axis[0].set_xlabel('x-S-Ma / [-]')
                axis[0].legend()

                line_1 = LineString(h.np.column_stack((_xin, _tcalcS)))
                line_2 = LineString(h.np.column_stack((_xinRac, _tcalcRSS)))
                intersection = line_1.intersection(line_2)

                axis[0].plot(*intersection.xy, 'ro')
                x, y = intersection.xy
                xout0 = x[0]
                xoutstr0 = "{:10.4f}".format(xout0)
                yout0 = y[0]
                youtstr0 = "{:10.4f}".format(yout0)
                stringout0 = str(f"SP({xoutstr0}"f"\\{youtstr0})")
                axis[0].text(x[0], y[0], 'SP')
                axis[0].text(x[0], 300, stringout0)
                print(x, y)

                axis[1].plot(_xin, _tcalcR, '-b', label='R-Ma-NRTL')
                axis[1].plot(_xinRac, _tcalcRSR, '--b', label='Rac-Ma-Prigogine')
                axis[1].set_title("R-Rac-NRTL")
                axis[1].set_ylabel('Temperatur / [K]')
                axis[1].set_xlabel('x-R-Ma / [-]')
                axis[1].legend()

                line_3 = LineString(h.np.column_stack((_xin, _tcalcR)))
                line_4 = LineString(h.np.column_stack((_xinRac, _tcalcRSR)))
                intersection = line_3.intersection(line_4)

                axis[1].plot(*intersection.xy, 'ro')
                x2, y2 = intersection.xy
                xout1 = x2[0]
                xoutstr1 = "{:10.4f}".format(xout1)
                yout1 = y2[0]
                youtstr1 = "{:10.4f}".format(yout1)
                stringout1 = str(f"SP({xoutstr1}"f"\\{youtstr1})")
                axis[1].text(x2[0], y2[0], 'SP')
                axis[1].text(x2[0], 300, stringout1)
                print(x2, y2)

                plt.show()

                return 0

            def _nrtl_nrtl():
                _xin = h.np.zeros(100)
                _xin[0] = 0
                _xinRac = h.np.zeros(100)
                _xinRac[0] = 0.5
                _tin = h.np.zeros(100)
                _tin[0] = 273.15

                _tcalcS = h.np.zeros(100)
                _tcalcRSS = h.np.zeros(100)
                _tcalcRSR = h.np.zeros(100)
                _tcalcR = h.np.zeros(100)

                _tTest = h.np.zeros(100)

                # loading up values
                for x in range(len(_xin)):
                    _xin[x] = _xin[x-1] + 0.01
                    _xinRac[x] = _xinRac[0]+0.005*x

                for x in range(len(_tin)):
                    _tin[x] = _tin[0] + (2.0 * x)

                # loading up calc points
                for x in range(len(_xin+1)):
                    _tcalcS[x] = h.spo.fsolve(EutFind.t_sle, _tin[x], args=(_xin[x], _h0S, _t0S, _gSa, _Alpha))
                    # _tcalcRSS[x] = h.spo.fsolve(EutFind.t_sle_rac_nrtl, _tin[x], args=(_xinRac[x], _h0RS, _t0RS,_gSa, _gRa, _Alpha))
                    # _tcalcRSR[x] = h.spo.fsolve(EutFind.t_sle_rac_nrtl, _tin[x], args=(_xinRac[x], _h0RS, _t0RS,_gSa, _gRa, _Alpha))
                    _tcalcR[x] = h.spo.fsolve(EutFind.t_sle, _tin[x], args=(_xin[x], _h0R, _t0R, _gRa, _Alpha))

                _tinitial = 393.35
                print(_xinRac)
                _initial = [_tinitial]
                _tcalcRSS = h.spi.solve_ivp(h.Eut.EutFind.num_equation123, [0.5, 1.0], _initial, method='RK45', args=(_Alpha, _gabS, _gbaS, _h0RS), t_eval=_xinRac, dense_output=True)
                _tcalcRSR = h.spi.solve_ivp(h.Eut.EutFind.num_equation1234, [0.5, 1.0], _initial, method='RK45', args=(_Alpha, _gabR, _gbaR, _h0RS), t_eval=_xinRac, dense_output=True)

                #print(_tcalcRSS)
                _tSLERSS = h.np.reshape(_tcalcRSS.y, 100)
                _tSLERSR = h.np.reshape(_tcalcRSR.y, 100)

                # plotting
                figure, axis = plt.subplots(2, constrained_layout=True)


                axis[0].plot(_xin, _tcalcS, '-g', label='S-Ma-NRTL')
                axis[0].plot(_xinRac, _tSLERSS, '--g', label='Rac-Ma-NRTL')
                axis[0].set_title("S-Rac-NRTL")
                axis[0].set_ylabel('Temperatur / [K]')
                axis[0].set_xlabel('x-S-Ma / [-]')
                axis[0].legend()

                line_1 = LineString(h.np.column_stack((_xin, _tcalcS)))
                line_2 = LineString(h.np.column_stack((_xinRac, _tSLERSS)))
                intersection = line_1.intersection(line_2)

                axis[0].plot(*intersection.xy, 'ro')
                x, y = intersection.xy
                xout0 = x[0]
                xoutstr0 = "{:10.4f}".format(xout0)
                yout0 = y[0]
                youtstr0 = "{:10.4f}".format(yout0)
                stringout0 = str(f"SP({xoutstr0}"f"\\{youtstr0})")
                axis[0].text(x[0], y[0], 'SP')
                axis[0].text(x[0], 300, stringout0)
                print(x, y)

                axis[1].plot(_xin, _tcalcR, '-b', label='R-Ma-NRTL')
                axis[1].plot(_xinRac, _tSLERSR, '--b', label='Rac-Ma-NRTL')
                axis[1].set_title("R-Rac-NRTL")
                axis[1].set_ylabel('Temperatur / [K]')
                axis[1].set_xlabel('x-R-Ma / [-]')
                axis[1].legend()

                line_3 = LineString(h.np.column_stack((_xin, _tcalcR)))
                line_4 = LineString(h.np.column_stack((_xinRac, _tSLERSR)))
                intersection = line_3.intersection(line_4)

                axis[1].plot(*intersection.xy, 'ro')
                x2, y2 = intersection.xy
                xout1 = x2[0]
                xoutstr1 = "{:10.4f}".format(xout1)
                yout1 = y2[0]
                youtstr1 = "{:10.4f}".format(yout1)
                stringout1 = str(f"SP({xoutstr1}"f"\\{youtstr1})")
                axis[1].text(x2[0], y2[0], 'SP')
                axis[1].text(x2[0], 300, stringout1)
                print(x2, y2)

                plt.show()

                return 0

            def _eq_test():
                _xSexp = [1-0, 1-0.0594, 1-0.0998, 1-0.1496, 1-0.2003, 1-0.2492, 1-0.3098]
                _tSexp = [404.65, 400.65, 399.15, 397.15, 393.65, 391.75, 387.75]

                _xRexp = [1-0, 1-0.1994, 1-0.3192, 1-.6902]

                _xin = h.np.zeros(100)
                _xinRac = h.np.zeros(100)
                _tin = h.np.zeros(100)
                _tcalcS = h.np.zeros(100)
                _tin[0] = 273.15

                for x in range(len(_xin)):
                    _xin[x] = _xin[x - 1] + 0.01
                    _xinRac[x] = _xinRac[x - 1] + 0.005
                    _tin[x] = _tin[0] + (2.0 * x)

                for x in range(len(_xin)):
                    _tcalcS[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0S, _t0S, _gSa, _Alpha))

                figure, axis = plt.subplots(4, constrained_layout=True)

                axis[0].plot(_xSexp, _tSexp, 'ro')
                axis[0].plot(_xin, _tcalcS, '-g')
                axis[0].set_title("S-Fit-Test")
                axis[0].set_ylabel('Temperature')
                axis[0].set_xlabel('xMa')

                axis[1].plot(_xSexp, _tSexp, 'ro')
                axis[1].plot(_xin, _tcalcS, '-g')
                axis[1].set_title("R-Fit-Test")
                axis[1].set_ylabel('Temperature')
                axis[1].set_xlabel('xMa')

                return 0

            def _eut_test():

                _guesses = h.np.array([388.0, 0.6])
                a = h.spo.root(EutFind.calc_eut_ideal, _guesses, args=(_h0S, _t0S, _h0RS, _t0RS), method='hybr')
                b = h.spo.root(EutFind.calc_eut_ideal, _guesses, args=(_h0R, _t0R, _h0RS, _t0RS),)
                print(a.x, 'S-Ma_ideal, Rac-Ma_Prigogine, delta_h-Neumann')
                print(b.x, 'R-Ma_ideal, Rac-Ma_Prigogine, delta_h-Neumann')
                c = h.spo.root(EutFind.calc_eut_nrtl1, _guesses, args=(_h0S, _t0S, _gSa, _h0RS, _t0RS, _Alpha))
                d = h.spo.root(EutFind.calc_eut_nrtl1, _guesses, args=(_h0R, _t0R, _gRa, _h0RS, _t0RS, _Alpha))
                print(c.x, 'S-Ma_NRTL, Rac-Ma_Prigogine, delta_h-Neumann')
                print(d.x, 'R-Ma_NRTL, Rac-Ma_Prigogine, delta_h-Neumann')

                e = h.spo.root(EutFind.calc_eut_ideal, _guesses, args=(_h0S_fit, _t0S, _h0RS_fit, _t0RS), method='hybr')
                f = h.spo.root(EutFind.calc_eut_ideal, _guesses, args=(_h0R_fit, _t0R, _h0RS_fit, _t0RS), )
                print(e.x, 'S-Ma_ideal, Rac-Ma_Prigogine, Thermochimica Acta 382 (2002) 129-142')
                print(f.x, 'R-Ma_ideal, Rac-Ma_Prigogine, Thermochimica Acta 382 (2002) 129-142')
                g = h.spo.root(EutFind.calc_eut_nrtl1, _guesses, args=(_h0S_fit, _t0S, _gSa, _h0RS_fit, _t0RS, _Alpha))
                _h = h.spo.root(EutFind.calc_eut_nrtl1, _guesses, args=(_h0R_fit, _t0R, _gRa, _h0RS_fit, _t0RS, _Alpha))
                print(g.x, 'S-Ma_NRTL, Rac-Ma_Prigogine, Thermochimica Acta 382 (2002) 129-142')
                print(_h.x, 'R-Ma_NRTL, Rac-Ma_Prigogine, Thermochimica Acta 382 (2002) 129-142')

                i = h.spo.root(EutFind.calc_eut_nrtl2, _guesses, args=(_h0S, _t0S, _gSa, _gRa, _h0RS, _t0RS, _Alpha))
                j = h.spo.root(EutFind.calc_eut_nrtl2, _guesses, args=(_h0R, _t0R, _gRa, _gSa, _h0RS, _t0RS, _Alpha))
                print(i.x, 'S-Ma_NRTL, Rac-Ma_Prigogine-NRTL, delta_h-Neumann')
                print(j.x, 'R-Ma_NRTL, Rac-Ma_Prigogine-NRTL, delta_h-Neumann')
                return 0

            def _print():
                _xin = h.np.zeros(100)
                _xinRac = h.np.zeros(100)
                _xinRac[0] = 0.5
                _tin = h.np.zeros(100)
                _tin[0] = 273.15

                _tcalcS = h.np.zeros(100)
                _tcalcRS = h.np.zeros(100)
                _tcalcR = h.np.zeros(100)

                for x in range(len(_xin)):
                    _xin[x] = _xin[x - 1] + 0.01
                    _xinRac[x] = _xinRac[0]+.005*x
                    _tin[x] = _tin[0] + (2.0 * x)

                for x in range(len(_xin)):
                    _tcalcS[x] = h.spo.fsolve(EutFind.t_sle, _tin[x], args=(_xin[x], _h0S, _t0S, _gSa, _Alpha))
                    _tcalcRS[x] = h.spo.fsolve(EutFind.t_sle_rac_nrtl, _tin[x], args=(_xinRac[x], _h0RS, _t0RS, _gSa, _gRa, _Alpha))
                    _tcalcR[x] = h.spo.fsolve(EutFind.t_sle, _tin[x], args=(_xin[x], _h0R, _t0R, _gRa, _Alpha))

                #print(_tcalcS)

                _guesses = h.np.array([388.0, 0.6])
                _EutS = h.spo.root(EutFind.calc_eut_nrtl2, _guesses, args=(_h0S, _t0S, _gSa, _gRa, _h0RS, _t0RS, _Alpha))
                _EutR = h.spo.root(EutFind.calc_eut_nrtl2, _guesses, args=(_h0R, _t0R, _gRa, _gSa, _h0RS, _t0RS, _Alpha))

                _return = (_xinRac, _tcalcRS, _xin, _tcalcS, _tcalcR)
                return _return

            def solve1():
                #_ideal_sle()
                #_nrtl_pure_comp_sle()

                _eut_test()
                _xinRac, _tcalcRS, _xin, _tcalcS, _tcalcR = _print()

                #_EutSDF = h.pd.DataFrame(_EutS, columns=['Param Eut S'])
                #_EutRDF = h.pd.DataFrame(_EutR, columns=['Param Eut R'])
                _xinRacDF = h.pd.DataFrame(_xinRac, columns=['xinRac'])
                _tcalcRSDF = h.pd.DataFrame(_tcalcRS, columns=['t RS'])
                _xinDF = h.pd.DataFrame(_xin, columns=['xin'])
                _tcalcSDF = h.pd.DataFrame(_tcalcS, columns=['t S'])
                _tcalcRDF = h.pd.DataFrame(_tcalcR, columns=['t R'])

                _dataout = [_xinRacDF, _tcalcRSDF, _xinDF, _tcalcSDF, _tcalcRDF]
                _output = h.pd.concat(_dataout, axis=1)
                return _output

            def solve():
                # _ideal_sle()
                #_nrtl_pure_comp_sle()
                #_eut_test()
                _nrtl_nrtl()

                return 0

            return solve()