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
    def t_sle_ideal(_texp, _xi, _h0, _t0, _g, _alpha):
        _yi = 1
        _gx = _yi * (_xi)
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx) - 1
        return _func

    @staticmethod
    def pure_ma():

            _gabS = -2827.97
            _gbaS = 9907.41
            _gSa = h.np.array([_gabS, _gbaS, ])

            _gabRS = -3143.16
            _gbaRS = 9443.10
            _gRSa = h.np.array([_gabRS, _gbaRS, ])

            _gabR = -3509.61
            _gbaR = 9841.69
            _gRa = h.np.array([_gabR, _gbaR, ])

            def _test():
                _xin = h.np.zeros(100)
                _xinRac = h.np.zeros(100)
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
                    _xinRac[x] = _xinRac[x-1]+0.005

                for x in range(len(_tin)):
                    _tin[x] = _tin[0]+(2.0*x)

                # loading up calc points
                for x in range(len(_xin)):
                    _tcalcS[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0S, _t0S, _gSa, _Alpha))
                    _tcalcRSS[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0RS, _t0RS, _gRSa, _Alpha))
                    _tcalcRSR[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0RS, _t0RS, _gRSa, _Alpha))
                    _tcalcR[x] = h.spo.fsolve(EutFind.t_sle_ideal, _tin[x], args=(_xin[x], _h0R, _t0R, _gRa, _Alpha))



                xpS = _xin
                ypS = _tcalcS
                plt.plot(xpS, ypS, '-g')

                xrsS = -_xinRac+.5
                yrsS = _tcalcRSS
                plt.plot(xrsS, yrsS, '--g')

                line_1 = LineString(h.np.column_stack((xpS, ypS)))
                line_2 = LineString(h.np.column_stack((xrsS, yrsS)))
                intersection = line_1.intersection(line_2)

                plt.plot(*intersection.xy, 'ro')


                x, y = intersection.xy
                print(x, y)

                xpR = _xin
                ypR = _tcalcR
                plt.plot(xpS, ypR, '-b')

                xrsR = -_xinRac + .5
                yrsR = _tcalcRSS
                plt.plot(xrsR, yrsR, '--b')

                line_3 = LineString(h.np.column_stack((xpR, ypR)))
                line_4 = LineString(h.np.column_stack((xrsR, yrsR)))
                intersection = line_3.intersection(line_4)

                plt.plot(*intersection.xy, 'ro')
                x2, y2 = intersection.xy
                plt.show()
                print(x2, y2)

                return 0


            def solve():
                _test()
                _output = 3
                return _output

            return solve()