import headers as h

_r= 8.31462618
_t0A = 404.65
_t0B = 404.75
_h0A = 25736.0
_h0B = 24500.0

class NrtlFit:
    @staticmethod
    def estimategamma(datain):
        _alpha = 0.4
        _gab31 = -1403.640
        _gba31 = 2972.166
        _gab69 = -23.298
        _gba69 = 7526.385
        _xi = datain["x - fit"]
        _ti = datain["T - fit"]
        _steps = _xi.count()
        _gammaEstimated = []
        _dataouttempA = []

        for x in range(_steps):
            # print(_xi[x])
            _xa = float(_xi[x])
            _xb = float((1 - _xa))

            if _xi[x] <= 0.31:
                _tauab = (_gab31 / (_r * _ti[x]))
                _tauba = (_gba31 / (_r * _ti[x]))
                _Gab = h.np.exp(-_alpha * _tauab)
                _Gba = h.np.exp(-_alpha * _tauba)
                _gamma = (_xa ** 2) * ((_tauab * (_Gab / (_xb + (_Gab * _xa))) ** 2) + (
                            _tauba * (_Gba / ((_xb * _Gba) + _xa) ** 2)))
                _gammaA = h.np.exp(_gamma)
                _gammaEstimated.append(_gammaA)

                _gxB = _gammaEstimated[x] * _xb
                _htB = _h0B / _t0B
                _tcalc = _h0B / (_htB - _r * h.np.log(_gxB))
                _dataouttempA.append(_tcalc)

                print(_tcalc)
                print(_gammaA)
                print("gamma < o.31")

            elif _xi[x] >= 0.69:
                _tauab = (_gab69 / (_r * _ti[x]))
                _tauba = (_gba69 / (_r * _ti[x]))
                _Gab = h.np.exp((-1 * _alpha) * _tauab)
                _Gba = h.np.exp((-1 * _alpha) * _tauba)
                _gamma = (_xb ** 2) * (
                            (_tauba * (_Gba / (_xa + _Gba * _xb)) ** 2) + (_tauab * (_Gab / (_xa * _Gab + _xb) ** 2)))
                _gammaA = h.np.exp(_gamma)
                _gammaEstimated.append(_gammaA)

                _gxA = _gammaEstimated[x] * _xi[x]
                _htA = _h0A / _t0A
                _tcalc = _h0A / (_htA - _r * h.np.log(_gxA))
                _dataouttempA.append(_tcalc)

                print(_tcalc)
                print(_gammaA)
                print("gamma > o.69")

            else:
                _gammaA = 0
                # _dataoutgammaA.append(_gammaA)

    @staticmethod
    def objectivefunction():
        asd = 3