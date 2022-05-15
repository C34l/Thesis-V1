import headers as h


class NrtlFit:
    @staticmethod
    def datasets(datain):
        _x = datain['x - fit']
        _t = datain['t - fit']
        return _x, _t

    @staticmethod
    def objective_function(_gab, _gba, _xa, _ti):
        _alpha = 0.4
        _xb = (1-_xa)
        _r = 8.31462618

        _tauab = (_gab / (_r * _ti))
        _tauba = (_gba / (_r * _ti))
        _Gab = h.np.exp(-_alpha * _tauab)
        _Gba = h.np.exp(-_alpha * _tauba)

        return (_xa ** 2) * ((_tauab * (_Gab / (_xb + (_Gab * _xa))) ** 2) + (_tauba * (_Gba / ((_xb * _Gba) + _xa) ** 2)))

    @staticmethod
    def popt(_objective_function,_x, _t):
        _popt = h.spo.minimize(NrtlFit.objective_function(_gab, _gba,))