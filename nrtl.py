import numpy as np
import headers as h


class Nrtl:
    @staticmethod
    def approximate(data):
        _alpha = 0.4
        _r = 8.3147
        _gab = 1
        _gba = 2
        _xi = data["x+"]
        _ti = data ["T+"]
        steps = len(data)
        for x in range(steps):
            print(_xi[x])
            _xa = _xi[x]
            _xb = (1-_xa)
            _t = _ti[x]
            _tauab = (_gab/(_r*_t))
            _tauba = (_gba/(_r*_t))
            if _xa <= 0.31:
                gamma = np.square(_xb)*

        return 0

