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

        for x in np.nditer(_xi):
            print(_xi)

        return 0

