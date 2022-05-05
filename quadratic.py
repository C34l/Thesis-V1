_a = 1
            _b = -(((np.log(_xi[x])*np.log(_dataoutgammaA[x]))*_t0i)/_h0i)
            _c = -(_t0i/_r)
            _texp = ((_b **2) - (4 * _a * _c))

            if _texp < 0:
                _x1 = 0
                _x2 = 0
            elif _texp == 0:
                _x1 = (-_b + np.sqrt(_texp)) / (2 * _a)
                _dataouttempA.append(_x1)
                print(_x1)
                print("x1 < eindeutig")
            else:
                _x1 = (-_b + np.sqrt(_texp)) / (2 * _a)
                _x2 = (-_b - np.sqrt(_texp)) / (2 * _a)
                if _x1 > 0:
                    _dataouttempA.append(_x1)
                    print(_x1)
                    print("x1  zweideutig")
                else:
                    _dataouttempA.append(_x2)
                    print(_x2)
                    print("x2  zweideutig")