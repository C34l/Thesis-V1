import numpy as np
import headers as h


class Nrtl:
    @staticmethod
    def approximategamma(datain):
        _alpha = 0.4
        _r = 8.3147
        _gab31 = -1403.640
        _gba31 = 2972.166
        _gab69 = -23.298
        _gba69 = 7526.385
        _t0A = 404.65
        _t0B = 404.75
        _h0A = 25736.26
        _h0B = 24500.0
        _xi = datain["x+"]
        _ti = datain["T+"]
        _steps = len(datain)
        _dataoutgammaA = []
        _dataouttempA = []


        for x in range(_steps):
            #print(_xi[x])
            _xa = _xi[x]
            _xb = (1-_xa)
            _t = _ti[x]

            if _xa <= 0.32:
                _tauab = (_gab31 / (_r * _t))
                _tauba = (_gba31 / (_r * _t))
                _Gab = np.exp((-1* _alpha) *_tauab)
                _Gba = np.exp((-1* _alpha) *_tauba)
                #_gamma = (_xb**2)*((_tauba*(_Gba/(_xa+(_Gba*_xb)))**2)+(_tauab+_Gab/((_xa*_Gab)+_xb)**2))
                _gamma = (_xa**2)*((_tauab*(_Gab/(_xb+(_Gab*_xa)))**2)+(_tauba+_Gba/((_xb*_Gba)+_xa)**2))
                _gammaA = np.exp(_gamma)
                _dataoutgammaA.append(_gammaA)

                _h0i = (_h0A * _xi[x]) + (_h0B * (1 - _xi[x]))
                _t0i = (_t0A * _xi[x]) + (_t0B * (1 - _xi[x]))
                _gx = _dataoutgammaA[x] * (1-_xi[x])
                _ht0 = _h0i / _t0i
                _tcalc = -_h0i / (_r * np.log(_gx) - _ht0)
                _dataouttempA.append(_tcalc)

                print(_tcalc)
                print(_gammaA)
                print("gamma < o.31")


            elif _xa >= 0.69:
                _tauab = (_gab69 / (_r * _t))
                _tauba = (_gba69 / (_r * _t))
                _Gab = np.exp((-1 * _alpha) * _tauab)
                _Gba = np.exp((-1 * _alpha) * _tauba)
                #_gamma = (_xb * _xb) * (_tauba * ((_Gba / (_xa + (_Gba + _xb))) * (_Gba / (_xa + (_Gba + _xb)))) + _tauab * ((_Gab / ((_xa * _Gab) * _xb)) * (_Gab / ((_xa * _Gab) * _xb))))
                _gamma = (_xb ** 2) * ((_tauba * (_Gba / (_xa + _Gba * _xb)) ** 2) + (_tauab + _Gab / (_xa * _Gab + _xb) ** 2))
                _gammaA = np.exp(_gamma)
                _dataoutgammaA.append(_gammaA)

                _h0i = (_h0A * _xi[x]) + (_h0B * (1 - _xi[x]))
                _t0i = (_t0A * _xi[x]) + (_t0B * (1 - _xi[x]))
                _gx = _dataoutgammaA[x] * _xi[x]
                _ht0 = _h0i / _t0i
                _tcalc = -_h0i / (_r * np.log(_gx) - _ht0)
                _dataouttempA.append(_tcalc)

                print(_tcalc)
                print(_gammaA)
                print("gamma > o.69")

            else:
                _gammaA=0
                #_dataoutgammaA.append(_gammaA)

        #for x in range(_steps):
            #_h0i = (_h0A * _xi[x])+(_h0B * (1-_xi[x]))
            #_t0i = (_t0A * _xi[x])+(_t0B * (1-_xi[x]))
            #_gx = _dataoutgammaA[x]*_xi[x]
            #_ht0 = _h0i/_t0i
            #_tcalc = -_h0i/(_r*np.log(_gx)-_ht0)
            #_dataouttempA.append(_tcalc)
            #print(_tcalc)

        #print(_gamma)
        _dataoutgammaDF = h.pd.DataFrame(_dataoutgammaA, columns=['gamma'])
        _dataouttempDF = h.pd.DataFrame(_dataouttempA, columns=['Temp'])

        _dataout = [_dataoutgammaDF, _dataouttempDF, datain["x+"]]
        _result = h.pd.concat(_dataout, axis=1)
        #print(_result)
        return _result

