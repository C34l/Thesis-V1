
import headers as h


# replicating values found by Neumann
# using datasets from excel and given constants to loop for activity coefficients via nrtl
# and then looping SLE equation for each gamma calculated to get tcalc
# splits depend on enriched phases, solvent systems and eutectic compositions
class Nrtl:
    @staticmethod
    def binary_pure_ma(datain):
        _alpha = 0.4
        _r = 8.31462618
        _gab31 = -1403.640
        _gba31 = 2972.166
        _gab69 = -23.298
        _gba69 = 7526.385
        _t0A = 404.65
        _t0B = 404.75
        _h0A = 25736.0
        _h0B = 24500.0
        _xi = datain["x - RMA pur"]
        _ti = datain["T - RMA pur"]
        _steps = _xi.count()
        _dataoutgammaA = []
        _dataouttempA = []

        for x in range(_steps):
            # print(_xi[x])
            _xa = float(_xi[x])
            _xb = float((1-_xa))

            if _xi[x] <= 0.31:
                _tauab = (_gab31 / (_r * _ti[x]))
                _tauba = (_gba31 / (_r * _ti[x]))
                _Gab = h.np.exp(-_alpha * _tauab)
                _Gba = h.np.exp(-_alpha * _tauba)
                _gamma = (_xa**2)*((_tauab*(_Gab/(_xb+(_Gab*_xa)))**2)+(_tauba * (_Gba/((_xb*_Gba)+_xa)**2)))
                _gammaA = h.np.exp(_gamma)
                _dataoutgammaA.append(_gammaA)

                _gxB = _dataoutgammaA[x] * _xb
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
                _gamma = (_xb ** 2) * ((_tauba * (_Gba / (_xa + _Gba * _xb)) ** 2) + (_tauab * (_Gab / (_xa * _Gab + _xb) ** 2)))
                _gammaA = h.np.exp(_gamma)
                _dataoutgammaA.append(_gammaA)

                _gxA = _dataoutgammaA[x] * _xi[x]
                _htA = _h0A / _t0A
                _tcalc = _h0A / (_htA - _r * h.np.log(_gxA))
                _dataouttempA.append(_tcalc)

                print(_tcalc)
                print(_gammaA)
                print("gamma > o.69")

            else:
                _gammaA = 0
                # _dataoutgammaA.append(_gammaA)

        # for x in range(_steps):
            # _h0i = (_h0A * _xi[x])+(_h0B * (1-_xi[x]))
            # _t0i = (_t0A * _xi[x])+(_t0B * (1-_xi[x]))
            # _gx = _dataoutgammaA[x]*_xi[x]
            # _ht0 = _h0i/_t0i
            # _tcalc = -_h0i/(_r*np.log(_gx)-_ht0)
            # _dataouttempA.append(_tcalc)
            # print(_tcalc)

        # print(_gamma)
        _dataoutgammaDF = h.pd.DataFrame(_dataoutgammaA, columns=['gamma'])
        _dataouttempDF = h.pd.DataFrame(_dataouttempA, columns=['Temp'])

        _dataout = [_dataoutgammaDF, _dataouttempDF, datain["x - RMA pur"]]
        _result = h.pd.concat(_dataout, axis=1)
        # print(_result)

        return _result

    @staticmethod
    def ma_water(datain):
        _alpha = 0.4
        _r = 8.31462618
        _gabS = -2827.97
        _gbaS = 9907.41
        _gabRS = -3143.16
        _gbaRS = 9443.10
        _gabR = -3509.61
        _gbaR = 9841.69
        _t0S = 404.75
        _t0RS = 393.35
        _t0R = 404.65
        _h0S = 24500.0
        _h0RS = 25600.0
        _h0R = 25736.0

        _dataoutgammaS = []
        _dataoutgammaRS = []
        _dataoutgammaR = []
        _dataouttempS = []
        _dataouttempRS = []
        _dataouttempR = []

        _xs = datain["x - SWasser"]
        _ts = datain["T - SWasser"]
        _xrs = datain["x - RSWasser"]
        _trs = datain["T - RSWasser"]
        _xr = datain["x - RWasser"]
        _tr = datain["T - RWasser"]
        def s_ma():
            _steps = _xs.count()
            for x in range(_steps):
                _xMA = _xs[x]
                _xSolvent = float((1 - _xs[x]))
                _tauabS = (_gabS / (_r * _ts[x]))
                _taubaS = (_gbaS / (_r * _ts[x]))
                _GabS = h.np.exp(-_alpha * _tauabS)
                _GbaS = h.np.exp(-_alpha * _taubaS)
                _gammaS = (_xSolvent**2)*((_taubaS*(_GbaS/(_xMA+(_GbaS*_xSolvent)))**2)+(_tauabS * (_GabS/((_xMA*_GabS)+_xSolvent)**2)))
                _gammaAS = h.np.exp(_gammaS)
                _dataoutgammaS.append(_gammaAS)
                _gxS = _dataoutgammaS[x] * _xMA
                _htS = _h0S / _t0S
                _tcalcS = _h0S / (_htS - _r * h.np.log(_gxS))
                _dataouttempS.append(_tcalcS)
                # print(_tcalcS)
                # print(_gammaS)
                # print("gamma S")
            return 0

        def rs_ma():
            _steps = _xrs.count()
            for x in range(_steps):
                _xMA = _xrs[x]
                _xSolvent = float((1 - _xrs[x]))
                _tauabRS = (_gabRS / (_r * _trs[x]))
                _taubaRS = (_gbaRS / (_r * _trs[x]))
                _GabRS = h.np.exp(-_alpha * _tauabRS)
                _GbaRS = h.np.exp(-_alpha * _taubaRS)
                _gammaRS = (_xSolvent ** 2) * ((_taubaRS * (_GbaRS / (_xMA + (_GbaRS * _xSolvent))) ** 2) + (
                            _tauabRS * (_GabRS / ((_xMA * _GabRS) + _xSolvent) ** 2)))
                _gammaARS = h.np.exp(_gammaRS)
                _dataoutgammaRS.append(_gammaARS)
                _gxRS = _dataoutgammaRS[x] * _xMA
                _htRS = _h0RS / _t0RS
                _tcalcRS = _h0RS / (_htRS - _r * h.np.log(_gxRS))
                _dataouttempRS.append(_tcalcRS)
                # print(_tcalcRS)
                # print(_gammaRS)
                # print("gamma RS")
            return 0

        def r_ma():
            _steps = _xr.count()
            for x in range(_steps):
                _xMA = _xr[x]
                _xSolvent = float((1 - _xr[x]))
                _tauabR = (_gabR / (_r * _tr[x]))
                _taubaR = (_gbaR / (_r * _tr[x]))
                _GabR = h.np.exp(-_alpha * _tauabR)
                _GbaR = h.np.exp(-_alpha * _taubaR)
                _gammaR = (_xSolvent ** 2) * ((_taubaR * (_GbaR / (_xMA + (_GbaR * _xSolvent))) ** 2) + (
                            _tauabR * (_GabR / ((_xMA * _GabR) + _xSolvent) ** 2)))
                _gammaAR = h.np.exp(_gammaR)
                _dataoutgammaR.append(_gammaAR)
                _gxR = _dataoutgammaR[x] * _xMA
                _htR = _h0R / _t0R
                _tcalcR = _h0R / (_htR - _r * h.np.log(_gxR))
                _dataouttempR.append(_tcalcR)
                # print(_tcalcR)
                # print(_gammaAR)
                # print("gamma R")
            return 0


        s_ma()
        rs_ma()
        r_ma()

        _dataoutgammaDFS = h.pd.DataFrame(_dataoutgammaS, columns=['gamma S'])
        _dataoutgammaDFRS = h.pd.DataFrame(_dataoutgammaRS, columns=['gamma RS'])
        _dataoutgammaDFR = h.pd.DataFrame(_dataoutgammaR, columns=['gamma R'])
        _dataouttempDFS = h.pd.DataFrame(_dataouttempS, columns=['Temp S'])
        _dataouttempDFRS = h.pd.DataFrame(_dataouttempRS, columns=['Temp RS'])
        _dataouttempDFR = h.pd.DataFrame(_dataouttempR, columns=['Temp R'])

        _dataout = [datain["x - SWasser"], _dataoutgammaDFS, _dataouttempDFS, datain["x - RSWasser"], _dataoutgammaDFRS, _dataouttempDFRS, datain["x - RWasser"], _dataoutgammaDFR, _dataouttempDFR]
        _result = h.pd.concat(_dataout, axis=1)

        return _result

    @staticmethod
    def ma_lactate(datain):
        _alpha = 0.4
        _r = 8.31462618
        _gabS = -2827.97
        _gbaS = 9907.41
        _gabRS = -3143.16
        _gbaRS = 9443.10
        _gabR = -3509.61
        _gbaR = 9841.69
        _t0S = 404.75
        _t0RS = 393.35
        _t0R = 404.65
        _h0S = 24500.0
        _h0RS = 25600.0
        _h0R = 25736.0

        _dataoutgammaS = []
        _dataoutgammaRS = []
        _dataoutgammaR = []
        _dataouttempS = []
        _dataouttempRS = []
        _dataouttempR = []

        _xs = datain["x - SLak"]
        _ts = datain["T - SLak"]
        _xrs = datain["x - RSLak"]
        _trs = datain["T - RSLak"]
        _xr = datain["x - RLak"]
        _tr = datain["T - RLak"]

        def s_ma():
            _steps = _xs.count()
            for x in range(_steps):
                _xMA = _xs[x]
                _xSolvent = float((1 - _xs[x]))
                _tauabS = (_gabS / (_r * _ts[x]))
                _taubaS = (_gbaS / (_r * _ts[x]))
                _GabS = h.np.exp(-_alpha * _tauabS)
                _GbaS = h.np.exp(-_alpha * _taubaS)
                _gammaS = (_xSolvent ** 2) * ((_taubaS * (_GbaS / (_xMA + (_GbaS * _xSolvent))) ** 2) + (
                            _tauabS * (_GabS / ((_xMA * _GabS) + _xSolvent) ** 2)))
                _gammaAS = h.np.exp(_gammaS)
                _dataoutgammaS.append(_gammaAS)
                _gxS = _dataoutgammaS[x] * _xMA
                _htS = _h0S / _t0S
                _tcalcS = _h0S / (_htS - _r * h.np.log(_gxS))
                _dataouttempS.append(_tcalcS)
            # print(_tcalcS)
            # print(_gammaS)
            # print("gamma S")
            return 0

        def rs_ma():
            _steps = _xrs.count()
            for x in range(_steps):
                _xMA = _xrs[x]
                _xSolvent = float((1 - _xrs[x]))
                _tauabRS = (_gabRS / (_r * _trs[x]))
                _taubaRS = (_gbaRS / (_r * _trs[x]))
                _GabRS = h.np.exp(-_alpha * _tauabRS)
                _GbaRS = h.np.exp(-_alpha * _taubaRS)
                _gammaRS = (_xSolvent ** 2) * ((_taubaRS * (_GbaRS / (_xMA + (_GbaRS * _xSolvent))) ** 2) + (
                        _tauabRS * (_GabRS / ((_xMA * _GabRS) + _xSolvent) ** 2)))
                _gammaARS = h.np.exp(_gammaRS)
                _dataoutgammaRS.append(_gammaARS)
                _gxRS = _dataoutgammaRS[x] * _xMA
                _htRS = _h0RS / _t0RS
                _tcalcRS = _h0RS / (_htRS - _r * h.np.log(_gxRS))
                _dataouttempRS.append(_tcalcRS)
            # print(_tcalcRS)
            # print(_gammaRS)
            # print("gamma RS")
            return 0

        def r_ma():
            _steps = _xr.count()
            for x in range(_steps):
                _xMA = _xr[x]
                _xSolvent = float((1 - _xr[x]))
                _tauabR = (_gabR / (_r * _tr[x]))
                _taubaR = (_gbaR / (_r * _tr[x]))
                _GabR = h.np.exp(-_alpha * _tauabR)
                _GbaR = h.np.exp(-_alpha * _taubaR)
                _gammaR = (_xSolvent ** 2) * ((_taubaR * (_GbaR / (_xMA + (_GbaR * _xSolvent))) ** 2) + (
                        _tauabR * (_GabR / ((_xMA * _GabR) + _xSolvent) ** 2)))
                _gammaAR = h.np.exp(_gammaR)
                _dataoutgammaR.append(_gammaAR)
                _gxR = _dataoutgammaR[x] * _xMA
                _htR = _h0R / _t0R
                _tcalcR = _h0R / (_htR - _r * h.np.log(_gxR))
                _dataouttempR.append(_tcalcR)
            # print(_tcalcR)
            # print(_gammaAR)
            # print("gamma R")
            return 0

        s_ma()
        rs_ma()
        r_ma()

        _dataoutgammaDFS = h.pd.DataFrame(_dataoutgammaS, columns=['gamma S'])
        _dataoutgammaDFRS = h.pd.DataFrame(_dataoutgammaRS, columns=['gamma RS'])
        _dataoutgammaDFR = h.pd.DataFrame(_dataoutgammaR, columns=['gamma R'])
        _dataouttempDFS = h.pd.DataFrame(_dataouttempS, columns=['Temp S'])
        _dataouttempDFRS = h.pd.DataFrame(_dataouttempRS, columns=['Temp RS'])
        _dataouttempDFR = h.pd.DataFrame(_dataouttempR, columns=['Temp R'])

        _dataout = [datain["x - SLak"], _dataoutgammaDFS, _dataouttempDFS, datain["x - RSLak"], _dataoutgammaDFRS,
                    _dataouttempDFRS, datain["x - RLak"], _dataoutgammaDFR, _dataouttempDFR]
        _result = h.pd.concat(_dataout, axis=1)

        return _result