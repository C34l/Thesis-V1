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
    def ma_alkyl_lactate(datain):
        _alpha = 0.4
        _r = 8.31462618

        _gabSM = -5145.15
        _gbaSM = 4077.76
        _gabRM = -5031.57
        _gbaRM = 4787.22

        _gabSE = -4883.75
        _gbaSE = 12111.15
        _gabRE = -4833.93
        _gbaRE = 10744.97

        _gabSP = -5013.79
        _gbaSP = 5442.11
        _gabRP = -4696.62
        _gbaRP = 4967.03

        _gabSB = -4621.74
        _gbaSB = 4763.98
        _gabRB = -4345.53
        _gbaRB = 4545.06

        _t0S = 404.75
        _t0R = 404.65
        _h0S = 24500.0
        _h0R = 25736.0

        _dataoutgammaSM = []
        _dataoutgammaRM = []
        _dataouttempSM = []
        _dataouttempRM = []

        _dataoutgammaSE = []
        _dataoutgammaRE = []
        _dataouttempSE = []
        _dataouttempRE = []

        _dataoutgammaSP = []
        _dataoutgammaRP = []
        _dataouttempSP = []
        _dataouttempRP = []

        _dataoutgammaSB = []
        _dataoutgammaRB = []
        _dataouttempSB = []
        _dataouttempRB = []

        _xsM = datain["x - SLakM"]
        _tsM = datain["T - SLakM"]
        _xrM = datain["x - RLakM"]
        _trM = datain["T - RLakM"]

        _xsE = datain["x - SLakE"]
        _tsE = datain["T - SLakE"]
        _xrE = datain["x - RLakE"]
        _trE = datain["T - RLakE"]

        _xsP = datain["x - SLakP"]
        _tsP = datain["T - SLakP"]
        _xrP = datain["x - RLakP"]
        _trP = datain["T - RLakP"]

        _xsB = datain["x - SLakB"]
        _tsB = datain["T - SLakB"]
        _xrB = datain["x - RLakB"]
        _trB = datain["T - RLakB"]

        def nrtl(_xi, _gij, _gji, _texp, _a):
            _xj = float(1.0 - _xi)
            _tauij = _gij / (_r * _texp)
            _tauji = _gji / (_r * _texp)
            _Gij = h.np.exp(-_a * _tauij)
            _Gji = h.np.exp(-_a * _tauji)
            _func = ((_xj ** 2) * ((_tauji * (_Gji / (_xi + _xj * _Gji)) ** 2 + _tauij * (_Gij / (_xi * _Gij + _xj) ** 2))))
            return _func

        def sle(_yi, _xi, _h0, _t0):
            _gxS = _yi * _xi
            _ht = _h0 / _t0
            _func = _h0 / (_ht - _r * h.np.log(_gxS))
            return _func

        def methyl():
            _stepsS = _xsM.count()
            _stepsR = _xrM.count()

            for x in range(_stepsS):
                _gammaS = nrtl(_xsM[x], _gabSM, _gbaSM, _tsM[x], _alpha)
                _gammaAS = h.np.exp(_gammaS)
                _dataoutgammaSM.append(_gammaAS)
                _tcalcS = sle(_gammaAS, _xsM[x], _h0S, _t0S)
                _dataouttempSM.append(_tcalcS)

            for x in range(_stepsR):
                _gammaR = nrtl(_xrM[x], _gabRM, _gbaRM, _trM[x], _alpha)
                _gammaAR = h.np.exp(_gammaR)
                _dataoutgammaRM.append(_gammaAR)
                _tcalcR = sle(_gammaAR, _xrM[x], _h0R, _t0R)
                _dataouttempRM.append(_tcalcR)
            return 0

        def ethyl():
            _stepsS = _xsE.count()
            _stepsR = _xrE.count()

            for x in range(_stepsS):
                _gammaS = nrtl(_xsE[x], _gabSE, _gbaSE, _tsE[x], _alpha)
                _gammaAS = h.np.exp(_gammaS)
                _dataoutgammaSE.append(_gammaAS)
                _tcalcS = sle(_gammaAS, _xsE[x], _h0S, _t0S)
                _dataouttempSE.append(_tcalcS)

            for x in range(_stepsR):
                _gammaR = nrtl(_xrE[x], _gabRE, _gbaRE, _trE[x], _alpha)
                _gammaAR = h.np.exp(_gammaR)
                _dataoutgammaRE.append(_gammaAR)
                _tcalcR = sle(_gammaAR, _xrE[x], _h0R, _t0R)
                _dataouttempRE.append(_tcalcR)
            return 0

        def propyl():
            _stepsS = _xsP.count()
            _stepsR = _xrP.count()

            for x in range(_stepsS):
                _gammaS = nrtl(_xsP[x], _gabSP, _gbaSP, _tsP[x], _alpha)
                _gammaAS = h.np.exp(_gammaS)
                _dataoutgammaSP.append(_gammaAS)
                _tcalcS = sle(_gammaAS, _xsP[x], _h0S, _t0S)
                _dataouttempSP.append(_tcalcS)

            for x in range(_stepsR):
                _gammaR = nrtl(_xrP[x], _gabRP, _gbaRP, _trP[x], _alpha)
                _gammaAR = h.np.exp(_gammaR)
                _dataoutgammaRP.append(_gammaAR)
                _tcalcR = sle(_gammaAR, _xrP[x], _h0R, _t0R)
                _dataouttempRP.append(_tcalcR)
            return 0

        def butyl():
            _stepsS = _xsB.count()
            _stepsR = _xrB.count()

            for x in range(_stepsS):
                _gammaS = nrtl(_xsB[x], _gabSB, _gbaSB, _tsB[x], _alpha)
                _gammaAS = h.np.exp(_gammaS)
                _dataoutgammaSB.append(_gammaAS)
                _tcalcS = sle(_gammaAS, _xsB[x], _h0S, _t0S)
                _dataouttempSB.append(_tcalcS)

            for x in range(_stepsR):
                _gammaR = nrtl(_xrB[x], _gabRB, _gbaRB, _trB[x], _alpha)
                _gammaAR = h.np.exp(_gammaR)
                _dataoutgammaRB.append(_gammaAR)
                _tcalcR = sle(_gammaAR, _xrB[x], _h0R, _t0R)
                _dataouttempRB.append(_tcalcR)
            return 0

        def solve():
            methyl()
            ethyl()
            propyl()
            butyl()

            _dataoutgammaDFSM = h.pd.DataFrame(_dataoutgammaSM, columns=['gamma SM'])
            _dataoutgammaDFRM = h.pd.DataFrame(_dataoutgammaRM, columns=['gamma RM'])
            _dataouttempDFSM = h.pd.DataFrame(_dataouttempSM, columns=['Temp SM'])
            _dataouttempDFRM = h.pd.DataFrame(_dataouttempRM, columns=['Temp RM'])

            _dataoutgammaDFSE = h.pd.DataFrame(_dataoutgammaSE, columns=['gamma SE'])
            _dataoutgammaDFRE = h.pd.DataFrame(_dataoutgammaRE, columns=['gamma RE'])
            _dataouttempDFSE = h.pd.DataFrame(_dataouttempSE, columns=['Temp SE'])
            _dataouttempDFRE = h.pd.DataFrame(_dataouttempRE, columns=['Temp RE'])

            _dataoutgammaDFSP = h.pd.DataFrame(_dataoutgammaSP, columns=['gamma SP'])
            _dataoutgammaDFRP = h.pd.DataFrame(_dataoutgammaRP, columns=['gamma RP'])
            _dataouttempDFSP = h.pd.DataFrame(_dataouttempSP, columns=['Temp SP'])
            _dataouttempDFRP = h.pd.DataFrame(_dataouttempRP, columns=['Temp RP'])

            _dataoutgammaDFSB = h.pd.DataFrame(_dataoutgammaSB, columns=['gamma SB'])
            _dataoutgammaDFRB = h.pd.DataFrame(_dataoutgammaRB, columns=['gamma RB'])
            _dataouttempDFSB = h.pd.DataFrame(_dataouttempSB, columns=['Temp SB'])
            _dataouttempDFRB = h.pd.DataFrame(_dataouttempRB, columns=['Temp RB'])

            _dataout = [datain["x - SLakM"], _dataoutgammaDFSM, _dataouttempDFSM,
                        datain["x - RLakM"], _dataoutgammaDFRM, _dataouttempDFRM,
                        datain["x - SLakE"], _dataoutgammaDFSE, _dataouttempDFSE,
                        datain["x - RLakE"], _dataoutgammaDFRE, _dataouttempDFRE,
                        datain["x - SLakP"], _dataoutgammaDFSP, _dataouttempDFSP,
                        datain["x - RLakP"], _dataoutgammaDFRP, _dataouttempDFRP,
                        datain["x - SLakB"], _dataoutgammaDFSB, _dataouttempDFSB,
                        datain["x - RLakB"], _dataoutgammaDFRB, _dataouttempDFRB,
                        ]
            _result = h.pd.concat(_dataout, axis=1)

            return _result

        return solve()
