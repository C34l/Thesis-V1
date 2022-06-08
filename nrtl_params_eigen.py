import numpy as np

import headers as h


# replicating values found by Neumann
# using datasets from excel and given constants to loop for activity coefficients via nrtl
# and then looping SLE equation for each gamma calculated to get tcalc
# splits depend on enriched phases, solvent systems and eutectic compositions
_r = 8.31462618
_t0S = 404.75
_t0RS = 393.35
_t0R = 404.65
_h0S = 24500.0
_h0RS = 25600.0
_h0R = 25736.0


class Nrtl:
    @staticmethod
    def nrtl(_xi, _texp, _g):
        _xj = 1.0 - _xi
        _tauij = _g[0] / (_r * _texp)
        _tauji = _g[1] / (_r * _texp)
        _Gij = h.np.exp(-_g[2] * _tauij)
        _Gji = h.np.exp(-_g[3] * _tauji)
        _func = h.np.exp((_xj ** 2) * ((_tauji * (_Gji / (_xi + _xj * _Gji)) ** 2 + _tauij * (_Gij / (_xi * _Gij + _xj) ** 2))))
        return _func

    @staticmethod
    def sle(_texp, _xi, _g, _h0, _t0):
        _yi = Nrtl.nrtl(_xi, _texp, _g)
        # _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gxS)
        _func = 1 / (_yi * _xi) * h.np.exp(-_h0 / (_r * _texp) * (1 - _texp / _t0)) - 1
        return _func

    @staticmethod
    def ma_water(datain):

        _gabS = -2984.75
        _gbaS = 9156.563
        _alpha_abS = 0.686433
        _alpha_baS = 0.217498
        _gS = h.np.array([_gabS, _gbaS, _alpha_abS, _alpha_baS])

        _gabRS = -3352.98
        _gbaRS = 9315.423
        _alpha_abRS = 0.398493
        _alpha_baRS = 0.399798
        _gRS = h.np.array([_gabRS, _gbaRS, _alpha_abRS, _alpha_baRS])

        _gabR = -3999.01
        _gbaR = 10387.36
        _alpha_abR = 0.363806
        _alpha_baR = 0.39205
        _gR = h.np.array([_gabR, _gbaR, _alpha_abR, _alpha_baR])

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
                _gamma = Nrtl.nrtl(_xs[x], _ts[x], _gS)
                _dataoutgammaS.append(_gamma)

                _tcalc = h.spo.fsolve(h.nrpe.Nrtl.sle, _ts[x], args=(_xs[x], _gS, _h0S, _t0S))
                _dataouttempS.append(_tcalc)
            return 0

        def rs_ma():
            _steps = _xrs.count()
            for x in range(_steps):
                _gamma = Nrtl.nrtl(_xrs[x], _trs[x], _gRS)
                _dataoutgammaRS.append(_gamma)

                _tcalc = h.spo.fsolve(h.nrpe.Nrtl.sle, _trs[x], args=(_xrs[x], _gRS, _h0RS, _t0RS))
                _dataouttempRS.append(_tcalc)
            return 0

        def r_ma():
            _steps = _xr.count()
            for x in range(_steps):
                _gamma = Nrtl.nrtl(_xr[x], _tr[x], _gR)
                _dataoutgammaR.append(_gamma)

                _tcalc = h.spo.fsolve(h.nrpe.Nrtl.sle, _tr[x], args=(_xr[x], _gR, _h0R, _t0R))
                _dataouttempR.append(_tcalc)
            return 0

        def solve():
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

        return solve()

    @staticmethod
    def ma_alkyl_lactate(datain):
        _gabSM = -2827.97
        _gbaSM = 9907.41
        _alpha_SM = 0.4
        _gSM = h.np.array([_gabSM, _gbaSM, _alpha_SM])
        _gabRM = -2827.97
        _gbaRM = 9907.41
        _alpha_RM = 0.4
        _gRM = h.np.array([_gabRM, _gbaRM, _alpha_RM])

        _gabSE = -2827.97
        _gbaSE = 9907.41
        _alpha_SE = 0.4
        _gSE = h.np.array([_gabSE, _gbaSE, _alpha_SE])
        _gabRE = -2827.97
        _gbaRE = 9907.41
        _alpha_RE = 0.4
        _gRE = h.np.array([_gabRE, _gbaRE, _alpha_RE])

        _gabSP = -2827.97
        _gbaSP = 9907.41
        _alpha_SP = 0.4
        _gSP = h.np.array([_gabSP, _gbaSP, _alpha_SP])
        _gabRP = -2827.97
        _gbaRP = 9907.41
        _alpha_RP = 0.4
        _gRP = h.np.array([_gabRP, _gbaRP, _alpha_RP])

        _gabSB = -2827.97
        _gbaSB = 9907.41
        _alpha_SB = 0.4
        _gSB = h.np.array([_gabSB, _gbaSB, _alpha_SB])
        _gabRB = -2827.97
        _gbaRB = 9907.41
        _alpha_RB = 0.4
        _gRB = h.np.array([_gabRB, _gbaRB, _alpha_RB])

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

        def methyl():
            _stepsS = _xsM.count
            _stepsR = _xrM.count

            for x in range(_stepsS):
                _gamma_log_S = h.nrpe.Nrtl.nrtl(_xsM[x], _tsM[x], _gSM)
                _gammaS = h.np.exp(_gamma_log_S)
                _dataoutgammaSM.append(_gammaS)
                _tcalcS = h.nrpe.Nrtl.sle(_xsM[x], _tsM[x], _gSM, _h0S, _t0S)
                _dataouttempSM.append(_tcalcS)

            for x in range(_stepsR):
                _gamma_log_R = h.nrpe.Nrtl.nrtl(_xrM[x], _trM[x], _gRM)
                _gammaR = h.np.exp(_gamma_log_R)
                _dataoutgammaSM.append(_gammaR)
                _tcalcR = h.nrpe.Nrtl.sle(_xrM[x], _trM[x], _gRM, _h0R, _t0R)
                _dataouttempSM.append(_tcalcR)
            return 0

        def ethyl():
            _stepsS = _xsE.count
            _stepsR = _xrE.count

            for x in range(_stepsS):
                _gamma_log_S = h.nrpe.Nrtl.nrtl(_xsE[x], _tsE[x], _gSE)
                _gammaS = h.np.exp(_gamma_log_S)
                _dataoutgammaSM.append(_gammaS)
                _tcalcS = h.nrpe.Nrtl.sle(_xsE[x], _tsE[x], _gSE, _h0S, _t0S)
                _dataouttempSM.append(_tcalcS)

            for x in range(_stepsR):
                _gamma_log_R = h.nrpe.Nrtl.nrtl(_xrE[x], _trE[x], _gRE)
                _gammaR = h.np.exp(_gamma_log_R)
                _dataoutgammaSM.append(_gammaR)
                _tcalcR = h.nrpe.Nrtl.sle(_xrE[x], _trE[x], _gRE, _h0R, _t0R)
                _dataouttempSM.append(_tcalcR)
            return 0

        def propyl():
            _stepsS = _xsP.count
            _stepsR = _xrP.count

            for x in range(_stepsS):
                _gamma_log_S = h.nrpe.Nrtl.nrtl(_xsP[x], _tsP[x], _gSP)
                _gammaS = h.np.exp(_gamma_log_S)
                _dataoutgammaSM.append(_gammaS)
                _tcalcS = h.nrpe.Nrtl.sle(_xsP[x], _tsP[x], _gSP, _h0S, _t0S)
                _dataouttempSM.append(_tcalcS)

            for x in range(_stepsR):
                _gamma_log_R = h.nrpe.Nrtl.nrtl(_xrP[x], _trP[x], _gRP)
                _gammaR = h.np.exp(_gamma_log_R)
                _dataoutgammaSM.append(_gammaR)
                _tcalcR = h.nrpe.Nrtl.sle(_xrP[x], _trP[x], _gRP, _h0R, _t0R)
                _dataouttempSM.append(_tcalcR)
            return 0

        def butyl():
            _stepsS = _xsB.count
            _stepsR = _xrB.count

            for x in range(_stepsS):
                _gamma_log_S = h.nrpe.Nrtl.nrtl(_xsB[x], _tsB[x], _gSB)
                _gammaS = h.np.exp(_gamma_log_S)
                _dataoutgammaSM.append(_gammaS)
                _tcalcS = h.nrpe.Nrtl.sle(_xsB[x], _tsB[x], _gSB, _h0S, _t0S)
                _dataouttempSM.append(_tcalcS)

            for x in range(_stepsR):
                _gamma_log_R = h.nrpe.Nrtl.nrtl(_xrB[x], _trB[x], _gRB)
                _gammaR = h.np.exp(_gamma_log_R)
                _dataoutgammaSM.append(_gammaR)
                _tcalcR = h.nrpe.Nrtl.sle(_xrB[x], _trB[x], _gRB, _h0R, _t0R)
                _dataouttempSM.append(_tcalcR)
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