
import scipy.optimize

import headers as h


_r = 8.31462618
_t0S = 404.75
_t0RS = 393.35
_t0R = 404.65
_h0S = 24500.0
_h0RS = 25600.0
_h0R = 25736.0


class NrtlFit:
    @staticmethod
    def y_nrtl(_g, _xi, _texp,):
        _xj = float(1.0 - _xi)
        _tauij = _g[0] / (_r * _texp)
        _tauji = _g[1] / (_r * _texp)
        _Gij = h.np.exp(-_g[2] * _tauij)
        _Gji = h.np.exp(-_g[3] * _tauji)
        _func = ((_xj ** 2) * (_tauji * (_Gji / (_xi + _xj * _Gji)) ** 2 + _tauij * (_Gij / (_xi * _Gij + _xj) ** 2)))
        return _func

    @staticmethod
    def t_sle(_texp, _xi, _h0, _t0, _g,):
        _yi = NrtlFit.y_nrtl(_g, _xi, _texp)
        _gxS = h.np.exp(_yi) * _xi
        _ht = _h0 / _t0
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gxS)
        return _func

    @staticmethod
    def min_fqs(_g, _xi, _texp, _h0, _t0, steps):
        _tmod = h.np.zeros(steps)
        for x in range(steps):
            _tmod0 = _texp[x]
            _tmod1 = scipy.optimize.fsolve(NrtlFit.t_sle, _tmod0, args=(_xi[x], _h0, _t0, _g))
            _tmod[x] = (_tmod1)

        # fqs_terme = (_texp[0:steps] - _tmod) ** 2
        fqs_norm = ((_texp[0:steps] - _tmod[0:steps]) / _texp[0:steps]) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    @staticmethod
    def parametersH2O(datain):

            _alphaabS = 0.4
            _alphabaS = 0.4
            _gabS = -2827.97
            _gbaS = 9907.41
            _gSa = h.np.array([_gabS, _gbaS, _alphaabS, _alphabaS])
            _alphaabRS = 0.4
            _alphabaRS = 0.4
            _gabRS = -3143.16
            _gbaRS = 9443.10
            _gRSa = h.np.array([_gabRS, _gbaRS, _alphaabRS, _alphabaRS])
            _alphaabR = 0.4
            _alphabaR = 0.4
            _gabR = -3509.61
            _gbaR = 9841.69
            _gRa = h.np.array([_gabR, _gbaR, _alphaabR, _alphabaR])

            _xs = datain["x - SWasser"]
            _tSexp = datain["T - SWasser"]
            _xrs = datain["x - RSWasser"]
            _tRSexp = datain["T - RSWasser"]
            _xr = datain["x - RWasser"]
            _tRexp = datain["T - RWasser"]

            _bnds1 = ((-3000.0, -2000.0), (9000.0, 10000.0), (0.0, 1.0), (0.0, 1.0))
            _bnds2 = ((-4000.0, -3000.0), (9000.0, 10000.0), (0.0, 1.0), (0.0, 1.0))
            _bnds3 = ((-4000.0, -3000.0), (9000.0, 12000.0), (0.0, 1.0), (0.0, 1.0))

            def s_ma():
                steps_x = _xs.count()
                steps_t = _tSexp.count()

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gSa, args=(_xs, _tSexp, _h0S, _t0S, steps_x), method='Nelder-Mead', bounds=_bnds1)
                print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_1.fun, _tSexp)))))
                print('deltag_AB s = ' + str(res_1.x[0]))
                print('deltag_BA s = ' + str(res_1.x[1]))
                print('alpha ab s =', res_1.x[2])
                print('alpha ba s =', res_1.x[3])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)

                _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])

                for x in range(steps_t):
                    _tcalc[x] = NrtlFit.t_sle(_tSexp[x], _xs[x], _h0S, _t0S, _gneu)

                for x in range(steps_t):
                    t_diff_norm[x] = abs((_tSexp[x] - _tcalc[x]) / _tSexp[x])
                    t_diff[x] = abs((_tSexp[x] - _tcalc[x]))
                ard_neu_norm = 100 / len(_tSexp) * sum(t_diff_norm)
                ard_neu = 100 / len(_tSexp) * sum(t_diff)
                print('ARD_neu_norm_bzgl.T [%] =', ard_neu_norm)
                print('ARD_neu_bzgl.T [%] =', ard_neu)
                _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                return _params

            def rs_ma():
                steps_x = _xrs.count()
                steps_t = _tRSexp.count()

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gRSa, args=(_xrs, _tRSexp, _h0RS, _t0RS, steps_x),
                                                method='Nelder-Mead', bounds=_bnds2)
                print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_1.fun, _tRSexp)))))
                print('deltag_AB [ ] = ' + str(res_1.x[0]))
                print('deltag_BA [ ] = ' + str(res_1.x[1]))
                print('alpha =', res_1.x[2])
                print('alpha ba s =', res_1.x[3])

                t_diff_normRS = h.np.zeros(steps_t)
                t_diffRS = h.np.zeros(steps_t)
                _tcalcRS = h.np.zeros(steps_t)

                _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])

                for x in range(steps_t):
                    _tcalcRS[x] = NrtlFit.t_sle(_tRSexp[x], _xrs[x], _h0RS, _t0RS, _gneu)

                for x in range(steps_t):
                    t_diff_normRS[x] = abs((_tRSexp[x] - _tcalcRS[x]) / _tRSexp[x])
                    t_diffRS[x] = abs((_tRSexp[x] - _tcalcRS[x]))
                ard_neu_norm_rs = 100 / len(_tRSexp) * sum(t_diff_normRS)
                ard_neu_rs = 100 / len(_tRSexp) * sum(t_diffRS)
                print('ARD_neu_norm_bzgl.T [%] =', ard_neu_norm_rs)
                print('ARD_neu_bzgl.T [%] =', ard_neu_rs)
                _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                return _params

            def r_ma():
                steps_x = _xr.count()
                steps_t = _tRexp.count()

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gRa, args=(_xr, _tRexp, _h0R, _t0R, steps_x),
                                                method='Nelder-Mead', bounds=_bnds3)
                print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_1.fun, _tRexp)))))
                print('deltag_AB [ ] = ' + str(res_1.x[0]))
                print('deltag_BA [ ] = ' + str(res_1.x[1]))
                print('alpha =', res_1.x[2])
                print('alpha ba s =', res_1.x[3])

                _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                for x in range(steps_t):
                    _tcalc[x] = NrtlFit.t_sle(_tRexp[x], _xr[x], _h0R, _t0R, _gneu)

                for x in range(steps_t):
                    t_diff_norm[x] = abs((_tRexp[x] - _tcalc[x]) / _tRexp[x])
                    t_diff[x] = abs((_tRexp[x] - _tcalc[x]))
                ard_neu_norm = 100 / len(_tRexp) * sum(t_diff_norm)
                ard_neu = 100 / len(_tRexp) * sum(t_diff)
                print('ARD_neu_norm_bzgl.T [%] =', ard_neu_norm)
                print('ARD_neu_bzgl.T [%] =', ard_neu)
                _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                return _params

            def solve():
                _dataoutParamsS = s_ma()
                _dataoutParamsRS = rs_ma()
                _dataoutParamsR = r_ma()

                _dataoutFrame = ('gAB', 'gBA', 'a')
                _dataoutFrameDF = h.pd.DataFrame(_dataoutFrame, columns=['Stoffsystem'])

                _dataoutParamsSDF = h.pd.DataFrame(_dataoutParamsS, columns=['Params S'])
                _dataoutParamsRSDF = h.pd.DataFrame(_dataoutParamsRS, columns=['Params RS'])
                _dataoutParamsRDF = h.pd.DataFrame(_dataoutParamsR, columns=['Params R'])

                _dataout = [_dataoutFrameDF, _dataoutParamsSDF, _dataoutParamsRSDF, _dataoutParamsRDF]

                _output = h.pd.concat(_dataout, axis=1)
                return _output

            return solve()

    @staticmethod
    def _parameters_lak(datain):
        _alpha = 0.4

        _gabSM = -5145.15
        _gbaSM = 4077.76
        _gSaM = h.np.array([_gabSM, _gbaSM, _alpha])
        _gabRM = -5031.57
        _gbaRM = 4787.22
        _gRaM = h.np.array([_gabRM, _gbaRM, _alpha])

        _gabSE = -4883.75
        _gbaSE = 12111.15
        _gSaE = h.np.array([_gabSE, _gbaSE, _alpha])
        _gabRE = -4833.93
        _gbaRE = 10744.97
        _gRaE = h.np.array([_gabRE, _gbaRE, _alpha])

        _gabSP = -5013.79
        _gbaSP = 5442.11
        _gSaP = h.np.array([_gabSP, _gbaSP, _alpha])
        _gabRP = -4696.62
        _gbaRP = 4967.03
        _gRaP = h.np.array([_gabRP, _gbaRP, _alpha])

        _gabSB = -4621.74
        _gbaSB = 4763.98
        _gSaB = h.np.array([_gabSB, _gbaSB, _alpha])
        _gabRB = -4345.53
        _gbaRB = 4545.06
        _gRaB = h.np.array([_gabRB, _gbaRB, _alpha])

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

        _bnds = ((-6000.0, None), (None, 8000.0), (0.1, 1.0))

        def _methyl():
            steps_x = _xsM.count()
            steps_t = _tsM.count()

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaM, args=(_xsM, _tsM, _h0S, _t0S, steps_x), method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_S.fun, _tsM)))))
            print('deltag_AB s [ ] = ' + str(res_S.x[0]))
            print('deltag_BA s [ ] = ' + str(res_S.x[1]))
            print('alpha s =', res_S.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_tsM[x], _xsM[x], _h0S, _t0S, _gSaM)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_tsM[x] - _tcalc[x]) / _tsM[x])
                t_diff[x] = abs((_tsM[x] - _tcalc[x]))
            ard_neu_norm = 100 / len(_tsM) * sum(t_diff_norm)
            ard_neu = 100 / len(_tsM) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ard_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ard_neu)

            steps_x = _xrM.count()
            steps_t = _trM.count()

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaM, args=(_xrM, _trM, _h0R, _t0R, steps_x),
                                            method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_R.fun, _trM)))))
            print('deltag_AB r [ ] = ' + str(res_R.x[0]))
            print('deltag_BA r [ ] = ' + str(res_R.x[1]))
            print('alpha r =', res_R.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_trM[x], _xrM[x], _h0R, _t0R, _gRaM)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_trM[x] - _tcalc[x]) / _trM[x])
                t_diff[x] = abs((_trM[x] - _tcalc[x]))
            ard_neu_norm = 100 / len(_trM) * sum(t_diff_norm)
            ard_neu = 100 / len(_trM) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ard_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ard_neu)
            _params = (res_S.x[0], res_S.x[1], res_S.x[2], res_R.x[0], res_R.x[1], res_R.x[2])
            return _params

        def _ethyl():
            steps_x = _xsE.count()
            steps_t = _tsE.count()

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaE, args=(_xsE, _tsE, _h0S, _t0S, steps_x), method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_S.fun, _tsE)))))
            print('deltag_AB s [ ] = ' + str(res_S.x[0]))
            print('deltag_BA s [ ] = ' + str(res_S.x[1]))
            print('alpha s =', res_S.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_tsE[x], _xsE[x], _h0S, _t0S, _gSaE)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_tsE[x] - _tcalc[x]) / _tsE[x])
                t_diff[x] = abs((_tsE[x] - _tcalc[x]))
            ARD_neu_norm = 100 / len(_tsE) * sum(t_diff_norm)
            ARD_neu = 100 / len(_tsE) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ARD_neu)

            steps_x = _xrE.count()
            steps_t = _trE.count()

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaE, args=(_xrE, _trE, _h0R, _t0R, steps_x),
                                            method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_R.fun, _trE)))))
            print('deltag_AB r [ ] = ' + str(res_R.x[0]))
            print('deltag_BA r [ ] = ' + str(res_R.x[1]))
            print('alpha r =', res_R.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_trE[x], _xrE[x], _h0R, _t0R, _gRaE)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_trE[x] - _tcalc[x]) / _trE[x])
                t_diff[x] = abs((_trE[x] - _tcalc[x]))
            ARD_neu_norm = 100 / len(_trE) * sum(t_diff_norm)
            ARD_neu = 100 / len(_trE) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ARD_neu)
            _params = (res_S.x[0], res_S.x[1], res_S.x[2], res_R.x[0], res_R.x[1], res_R.x[2])
            return _params

        def _propyl():
            steps_x = _xsP.count()
            steps_t = _tsP.count()

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaP, args=(_xsP, _tsP, _h0S, _t0S, steps_x), method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_S.fun, _tsP)))))
            print('deltag_AB s [ ] = ' + str(res_S.x[0]))
            print('deltag_BA s [ ] = ' + str(res_S.x[1]))
            print('alpha s =', res_S.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_tsP[x], _xsP[x], _h0S, _t0S, _gSaP)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_tsP[x] - _tcalc[x]) / _tsP[x])
                t_diff[x] = abs((_tsP[x] - _tcalc[x]))
            ARD_neu_norm = 100 / len(_tsP) * sum(t_diff_norm)
            ARD_neu = 100 / len(_tsP) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ARD_neu)

            steps_x = _xrP.count()
            steps_t = _trP.count()

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaP, args=(_xrP, _trP, _h0R, _t0R, steps_x),
                                            method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_R.fun, _trP)))))
            print('deltag_AB r [ ] = ' + str(res_R.x[0]))
            print('deltag_BA r [ ] = ' + str(res_R.x[1]))
            print('alpha r =', res_R.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_trP[x], _xrP[x], _h0R, _t0R, _gRaP)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_trP[x] - _tcalc[x]) / _trP[x])
                t_diff[x] = abs((_trP[x] - _tcalc[x]))
            ARD_neu_norm = 100 / len(_trP) * sum(t_diff_norm)
            ARD_neu = 100 / len(_trP) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ARD_neu)
            _params = (res_S.x[0], res_S.x[1], res_S.x[2], res_R.x[0], res_R.x[1], res_R.x[2])
            return _params

        def _butyl():
            steps_x = _xsB.count()
            steps_t = _tsB.count()

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaB, args=(_xsB, _tsB, _h0S, _t0S, steps_x), method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_S.fun, _tsB)))))
            print('deltag_AB s [ ] = ' + str(res_S.x[0]))
            print('deltag_BA s [ ] = ' + str(res_S.x[1]))
            print('alpha s =', res_S.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_tsB[x], _xsB[x], _h0S, _t0S, _gSaB)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_tsB[x] - _tcalc[x]) / _tsB[x])
                t_diff[x] = abs((_tsB[x] - _tcalc[x]))
            ARD_neu_norm = 100 / len(_tsB) * sum(t_diff_norm)
            ARD_neu = 100 / len(_tsB) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ARD_neu)

            steps_x = _xrB.count()
            steps_t = _trB.count()

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaB, args=(_xrB, _trB, _h0R, _t0R, steps_x),
                                            method='Nelder-Mead',bounds=_bnds)
            print('ARD [%] = ' + str(100 * h.np.sum(h.np.abs(h.np.divide(-res_R.fun, _trB)))))
            print('deltag_AB r [ ] = ' + str(res_R.x[0]))
            print('deltag_BA r [ ] = ' + str(res_R.x[1]))
            print('alpha r =', res_R.x[2])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)
            for x in range(steps_t):
                _tcalc[x] = NrtlFit.t_sle(_trB[x], _xrB[x], _h0R, _t0R, _gRaB)

            for x in range(steps_t):
                t_diff_norm[x] = abs((_trB[x] - _tcalc[x]) / _trB[x])
                t_diff[x] = abs((_trB[x] - _tcalc[x]))
            ARD_neu_norm = 100 / len(_trB) * sum(t_diff_norm)
            ARD_neu = 100 / len(_trB) * sum(t_diff)
            print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
            print('ARD_neu_bzgl.T [%] =', ARD_neu)
            _params = (res_S.x[0], res_S.x[1], res_S.x[2], res_R.x[0], res_R.x[1], res_R.x[2])
            return _params

        def _solve():
            _methyl()
            _ethyl()
            _propyl()
            _butyl()
            _dataoutParamsM = _methyl()
            _dataoutParamsE = _ethyl()
            _dataoutParamsP = _propyl()
            _dataoutParamsB = _butyl()

            _dataoutFrame = ('gAB_S','gBA_S','a_S','gAB_R','gBA_R','a_R')
            _dataoutFrameDF = h.pd.DataFrame(_dataoutFrame, columns=['Stoffsystem'])

            _dataoutParamsM = h.pd.DataFrame(_dataoutParamsM, columns=['Params M'])
            _dataoutParamsE = h.pd.DataFrame(_dataoutParamsE, columns=['Params E'])
            _dataoutParamsP = h.pd.DataFrame(_dataoutParamsP, columns=['Params P'])
            _dataoutParamsB = h.pd.DataFrame(_dataoutParamsB, columns=['Params B'])

            _dataout = [_dataoutFrameDF, _dataoutParamsM, _dataoutParamsE, _dataoutParamsP, _dataoutParamsB]

            _output = h.pd.concat(_dataout, axis=1)

            return _output
        return _solve()
