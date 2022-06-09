
import scipy.optimize

import headers as h


_r = 8.31462618
_t0S = 404.75
_t0RS = 393.35
_t0R = 404.65
_h0S = 24500.0
_h0RS = 25600.0
_h0R = 25736.0
_Alpha = 0.4

class NrtlFit:
    @staticmethod
    def y_nrtl(_g, _xi, _texp, _alpha):
        _xj = float(1.0 - _xi)
        _tauij = _g[0] / (_r * _texp)
        _tauji = _g[1] / (_r * _texp)
        _Gij = h.np.exp(-_alpha * _tauij)
        _Gji = h.np.exp(-_alpha * _tauji)
        _func = ((_xj ** 2) * (_tauji * (_Gji / (_xi + _xj * _Gji)) ** 2 + _tauij * (_Gij / (_xi * _Gij + _xj) ** 2)))
        return _func

    @staticmethod
    def t_sle(_texp, _xi, _h0, _t0, _g, _alpha):
        _yi = NrtlFit.y_nrtl(_g, _xi, _texp, _alpha)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx)-1
        return _func

    @staticmethod
    def min_fqs(_g, _xi, _texp, _h0, _t0, _steps, _alpha):
        _tcalc = h.np.zeros(_steps)
        _tdiff = h.np.zeros(_steps)

        for x in range(_steps):
            _tcalc1 = h.spo.fsolve(NrtlFit.t_sle, _texp[x], args=(_xi[x], _h0, _t0, _g, _alpha), full_output=True)
            _tcalc[x] = _tcalc1[0]
            # _tmod[x] = scipy.optimize.fsolve(NrtlFit.t_sle, _tmod0, args=(_xi[x], _h0, _t0, _g))
            _tdiff[x] = _texp[x]-_tcalc[x]
        # fqs_terme = (_texp[0:steps] - _tmod) ** 2
        fqs_norm1 = (h.np.abs(h.np.divide(_tdiff, _texp)))
        fqs_norm = (h.np.abs(h.np.divide(_tdiff, _texp))) ** 2
        # fqs_norm = ((_texp[0:steps] - _tmod[0:steps]) / _texp[0:steps]) ** 2

        fqs_summe = h.np.sum(fqs_norm)

        return fqs_summe

    @staticmethod
    def parametersH2O(datain):
            _alpha = 0.4
            _alphaabS = 0.3
            _alphabaS = 0.3
            _gabS = -2000.0
            _gbaS = 8000.0
            # _gSa = h.np.array([_gabS, _gbaS, _alphaabS, _alphabaS])
            _gSa = h.np.array([_gabS, _gbaS, ])
            _alphaabRS = 0.3
            _alphabaRS = 0.3
            _gabRS = -3000.0
            _gbaRS = 9000.0
            # _gRSa = h.np.array([_gabRS, _gbaRS, _alphaabRS, _alphabaRS])
            _gRSa = h.np.array([_gabS, _gbaS, ])
            _alphaabR = 0.4
            _alphabaR = 0.4
            _gabR = -2000.0
            _gbaR = 8000.0
            # _gRa = h.np.array([_gabR, _gbaR, _alphaabR, _alphabaR])
            _gRa = h.np.array([_gabS, _gbaS, ])

            _xs = datain["x - SWasser"]
            _tSexp = datain["T - SWasser"]
            _lens = _tSexp.count()
            _lens1 = _lens.astype(int)
            _ts = h.np.array(_tSexp[0:_lens1])

            _xrs = datain["x - RSWasser"]
            _tRSexp = datain["T - RSWasser"]
            _lenrs = _tRSexp.count()
            _lenrs1 = _lenrs.astype(int)
            _trs = h.np.array(_tRSexp[0:_lenrs1])

            _xr = datain["x - RWasser"]
            _tRexp = datain["T - RWasser"]
            _lenr = _tRexp.count()
            _lenr1 = _lenr.astype(int)
            _tr = h.np.array(_tRexp[0:_lenr1])

            _bnds1 = ((None, None), (None, None),)
            _bnds2 = ((None, None), (None, None),)
            _bnds3 = ((None, None), (None, None),)

            def s_ma():
                steps_t = len(_ts)

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gSa, args=(_xs, _ts, _h0S, _t0S, steps_t, _Alpha), method='Nelder-Mead', bounds=_bnds1)
                print('deltag_AB s = ' + str(res_1.x[0]))
                print('deltag_BA s = ' + str(res_1.x[1]))
                # print('alpha ab s =', res_1.x[2])
                # print('alpha ba s =', res_1.x[3])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)

                # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _gneu = (res_1.x[0], res_1.x[1])

                for x in range(steps_t):
                    _tcalc1 = h.spo.fsolve(NrtlFit.t_sle, _ts[x], args=(_xs[x], _h0S, _t0S, _gneu, _Alpha), full_output=True)
                    _tcalc[x] = _tcalc1[0]
                    t_diff[x] = abs((_ts[x] - _tcalc[x]))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _ts))
                ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (res_1.x[0], res_1.x[1])
                return _params

            def rs_ma():
                steps_t = len(_trs)

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gRSa, args=(_xrs, _trs, _h0RS, _t0RS, steps_t, _Alpha),
                                                method='Nelder-Mead', bounds=_bnds2)
                print('deltag_AB s = ' + str(res_1.x[0]))
                print('deltag_BA s = ' + str(res_1.x[1]))
                # print('alpha ab s =', res_1.x[2])
                # print('alpha ba s =', res_1.x[3])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)

                # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _gneu = (res_1.x[0], res_1.x[1])

                for x in range(steps_t):
                    _tcalc[x] = h.spo.fsolve(NrtlFit.t_sle, _trs[x], args=(_xrs[x], _h0RS, _t0RS, _gneu, _Alpha))
                    t_diff[x] = abs((_trs[x] - _tcalc[x]))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _trs))
                ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (res_1.x[0], res_1.x[1])
                return _params

            def r_ma():
                steps_t = len(_tr)

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gRa, args=(_xr, _tr, _h0R, _t0R, steps_t, _Alpha),
                                                method='Powell', bounds=_bnds3)
                print('deltag_AB [ ] = ' + str(res_1.x[0]))
                print('deltag_BA [ ] = ' + str(res_1.x[1]))
                #print('alpha =', res_1.x[2])
                #print('alpha ba s =', res_1.x[3])

                # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _gneu = (res_1.x[0], res_1.x[1])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                for x in range(steps_t):
                    _tcalc[x] = h.spo.fsolve(NrtlFit.t_sle, _tr[x], args=(_xr[x], _h0R, _t0R, _gneu, _Alpha))

                for x in range(steps_t):
                    t_diff[x] = abs((_tr[x] - _tcalc[x]))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _tr))
                ard_neu_norm = (100 / (steps_t)) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (res_1.x[0], res_1.x[1])
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
        _gSaM = h.np.array([_gabSM, _gbaSM,])
        _gabRM = -5031.57
        _gbaRM = 4787.22
        _gRaM = h.np.array([_gabRM, _gbaRM,])

        _gabSE = -4883.75
        _gbaSE = 12111.15
        _gSaE = h.np.array([_gabSE, _gbaSE,])
        _gabRE = -4833.93
        _gbaRE = 10744.97
        _gRaE = h.np.array([_gabRE, _gbaRE,])

        _gabSP = -5013.79
        _gbaSP = 5442.11
        _gSaP = h.np.array([_gabSP, _gbaSP,])
        _gabRP = -4696.62
        _gbaRP = 4967.03
        _gRaP = h.np.array([_gabRP, _gbaRP,])

        _gabSB = -4621.74
        _gbaSB = 4763.98
        _gSaB = h.np.array([_gabSB, _gbaSB,])
        _gabRB = -4345.53
        _gbaRB = 4545.06
        _gRaB = h.np.array([_gabRB, _gbaRB,])

        _xsM = datain["x - SLakM"]
        _tsM = datain["T - SLakM"]
        _lensm = _tsM.count()
        _lensm1 = _lensm.astype(int)
        _tsm = h.np.array(_tsM[0:_lensm1])

        _xrM = datain["x - RLakM"]
        _trM = datain["T - RLakM"]
        _lenrm = _trM.count()
        _lenrm1 = _lenrm.astype(int)
        _trm = h.np.array(_trM[0:_lenrm1])

        _xsE = datain["x - SLakE"]
        _tsE = datain["T - SLakE"]
        _lense = _tsE.count()
        _lense1 = _lense.astype(int)
        _tse = h.np.array(_tsE[0:_lense1])

        _xrE = datain["x - RLakE"]
        _trE = datain["T - RLakE"]
        _lenre = _trE.count()
        _lenre1 = _lenre.astype(int)
        _tre = h.np.array(_trE[0:_lenre1])

        _xsP = datain["x - SLakP"]
        _tsP = datain["T - SLakP"]
        _lensp = _tsP.count()
        _lensp1 = _lensp.astype(int)
        _tsp = h.np.array(_tsP[0:_lensp1])

        _xrP = datain["x - RLakP"]
        _trP = datain["T - RLakP"]
        _lenrp = _trP.count()
        _lenrp1 = _lenrp.astype(int)
        _trp = h.np.array(_trP[0:_lenrp1])

        _xsB = datain["x - SLakB"]
        _tsB = datain["T - SLakB"]
        _lensb = _tsB.count()
        _lensb1 = _lensb.astype(int)
        _tsb = h.np.array(_tsB[0:_lensb1])

        _xrB = datain["x - RLakB"]
        _trB = datain["T - RLakB"]
        _lenrb = _trB.count()
        _lenrb1 = _lenrb.astype(int)
        _trb = h.np.array(_trB[0:_lenrb1])

        _bnds = ((-6000.0, None), (None, 8000.0), (0.1, 1.0))

        def _methyl():
            steps_tS = len(_tsm)

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaM, args=(_xsM, _tsm, _h0S, _t0S, steps_tS, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB s = ' + str(res_S.x[0]))
            print('deltag_BA s = ' + str(res_S.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normS = h.np.zeros(steps_tS)
            t_diffS = h.np.zeros(steps_tS)
            _tcalcS = h.np.zeros(steps_tS)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneu = (res_S.x[0], res_S.x[1])

            for x in range(steps_tS):
                _tcalc1S = h.spo.fsolve(NrtlFit.t_sle, _tsm[x], args=(_xsM[x], _h0S, _t0S, _gneu, _Alpha),
                                       full_output=True)
                _tcalcS[x] = _tcalc1S[0]
                t_diffS[x] = abs((_tsm[x] - _tcalcS[x]))

            t_diff_normS= h.np.abs(h.np.divide(t_diffS, _tsm))
            ard_neu_normS = (100 / steps_tS) * sum(t_diff_normS)
            print('ARD_normiert [%] =', ard_neu_normS)
            print(_tcalcS)

            steps_tR = len(_trm)

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaM, args=(_xrM, _tsm, _h0R, _t0R, steps_tR, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB r = ' + str(res_R.x[0]))
            print('deltag_BA r = ' + str(res_R.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normR = h.np.zeros(steps_tS)
            t_diffR = h.np.zeros(steps_tR)
            _tcalcR = h.np.zeros(steps_tR)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneu = (res_R.x[0], res_R.x[1])

            for x in range(steps_tR):
                _tcalc1R = h.spo.fsolve(NrtlFit.t_sle, _trm[x], args=(_xrM[x], _h0R, _t0R, _gneu, _Alpha),
                                        full_output=True)
                _tcalcR[x] = _tcalc1R[0]
                t_diffR[x] = abs((_trm[x] - _tcalcR[x]))

            t_diff_normR = h.np.abs(h.np.divide(t_diffR, _trm))
            ard_neu_normR = (100 / steps_tR) * sum(t_diff_normR)
            print('ARD_normiert [%] =', ard_neu_normR)
            print(_tcalcR)
            # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _params = (res_S.x[0], res_S.x[1], res_R.x[0], res_R.x[1])
            return _params

        def _ethyl():
            steps_t = len(_tse)

            res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaE, args=(_xsE, _tse, _h0S, _t0S, steps_t, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB s = ' + str(res_1.x[0]))
            print('deltag_BA s = ' + str(res_1.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_norm = h.np.zeros(steps_t)
            t_diff = h.np.zeros(steps_t)
            _tcalc = h.np.zeros(steps_t)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneu = (res_1.x[0], res_1.x[1])

            for x in range(steps_t):
                _tcalc1 = h.spo.fsolve(NrtlFit.t_sle, _tsm[x], args=(_xsM[x], _h0S, _t0S, _gneu, _Alpha),
                                       full_output=True)
                _tcalc[x] = _tcalc1[0]
                t_diff[x] = abs((_tsm[x] - _tcalc[x]))

            t_diff_norm = h.np.abs(h.np.divide(t_diff, _tsm))
            ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
            print('ARD_normiert [%] =', ard_neu_norm)
            print(_tcalc)
            # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _params = (res_1.x[0], res_1.x[1])
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
