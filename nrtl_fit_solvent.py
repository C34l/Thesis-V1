
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
            _gabS = -3000.0
            _gbaS = 9000.0
            # _gSa = h.np.array([_gabS, _gbaS, _alphaabS, _alphabaS])
            _gSa = h.np.array([_gabS, _gbaS, ])
            _alphaabRS = 0.3
            _alphabaRS = 0.3
            _gabRS = -3000.0
            _gbaRS = 9000.0
            # _gRSa = h.np.array([_gabRS, _gbaRS, _alphaabRS, _alphabaRS])
            _gRSa = h.np.array([_gabRS, _gbaRS,])
            _alphaabR = 0.4
            _alphabaR = 0.4
            _gabR = -3500
            _gbaR = 9800.0
            # _gRa = h.np.array([_gabR, _gbaR, _alphaabR, _alphabaR])
            _gRa = h.np.array([_gabR, _gbaR,])

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

            _bnds1 = ((None, None), (0, None),)
            _bnds2 = ((None, None), (0, None),)
            _bnds3 = ((None, None), (0, None),)

            def s_ma():
                steps_t = len(_ts)

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gSa, args=(_xs, _ts, _h0S, _t0S, steps_t,_Alpha), method='Powell', bounds=_bnds1)
                print('deltag_AB s = ' + str(res_1.x[0]))
                print('deltag_BA s = ' + str(res_1.x[1]))
                # print('alpha ab s =', res_1.x[2])
                # print('alpha ba s =', res_1.x[3])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                _gamma = h.np.zeros(steps_t)

                # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _gneu = (res_1.x[0], res_1.x[1],)

                for x in range(steps_t):
                    _tcalc1 = h.spo.fsolve(NrtlFit.t_sle, _ts[x], args=(_xs[x], _h0S, _t0S, _gneu, _Alpha), full_output=True)
                    _tcalc[x] = _tcalc1[0]
                    t_diff[x] = abs((_ts[x] - _tcalc[x]))
                    _gamma[x] = h.np.exp(NrtlFit.y_nrtl(_gneu, _xs[x], _tcalc[x], _Alpha))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _ts))
                ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (res_1.x[0], res_1.x[1], ard_neu_norm, _tcalc, _gamma, )
                return _params

            def rs_ma():
                steps_t = len(_trs)

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gRSa, args=(_xrs, _trs, _h0RS, _t0RS, steps_t, _Alpha),
                                                method='Powell', bounds=_bnds2)
                print('deltag_AB s = ' + str(res_1.x[0]))
                print('deltag_BA s = ' + str(res_1.x[1]))
                # print('alpha ab s =', res_1.x[2])
                # print('alpha ba s =', res_1.x[3])

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                _gamma = h.np.zeros(steps_t)

                # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _gneu = (res_1.x[0], res_1.x[1],)

                for x in range(steps_t):
                    _tcalc[x] = h.spo.fsolve(NrtlFit.t_sle, _trs[x], args=(_xrs[x], _h0RS, _t0RS, _gneu, _Alpha))
                    t_diff[x] = abs((_trs[x] - _tcalc[x]))
                    _gamma[x] = h.np.exp(NrtlFit.y_nrtl(_gneu, _xrs[x], _tcalc[x], _Alpha))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _trs))
                ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (res_1.x[0], res_1.x[1], ard_neu_norm, _tcalc, _gamma,)
                return _params

            def r_ma():
                steps_t = len(_tr)

                res_1 = scipy.optimize.minimize(NrtlFit.min_fqs, _gRa, args=(_xr, _tr, _h0R, _t0R, steps_t,_Alpha),
                                                method='Nelder-Mead', bounds=_bnds3)
                print('deltag_AB [ ] = ' + str(res_1.x[0]))
                print('deltag_BA [ ] = ' + str(res_1.x[1]))
                #print('alpha =', res_1.x[2])
                #print('alpha ba s =', res_1.x[3])

                # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _gneu = (res_1.x[0], res_1.x[1],)

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                _gamma = h.np.zeros(steps_t)

                for x in range(steps_t):
                    _tcalc[x] = h.spo.fsolve(NrtlFit.t_sle, _tr[x], args=(_xr[x], _h0R, _t0R, _gneu, _Alpha))
                    t_diff[x] = abs((_tr[x] - _tcalc[x]))
                    _gamma[x] = h.np.exp(NrtlFit.y_nrtl(_gneu, _xr[x], _tcalc[x], _Alpha))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _tr))
                ard_neu_norm = (100 / (steps_t)) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (res_1.x[0], res_1.x[1], ard_neu_norm, _tcalc, _gamma, )
                return _params

            def solve():
                _GABS, _GBAS, _ARDS, _TCALCS, _GAMMAS = s_ma()
                _GABRS, _GBARS, _ARDRS, _TCALCRS, _GAMMARS = rs_ma()
                _GABR, _GBAR, _ARDR, _TCALCR, _GAMMAR = r_ma()

                _dataoutParamsS = [_GABS, _GBAS, _ARDS]
                _dataoutParamsRS = [_GABRS, _GBARS, _ARDRS]
                _dataoutParamsR = [_GABR, _GBAR, _ARDR]


                _dataoutFrame = ('gAB', 'gBA', 'ARD', 'a')
                _dataoutFrameDF = h.pd.DataFrame(_dataoutFrame, columns=['Stoffsystem/ Temperaturen / Aktivitäten / ARD'])

                _dataoutParamsSDF = h.pd.DataFrame(_dataoutParamsS, columns=['Params S'])
                _dataoutParamsSDFt = h.pd.DataFrame(_TCALCS, columns=['T S'])
                _dataoutParamsSDFy = h.pd.DataFrame(_GAMMAS, columns=['y S'])

                _dataoutParamsRSDF = h.pd.DataFrame(_dataoutParamsRS, columns=['Params RS'])
                _dataoutParamsRSDFt = h.pd.DataFrame(_TCALCRS, columns=['T RS'])
                _dataoutParamsRSDFy = h.pd.DataFrame(_GAMMARS, columns=['y RS'])

                _dataoutParamsRDF = h.pd.DataFrame(_dataoutParamsR, columns=['Params R'])
                _dataoutParamsRDFt = h.pd.DataFrame(_TCALCR, columns=['T S'])
                _dataoutParamsRDFy = h.pd.DataFrame(_GAMMAR, columns=['y S'])

                _dataout = [_dataoutFrameDF,
                            _dataoutParamsSDF, _dataoutParamsSDFt, _dataoutParamsSDFy,
                            _dataoutParamsRSDF, _dataoutParamsRSDFt, _dataoutParamsRSDFy,
                            _dataoutParamsRDF, _dataoutParamsRDFt, _dataoutParamsRDFy]

                _output = h.pd.concat(_dataout, axis=1)
                return _output

            return solve()

    @staticmethod
    def _parameters_lak(datain):
        _alpha = 0.4

        _gabSM = -5145.15
        _gbaSM = 4077.76
        _gSaM = h.np.array([_gabSM, _gbaSM, ])
        _gabRM = -5031.57
        _gbaRM = 4787.22
        _gRaM = h.np.array([_gabRM, _gbaRM, ])

        _gabSE = -4883.75
        _gbaSE = 12111.15
        _gSaE = h.np.array([_gabSE, _gbaSE, ])
        _gabRE = -4833.93
        _gbaRE = 10744.97
        _gRaE = h.np.array([_gabRE, _gbaRE, ])

        _gabSP = -5013.79
        _gbaSP = 5442.11
        _gSaP = h.np.array([_gabSP, _gbaSP, ])
        _gabRP = -4696.62
        _gbaRP = 4967.03
        _gRaP = h.np.array([_gabRP, _gbaRP, ])

        _gabSB = -4621.74
        _gbaSB = 4763.98
        _gSaB = h.np.array([_gabSB, _gbaSB, ])
        _gabRB = -4345.53
        _gbaRB = 4545.06
        _gRaB = h.np.array([_gabRB, _gbaRB, ])

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

        _bnds = ((-6000.0, None), (None, None),)

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
            _gammaS = h.np.zeros(steps_tS)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuS = (res_S.x[0], res_S.x[1])

            for x in range(steps_tS):
                _tcalc1S = h.spo.fsolve(NrtlFit.t_sle, _tsm[x], args=(_xsM[x], _h0S, _t0S, _gneuS, _Alpha),
                                       full_output=True)
                _tcalcS[x] = _tcalc1S[0]
                t_diffS[x] = abs((_tsm[x] - _tcalcS[x]))
                _gammaS[x] = h.np.exp(NrtlFit.y_nrtl(_gneuS, _xsM[x], _tcalcS[x], _Alpha))

            t_diff_normS= h.np.abs(h.np.divide(t_diffS, _tsm))
            ard_neu_normS = (100 / steps_tS) * sum(t_diff_normS)
            print('ARD_normiert [%] =', ard_neu_normS)
            print(_tcalcS)

            steps_tR = len(_trm)

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaM, args=(_xrM, _trm, _h0R, _t0R, steps_tR, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB r = ' + str(res_R.x[0]))
            print('deltag_BA r = ' + str(res_R.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normR = h.np.zeros(steps_tR)
            t_diffR = h.np.zeros(steps_tR)
            _tcalcR = h.np.zeros(steps_tR)
            _gammaR = h.np.zeros(steps_tR)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuR = (res_R.x[0], res_R.x[1])

            for x in range(steps_tR):
                _tcalc1R = h.spo.fsolve(NrtlFit.t_sle, _trm[x], args=(_xrM[x], _h0R, _t0R, _gneuR, _Alpha),
                                        full_output=True)
                _tcalcR[x] = _tcalc1R[0]
                t_diffR[x] = abs((_trm[x] - _tcalcR[x]))
                _gammaR[x] = h.np.exp(NrtlFit.y_nrtl(_gneuR, _xrM[x], _tcalcR[x], _Alpha))

            t_diff_normR = h.np.abs(h.np.divide(t_diffR, _trm))
            ard_neu_normR = (100 / steps_tR) * sum(t_diff_normR)
            print('ARD_normiert [%] =', ard_neu_normR)
            print(_tcalcR)
            # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _params = (res_S.x[0], res_S.x[1], ard_neu_normS, _tcalcS, _gammaS, res_R.x[0], res_R.x[1], ard_neu_normR, _tcalcR, _gammaR)
            return _params

        def _ethyl():
            steps_tS = len(_tse)

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaE, args=(_xsE, _tse, _h0S, _t0S, steps_tS, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB s = ' + str(res_S.x[0]))
            print('deltag_BA s = ' + str(res_S.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normS = h.np.zeros(steps_tS)
            t_diffS = h.np.zeros(steps_tS)
            _tcalcS = h.np.zeros(steps_tS)
            _gammaS = h.np.zeros(steps_tS)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuS = (res_S.x[0], res_S.x[1])

            for x in range(steps_tS):
                _tcalc1S = h.spo.fsolve(NrtlFit.t_sle, _tse[x], args=(_xsE[x], _h0S, _t0S, _gneuS, _Alpha),
                                        full_output=True)
                _tcalcS[x] = _tcalc1S[0]
                t_diffS[x] = abs((_tse[x] - _tcalcS[x]))
                _gammaS[x] = h.np.exp(NrtlFit.y_nrtl(_gneuS, _xsE[x], _tcalcS[x], _Alpha))

            t_diff_normS = h.np.abs(h.np.divide(t_diffS, _tse))
            ard_neu_normS = (100 / steps_tS) * sum(t_diff_normS)
            print('ARD_normiert [%] =', ard_neu_normS)
            print(_tcalcS)

            steps_tR = len(_tre)

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaE, args=(_xrE, _tre, _h0R, _t0R, steps_tR, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB r = ' + str(res_R.x[0]))
            print('deltag_BA r = ' + str(res_R.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normR = h.np.zeros(steps_tR)
            t_diffR = h.np.zeros(steps_tR)
            _tcalcR = h.np.zeros(steps_tR)
            _gammaR = h.np.zeros(steps_tR)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuR = (res_R.x[0], res_R.x[1])

            for x in range(steps_tR):
                _tcalc1R = h.spo.fsolve(NrtlFit.t_sle, _tre[x], args=(_xrE[x], _h0R, _t0R, _gneuR, _Alpha),
                                        full_output=True)
                _tcalcR[x] = _tcalc1R[0]
                t_diffR[x] = abs((_tre[x] - _tcalcR[x]))
                _gammaR[x] = h.np.exp(NrtlFit.y_nrtl(_gneuR, _xrE[x], _tcalcR[x], _Alpha))

            t_diff_normR = h.np.abs(h.np.divide(t_diffR, _tre))
            ard_neu_normR = (100 / steps_tR) * sum(t_diff_normR)
            print('ARD_normiert [%] =', ard_neu_normR)
            print(_tcalcR)
            # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _params = (
            res_S.x[0], res_S.x[1], ard_neu_normS, _tcalcS, _gammaS, res_R.x[0], res_R.x[1], ard_neu_normR, _tcalcR,
            _gammaR)
            return _params

        def _propyl():
            steps_tS = len(_tsp)

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaP, args=(_xsP, _tsp, _h0S, _t0S, steps_tS, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB s = ' + str(res_S.x[0]))
            print('deltag_BA s = ' + str(res_S.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normS = h.np.zeros(steps_tS)
            t_diffS = h.np.zeros(steps_tS)
            _tcalcS = h.np.zeros(steps_tS)
            _gammaS = h.np.zeros(steps_tS)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuS = (res_S.x[0], res_S.x[1])

            for x in range(steps_tS):
                _tcalc1S = h.spo.fsolve(NrtlFit.t_sle, _tsp[x], args=(_xsP[x], _h0S, _t0S, _gneuS, _Alpha),
                                        full_output=True)
                _tcalcS[x] = _tcalc1S[0]
                t_diffS[x] = abs((_tsp[x] - _tcalcS[x]))
                _gammaS[x] = h.np.exp(NrtlFit.y_nrtl(_gneuS, _xsP[x], _tcalcS[x], _alpha))

            t_diff_normS = h.np.abs(h.np.divide(t_diffS, _tsp))
            ard_neu_normS = (100 / steps_tS) * sum(t_diff_normS)
            print('ARD_normiert [%] =', ard_neu_normS)
            print(_tcalcS)

            steps_tR = len(_trp)

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaP, args=(_xrP, _trp, _h0R, _t0R, steps_tR, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB r = ' + str(res_R.x[0]))
            print('deltag_BA r = ' + str(res_R.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normR = h.np.zeros(steps_tR)
            t_diffR = h.np.zeros(steps_tR)
            _tcalcR = h.np.zeros(steps_tR)
            _gammaR = h.np.zeros(steps_tR)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuR = (res_R.x[0], res_R.x[1])

            for x in range(steps_tR):
                _tcalc1R = h.spo.fsolve(NrtlFit.t_sle, _trp[x], args=(_xrP[x], _h0R, _t0R, _gneuR, _Alpha),
                                        full_output=True)
                _tcalcR[x] = _tcalc1R[0]
                t_diffR[x] = abs((_trp[x] - _tcalcR[x]))
                _gammaR[x] = h.np.exp(NrtlFit.y_nrtl(_gneuR, _xrP[x], _tcalcR[x], _Alpha))

            t_diff_normR = h.np.abs(h.np.divide(t_diffR, _trp))
            ard_neu_normR = (100 / steps_tR) * sum(t_diff_normR)
            print('ARD_normiert [%] =', ard_neu_normR)
            print(_tcalcR)
            # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _params = (
                res_S.x[0], res_S.x[1], ard_neu_normS, _tcalcS, _gammaS, res_R.x[0], res_R.x[1], ard_neu_normR, _tcalcR,
                _gammaR)
            return _params

        def _butyl():
            steps_tS = len(_tsb)

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gSaB, args=(_xsB, _tsb, _h0S, _t0S, steps_tS, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB s = ' + str(res_S.x[0]))
            print('deltag_BA s = ' + str(res_S.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normS = h.np.zeros(steps_tS)
            t_diffS = h.np.zeros(steps_tS)
            _tcalcS = h.np.zeros(steps_tS)
            _gammaS = h.np.zeros(steps_tS)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuS = (res_S.x[0], res_S.x[1])

            for x in range(steps_tS):
                _tcalc1S = h.spo.fsolve(NrtlFit.t_sle, _tsb[x], args=(_xsB[x], _h0S, _t0S, _gneuS, _Alpha),
                                        full_output=True)
                _tcalcS[x] = _tcalc1S[0]
                t_diffS[x] = abs((_tsb[x] - _tcalcS[x]))
                _gammaS[x] = h.np.exp(NrtlFit.y_nrtl(_gneuS, _xsB[x], _tcalcS[x], _Alpha))

            t_diff_normS = h.np.abs(h.np.divide(t_diffS, _tsb))
            ard_neu_normS = (100 / steps_tS) * sum(t_diff_normS)
            print('ARD_normiert [%] =', ard_neu_normS)
            print(_tcalcS)

            steps_tR = len(_trb)

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRaB, args=(_xrB, _trb, _h0R, _t0R, steps_tR, _Alpha),
                                            method='Nelder-Mead', bounds=_bnds)
            print('deltag_AB r = ' + str(res_R.x[0]))
            print('deltag_BA r = ' + str(res_R.x[1]))
            # print('alpha ab s =', res_1.x[2])
            # print('alpha ba s =', res_1.x[3])

            t_diff_normR = h.np.zeros(steps_tR)
            t_diffR = h.np.zeros(steps_tR)
            _tcalcR = h.np.zeros(steps_tR)
            _gammaR = h.np.zeros(steps_tR)

            # _gneu = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _gneuR = (res_R.x[0], res_R.x[1])

            for x in range(steps_tR):
                _tcalc1R = h.spo.fsolve(NrtlFit.t_sle, _trb[x], args=(_xrB[x], _h0R, _t0R, _gneuR, _Alpha),
                                        full_output=True)
                _tcalcR[x] = _tcalc1R[0]
                t_diffR[x] = abs((_trb[x] - _tcalcR[x]))
                _gammaR[x] = h.np.exp(NrtlFit.y_nrtl(_gneuR, _xrB[x], _tcalcR[x], _Alpha))

            t_diff_normR = h.np.abs(h.np.divide(t_diffR, _trb))
            ard_neu_normR = (100 / steps_tR) * sum(t_diff_normR)
            print('ARD_normiert [%] =', ard_neu_normR)
            print(_tcalcR)
            # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
            _params = (
                res_S.x[0], res_S.x[1], ard_neu_normS, _tcalcS, _gammaS, res_R.x[0], res_R.x[1], ard_neu_normR, _tcalcR,
                _gammaR)
            return _params

        def _solve():
            _GABSM, _GBASM, _ARDSM, _TCALCSM, _GAMMASM, _GABRM, _GBARM, _ARDRM, _TCALCRM, _GAMMARM = _methyl()
            _GABSE, _GBASE, _ARDSE, _TCALCSE, _GAMMASE, _GABRE, _GBARE, _ARDRE, _TCALCRE, _GAMMARE = _ethyl()
            _GABSP, _GBASP, _ARDSP, _TCALCSP, _GAMMASP, _GABRP, _GBARP, _ARDRP, _TCALCRP, _GAMMARP = _propyl()
            _GABSB, _GBASB, _ARDSB, _TCALCSB, _GAMMASB, _GABRB, _GBARB, _ARDRB, _TCALCRB, _GAMMARB = _butyl()

            _dataoutParamsSM = [_GABSM, _GBASM, _ARDSM]
            _dataoutParamsRM = [_GABRM, _GBARM, _ARDRM]
            _dataoutParamsSE = [_GABSE, _GBASE, _ARDSE]
            _dataoutParamsRE = [_GABRE, _GBARE, _ARDRE]
            _dataoutParamsSP = [_GABSP, _GBASP, _ARDSP]
            _dataoutParamsRP = [_GABRP, _GBARP, _ARDRP]
            _dataoutParamsSB = [_GABSB, _GBASB, _ARDSB]
            _dataoutParamsRB = [_GABRB, _GBARB, _ARDRB]


            _dataoutFrame = ('gAB','gBA','ARD','a')
            _dataoutFrameDF = h.pd.DataFrame(_dataoutFrame, columns=['Stoffsystem'])

            _dataoutParamsSMDF = h.pd.DataFrame(_dataoutParamsSM, columns=['Params SM'])
            _dataoutParamsRMDF = h.pd.DataFrame(_dataoutParamsRM, columns=['Params RM'])
            _TCALCSMDF = h.pd.DataFrame(_TCALCSM, columns=['T S M'])
            _GAMMASMDF = h.pd.DataFrame(_GAMMASM, columns=['y S M'])
            _TCALCRMDF = h.pd.DataFrame(_TCALCRM, columns=['T R M'])
            _GAMMARMDF = h.pd.DataFrame(_GAMMARM, columns=['y R M'])

            _dataoutParamsSEDF = h.pd.DataFrame(_dataoutParamsSE, columns=['Params SE'])
            _dataoutParamsREDF = h.pd.DataFrame(_dataoutParamsRE, columns=['Params RE'])
            _TCALCSEDF = h.pd.DataFrame(_TCALCSE, columns=['T S E'])
            _GAMMASEDF = h.pd.DataFrame(_GAMMASE, columns=['y S E'])
            _TCALCREDF = h.pd.DataFrame(_TCALCRE, columns=['T R E'])
            _GAMMAREDF = h.pd.DataFrame(_GAMMARE, columns=['y R E'])

            _dataoutParamsSPDF = h.pd.DataFrame(_dataoutParamsSP, columns=['Params SP'])
            _dataoutParamsRPDF = h.pd.DataFrame(_dataoutParamsRP, columns=['Params RP'])
            _TCALCSPDF = h.pd.DataFrame(_TCALCSP, columns=['T S P'])
            _GAMMASPDF = h.pd.DataFrame(_GAMMASP, columns=['y S P'])
            _TCALCRPDF = h.pd.DataFrame(_TCALCRP, columns=['T R P'])
            _GAMMARPDF = h.pd.DataFrame(_GAMMARP, columns=['y R P'])

            _dataoutParamsSBDF = h.pd.DataFrame(_dataoutParamsSB, columns=['Params SB'])
            _dataoutParamsRBDF = h.pd.DataFrame(_dataoutParamsRB, columns=['Params RB'])
            _TCALCSBDF = h.pd.DataFrame(_TCALCSB, columns=['T S B'])
            _GAMMASBDF = h.pd.DataFrame(_GAMMASB, columns=['y S B'])
            _TCALCRBDF = h.pd.DataFrame(_TCALCRB, columns=['T R B'])
            _GAMMARBDF = h.pd.DataFrame(_GAMMARB, columns=['y R B'])

            _dataout = [_dataoutFrameDF,
                        _dataoutParamsSMDF, _TCALCSMDF, _GAMMASMDF,
                        _dataoutParamsRMDF, _TCALCRMDF, _GAMMARMDF,
                        _dataoutParamsSEDF, _TCALCSEDF, _GAMMASEDF,
                        _dataoutParamsREDF, _TCALCREDF, _GAMMAREDF,
                        _dataoutParamsSPDF, _TCALCSPDF, _GAMMASPDF,
                        _dataoutParamsRPDF, _TCALCRPDF, _GAMMARPDF,
                        _dataoutParamsSBDF, _TCALCSBDF, _GAMMASBDF,
                        _dataoutParamsRBDF, _TCALCRBDF, _GAMMARBDF]

            _output = h.pd.concat(_dataout, axis=1)

            return _output
        return _solve()

    @staticmethod
    def parameters_Pure():
        def _rac():

            _alphaabRS = 0.4
            _alphabaRS = 0.4
            _gabRS = -3000.0
            _gbaRS = 9000.0
            _gRSa = h.np.array([_gabRS, _gbaRS, ])

            _xSRS = []
            _tSRS = []

            _xRRS = []
            _tRRS = []


            steps_tSRS = len(_tSRS)

            res_S = scipy.optimize.minimize(NrtlFit.min_fqs, _gRSa, args=(_xSRS, _tSRS, _h0RS, _t0RS, steps_tSRS, _Alpha),
                                            method='Nelder-Mead',)
            print('deltag_AB s = ' + str(res_S.x[0]))
            print('deltag_BA s = ' + str(res_S.x[1]))

            t_diff_normS = h.np.zeros(steps_tSRS)
            t_diffS = h.np.zeros(steps_tSRS)
            _tcalcS = h.np.zeros(steps_tSRS)
            _gammaS = h.np.zeros(steps_tSRS)

            _gneuS = (res_S.x[0], res_S.x[1])

            for x in range(steps_tSRS):
                _tcalc1S = h.spo.fsolve(NrtlFit.t_sle, _tSRS[x], args=(_xSRS[x], _h0RS, _t0RS, _gneuS, _Alpha),
                                        full_output=True)
                _tcalcS[x] = _tcalc1S[0]
                t_diffS[x] = abs((_tSRS[x] - _tcalcS[x]))
                _gammaS[x] = h.np.exp(NrtlFit.y_nrtl(_gneuS, _xSRS[x], _tcalcS[x], _Alpha))

            t_diff_normS = h.np.abs(h.np.divide(t_diffS, _tSRS))
            ard_neu_normS = (100 / steps_tSRS) * sum(t_diff_normS)
            print('ARD_normiert [%] =', ard_neu_normS)
            print(_tcalcS)

            steps_tRRS = len(_tRRS)

            res_R = scipy.optimize.minimize(NrtlFit.min_fqs, _gRSa, args=(_xRRS, _tRRS, _h0RS, _t0RS, steps_tRRS, _Alpha),
                                            method='Nelder-Mead',)
            print('deltag_AB r = ' + str(res_R.x[0]))
            print('deltag_BA r = ' + str(res_R.x[1]))


            t_diff_normR = h.np.zeros(steps_tRRS)
            t_diffR = h.np.zeros(steps_tRRS)
            _tcalcR = h.np.zeros(steps_tRRS)
            _gammaR = h.np.zeros(steps_tRRS)

            _gneuR = (res_R.x[0], res_R.x[1])

            for x in range(steps_tRRS):
                _tcalc1R = h.spo.fsolve(NrtlFit.t_sle, _tRRS[x], args=(_xRRS[x], _h0RS, _t0RS, _gneuR, _Alpha),
                                        full_output=True)
                _tcalcR[x] = _tcalc1R[0]
                t_diffR[x] = abs((_tRRS[x] - _tcalcR[x]))
                _gammaR[x] = h.np.exp(NrtlFit.y_nrtl(_gneuR, _xRRS[x], _tcalcR[x], _Alpha))

            t_diff_normR = h.np.abs(h.np.divide(t_diffR, _tRRS))
            ard_neu_normR = (100 / steps_tRRS) * sum(t_diff_normR)
            print('ARD_normiert [%] =', ard_neu_normR)
            print(_tcalcR)
            _params = (res_S.x[0], res_S.x[1], ard_neu_normS, _tcalcS, _gammaS, res_R.x[0], res_R.x[1], ard_neu_normR, _tcalcR,_gammaR)

            return 0

        return _rac()