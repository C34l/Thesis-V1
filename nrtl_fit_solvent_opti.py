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

class params:
    _paramset1 = h.np.array([])
    _paramset2 = h.np.array([])
    _paramset3 = h.np.array([])
    _params = h.np.array([_paramset1, _paramset2, _paramset3])

    _h0 = ([24500.0,25600.0,25736.0])
    _t0 = ([404.75,393.35,404.65])

    _alphaabS = 0.3
    _alphabaS = 0.3
    _gabS = -2000.0
    _gbaS = 8000.0
    # _gSa = h.np.array([_gabS, _gbaS, _alphaabS, _alphabaS])
    _gSa = h.np.array([_gabS, _gbaS, ])
    _t0S = 404.75
    _h0S = 24500.0
    _xs = datain["x - SWasser"]
    _tSexp = datain["T - SWasser"]
    _lens = _tSexp.count()
    _lens1 = _lens.astype(int)
    _ts = h.np.array(_tSexp[0:_lens1])
    # objekt beinhaltet: t0S, h0S, gab, gba, alpha, texp, xi,

    @staticmethod
    def _paramfill(datain,_params, _h0, _t0):
        _steps = len(_params)
        for x in range (_steps):
            _params[x] =

        return 0
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