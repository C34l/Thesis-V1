import headers as h
import scipy.optimize

_r = 8.3141
_t0S = 380
_t0RS = 380
_t0R = 380
_h0S = 25000
_h0RS = 25000.0
_h0R = 25000.0
_Alpha = 0.4


class NrtlFit:
    @staticmethod
    def y_nrtl(_g, _xi, _texp, _alpha):

        _xj = float(1.0 - _xi)
        _tauij = (_g[0] / (_r * _texp))
        _tauji = (_g[1] / (_r * _texp))
        _Gij = (h.math.exp(-_alpha * _tauij))
        _Gji = (h.math.exp(-_alpha * _tauji))
        _func = ((_xj ** 2) * (_tauji * (_Gji / (_xi + _xj * _Gji)) ** 2 + _tauij * (_Gij / (_xi * _Gij + _xj) ** 2)))

        return _func

    @staticmethod
    def t_sle(_texp, _xi, _h0, _t0, _g, _alpha):

        _yi = NrtlFit.y_nrtl(_g, _xi, _texp, _alpha)
        _gx = h.np.exp(_yi) * _xi
        _func = (h.np.exp((-_h0 / (_r * _texp)) * (1 - (_texp / _t0))) / _gx)-1
        return _func

    @staticmethod
    def parametersH2O(datain):

            _gabS = -3000.0
            _gbaS = 9000
            _gSa = h.np.array([_gabS, _gbaS, ])

            _gabRS = -3000.0
            _gbaRS = 9000
            _gRSa = h.np.array([_gabRS, _gbaRS, ])

            _gabR = -3000
            _gbaR = 9000
            _gRa = h.np.array([_gabR, _gbaR, ])

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

            def s_ma():
                steps_t = len(_ts)

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                _gamma = h.np.zeros(steps_t)

                for x in range(steps_t):

                    _tcalc1 = h.spo.fsolve(NrtlFit.t_sle, _ts[x], args=(_xs[x], _h0S, _t0S, _gSa, _Alpha), full_output=True)
                    _tcalc[x] = _tcalc1[0]
                    t_diff[x] = abs((_ts[x] - _tcalc[x]))
                    _gamma[x] = h.math.exp(NrtlFit.y_nrtl(_gSa, _xs[x], _tcalc[x], _Alpha))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _ts))
                ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                _params = (_gSa[0], _gSa[1], ard_neu_norm, _tcalc, _gamma, )
                return _params

            def rs_ma():
                steps_t = len(_trs)

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                _gamma = h.np.zeros(steps_t)

                for x in range(steps_t):
                    _tcalc[x] = h.spo.fsolve(NrtlFit.t_sle, _trs[x], args=(_xrs[x], _h0RS, _t0RS, _gRSa, _Alpha))
                    t_diff[x] = abs((_trs[x] - _tcalc[x]))
                    _gamma[x] = h.np.exp(NrtlFit.y_nrtl(_gRSa, _xrs[x], _tcalc[x], _Alpha))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _trs))
                ard_neu_norm = (100 / steps_t) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (_gRSa[0], _gRSa[1], ard_neu_norm, _tcalc, _gamma, )
                return _params

            def r_ma():
                steps_t = len(_tr)

                t_diff_norm = h.np.zeros(steps_t)
                t_diff = h.np.zeros(steps_t)
                _tcalc = h.np.zeros(steps_t)
                _gamma = h.np.zeros(steps_t)

                for x in range(steps_t):
                    _tcalc[x] = h.spo.fsolve(NrtlFit.t_sle, _tr[x], args=(_xr[x], _h0R, _t0R, _gRa, _Alpha))
                    t_diff[x] = abs((_tr[x] - _tcalc[x]))
                    _gamma[x] = h.np.exp(NrtlFit.y_nrtl(_gRa, _xr[x], _tcalc[x], _Alpha))

                t_diff_norm = h.np.abs(h.np.divide(t_diff, _tr))
                ard_neu_norm = (100 / (steps_t)) * sum(t_diff_norm)
                print('ARD_normiert [%] =', ard_neu_norm)
                print(_tcalc)
                # _params = (res_1.x[0], res_1.x[1], res_1.x[2], res_1.x[3])
                _params = (_gRa[0], _gRa[1], ard_neu_norm, _tcalc, _gamma, )
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