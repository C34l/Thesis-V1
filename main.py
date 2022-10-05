import pandas as pd

import headers as h

# INPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Input.xlsx
# OUTPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Output_pure.xlsx
# OUTPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Output_h2o.xlsx

# INPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Input.xlsx

# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_pure.xlsx

# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_h2o.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_h2o_param.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_h2o_refit.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_h2o_param1a.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_h2o_param2a.xlsx

# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_LAK.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_LAK_param.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_LAK_param1a.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_LAK_param2a.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_LAK_refit.xlsx

# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_Programmtest.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Calc_Eut.xlsx

#df = h.rw.IO._read()

# output = h.pd.DataFrame(h.nrp.Nrtl.binary_pure_ma(df))
# output = h.pd.DataFrame(h.nrp.Nrtl.ma_water(df))
# output = h.pd.DataFrame(h.nrp.Nrtl.ma_alkyl_lactate(df))

#output = h.pd.DataFrame(h.nrfs.NrtlFit.parametersH2O(df))
# output = h.pd.DataFrame(h.nrfs.NrtlFit._parameters_lak(df))

# output = h.pd.DataFrame(h.nrpe.Nrtl.ma_water(df))
# output = h.pd.DataFrame(h.nrpe.Nrtl.ma_alkyl_lactate(df))

# output = h.pd.DataFrame(h.nrfs1a.NrtlFit.parametersH2O(df))
# output = h.pd.DataFrame(h.nrfs1a.NrtlFit._parameters_lak(df))

# output = h.pd.DataFrame(h.nrfs2a.NrtlFit.parametersH2O(df))
# output = h.pd.DataFrame(h.nrfs2a.NrtlFit._parameters_lak(df))

# output = h.pd.DataFrame(h.nrtlt.NrtlFit.parametersH2O(df))
#output = h.pd.DataFrame(h.Eut.EutFind.pure_ma())

#output = h.Eut.EutFind.pure_ma()
#output = h.bt.batterytest.test()

output = h.diag.Diagrams.Ansatz_C()

#h.rw.IO._write(output)


