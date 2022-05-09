import pandas as pd

import headers as h

# INPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Input.xlsx
# OUTPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Output_pure.xlsx
# OUTPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Output_h2o.xlsx

# INPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Input.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_pure.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output_h2o.xlsx

df = h.rw.IO._read()

# output = h.pd.DataFrame(h.nrp.Nrtl.binary_pure_ma(df))
output = h.pd.DataFrame(h.nrp.Nrtl.ma_water(df))
h.rw.IO._write(output)
