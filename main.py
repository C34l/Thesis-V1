import pandas as pd

import headers as h

# INPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Input.xlsx
# OUTPUT-Laptop:- C:\\Users\\Admin\\Desktop\\Output.xlsx

# INPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Input.xlsx
# OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output.xlsx

df = h.rw.IO._read()

output = h.pd.DataFrame(h.nrp.Nrtl.binarypureMA(df))
h.rw.IO._write(output)
