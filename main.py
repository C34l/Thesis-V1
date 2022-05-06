import pandas as pd

import headers as h


#INPUT-Lappi:- C:\\Users\\Admin\\Desktop\\Input.xlsx
#OUTPUT-Lappi:- C:\\Users\\Admin\\Desktop\\Output.xlsx

#INPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Input.xlsx
#OUTPUT-Desk:- C:\\Users\\Ulf\\Desktop\\Output.xlsx

df = h.rw.IO._read()

output = h.pd.DataFrame(h.nr.Nrtl.approximategamma(df))
h.rw.IO._write(output)
