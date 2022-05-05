import pandas as pd

import headers as h


#INPUT:- C:\\Users\\Admin\\Desktop\\Input.xlsx
#OUTPUT:- C:\\Users\\Admin\\Desktop\\Output.xlsx

df = h.rw.IO._read()
#print(df)
#test = df["x+"]
#print(df.iloc[1:1,:])
#print(test)
#h.nr.Nrtl.approximate(df)
output = h.pd.DataFrame(h.nr.Nrtl.approximategamma(df))
#output = pd.DataFrame(columns=['gamma'])
#h.rw.IO._write(output)
