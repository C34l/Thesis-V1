import headers as h


#INPUT:- C:\\Users\\Admin\\Desktop\\Input.xlsx
#OUTPUT:- C:\\Users\\Admin\\Desktop\\Output.xlsx

df = h.rw.IO._read()
print(df)
test = df["x+"]
#print(df.iloc[1:1,:])
print(test)
h.nr.Nrtl.approximate(df)
