import headers as h

class IO:
    def _read():
        _InputFilePath = input("Dateipfad des Inputs eingeben")
        _df = h.pd.read_excel((_InputFilePath), skiprows = 0, usecols = [0, 1] )
        return _df

    def _write(_df):
        _OutputFilePath = input("Dateipfad für Output eingeben")
        _df.to_excel(_OutputFilePath)
        return 0
