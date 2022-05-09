import headers as h

# #usecols=[0, 1,2,3,4,5,6]
# read from and write to properly formatted Excel files
class IO:
    @staticmethod
    def _read():
        _InputFilePath = input("Dateipfad des Inputs eingeben")
        _df = h.pd.read_excel(_InputFilePath, skiprows=0, dtype={0: float, 1: float})
        return _df

    @staticmethod
    def _write(_df):
        _OutputFilePath = input("Dateipfad für Output eingeben")
        _df.to_excel(_OutputFilePath)
        return 0
