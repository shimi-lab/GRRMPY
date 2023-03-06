from xlsxwriter.worksheet import Worksheet
from pandas import DataFrame

class Worksheet(Worksheet):
    def __init__(self):
        super().__init__()
        
    def add_table(self, first_row=None, first_col=None, last_row=None, last_col=None,options=None):
        """tableを挿入する
        
        | 基本はxlsxwriterのadd_tableと同じ.
        | 第一引数にDataFrameをとった時のみA1の位置にTableを作成する
        """
        if type(first_row) == DataFrame:
            df = first_row
            data = df.values.tolist()
            last_row = len(data)
            last_col = len(data[0]) - 1
            columns = [{'header': text} for text in df.columns]
            option_dict = {"data":data,"columns":columns}
            if options:
                option_dict.update(options)
            super().add_table(0,0,last_row,last_col,option_dict)
        else:
            super().add_table(first_row, first_col, last_row, last_col,options)