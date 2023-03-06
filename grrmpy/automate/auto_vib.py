from ase.vibrations import Vibrations
from pathlib import Path
from tqdm.notebook import tqdm_notebook as tqdm
from plotly.subplots import make_subplots
# USER
from grrmpy.calculator import pfp_calculator
from grrmpy.vibrations.functions import get_vibdf,_to_table_and_imode

class AutoVib():
    def __init__(self,ts_list,indices=None,calc_func=pfp_calculator,name="vib",**kwargs):
        """振動数計算を行なう

        Paramters:
        
        ts_list: list of Atoms
            | 計算するTSのAtomsのリスト.
            | listの要素にNoneを指定する事もできる.Noneを指定した場合計算されない
        indice: list of int
            動かす原子のindex番号
        calc_func: functions object
            claculatorを返す関数
        name:
            指定したフォルダ名にキャッシュが作成される
        """
        self.ts_list = ts_list
        self.calc_func = calc_func
        for ts in self.ts_list:
            if not ts is None:
                ts.calc = self.calc_func()
        self.indices = indices
        self.folder = Path(name)
        self.kwargs = kwargs
        if not self.folder.exists():
            self.folder.mkdir()
        self.ivib = [] # 虚振動のindex番号のリストが入る
        self.vib = [] # 計算が終わったVibrationsオブジェクトのリストが入る
        
    def irun(self,ts,name,**kwargs):
        vib = Vibrations(ts,self.indices,name=f"{self.folder}/{name}") 
        vib.run()
        return vib
        
    def run(self):
        """計算実行"""
        for i,ts in enumerate(tqdm(self.ts_list)):
            if ts is None:
                self.vib.append(None)
                continue
            vib = self.irun(ts,i,**self.kwargs)
            self.vib.append(vib)
    
    def summary(self,html="Summary.html",cols=3):
        """
        振動数の表をhtmlへ出力する
        
        Parameters:
        
        html:str
            出力先のhtml
        """
        table_and_i = [_to_table_and_imode(vib) if vib else None
                        for vib in self.vib]
        splited_data = [table_and_i[i: i+cols] for i in range(0, len(table_and_i), cols)]
        subplot_titles = [f"TS{i}" for i in range(len(table_and_i))]
        
        rows = len(table_and_i)//cols + 1
        specs = [[{"type": "table"} for _ in range(cols)] for _ in range(rows)]
        fig = make_subplots(rows=rows, cols=cols,specs=specs,subplot_titles=subplot_titles)
        for i,row_data in enumerate(splited_data,start=1):
            for j,vib_data in enumerate(row_data,start=1):
                if vib_data is None:
                    self.ivib.append(None)
                    continue
                fig.add_trace(vib_data[0],row=i, col=j)
                self.ivib.append(vib_data[1]) #虚振動の振動モード番号
                
        fig.update_layout(height=700*rows, width=1000,title_text="振動数")
        fig.write_html(html)
    
    def write_mode(self,i=None):
        """虚振動モードのtrajファイルを作成する
        
        Parameters:
        
        i: int
            | TS{i}の虚振動の振動モードでtrajファイルを作成する
            | Noneを指定した場合,全てのTSの虚振動のtrajが作成される
            | 虚振動が確認できなかったTSのtrajファイルは作成されない
        forlder:str
            | 指定したフォルダ内にtrajファイルが作成される
            
        Note:
        
            事前にsummary()を実行しておく必要がある.
            
            >>> vib.summary()
            >>> vib.write_mode()
            
        Note:

            虚振動以外のモードを出力したい時,例えばTS5のmode=8の振動のtrajを作成する場合は
            
            >>> vib[5].write_mode(8)
            
        """
        if i is None:
            for i_mode,vib in zip(self.ivib,self.vib):
                if not i_mode is None:
                    vib.write_mode(i_mode)
        else:
            if not self.ivib[i] is None:
                self.vib[i].write_mode(self.ivib[i])    
        
    def __iter__(self):
        return iter(self.vib)
    
    def __len__(self):
        return len(self.vib)
    
    def __getitem__(self, item):
        return self.vib[item]