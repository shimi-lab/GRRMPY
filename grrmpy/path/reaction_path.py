import pandas as pd
from ase.units import kJ, mol, Hartree
import pickle
from grrmpy.path.reaction_paths import ReactPaths

#user
from grrmpy.path.functions import to_excell,to_fig, to_html, to_plotly
from grrmpy.calculator import pfp_calculator

class ReactPath():
    """
    
        Parameters:
        
        data: dict or DataFrame
            energyとnameまたはatomsとnameの2つのキーを持つ辞書.

            - 'energy': list of float
                | エネルギーをリストで与える
                | 単位はkJ/mol, eV単位で与える場合,unit='eV'とする
            - 'atoms': list of Atoms
                | Atomsのリスト
            - 'name' : list of str
                | Ex) ['EQ0','TS1','EQ3','TS3','EQ2','PT2','EQ5']
                | 'TS','PT'の文字列を含んでいる場合,それぞれグラフ描写時に赤色,緑色で表示される
                
        title: srt
            | Pathを識別するためにタイトルを付ける事が可能.
            | デフォルトはNone
        positions: list of bool
            | solidラインの位置をカスタマイズする
            | 使用方法は   を参照
            | 指定しない場合はすべてTrueリストが作成される
        unit: string
            dataの設定の際に'energy'で設定した場合,与えたエネルギーの単位を記入する
            'kJ/mol','Hartree','eV'のいずれか
        clac_funct: function object
            | dataの設定の際に'atoms'で設定した場合,calculatorを与える関数が必要.
            | デフォルトではPFP.
        mask: list of bool
            dot線を書かない部分をFalseにする.Noneの場合全てTrueのリスト
        barrier_less_index: list of int
            | バリアレスの部分のindex番号(設定は任意)
            | EQ1-TS0-EQ2-EQ3でEQ2-EQ3の部分がバリアレスの場合,barrier_less_index=[2]
            
        Note:
            constraintsはFixAtomsだけfromdict,frompkl変換できる
    """
    def __init__(self,data:dict,positions:list=None,mask=None,barrier_less_index=[],title=None,unit="eV",calc_func=pfp_calculator):
        """

        Parameters:
        
        data: dict or DataFrame
            energyとnameまたはatomsとnameの2つのキーを持つ辞書.
            
            - 'energy': list of float
                | エネルギーをリストで与える
                | 単位はkJ/mol, eV単位で与える場合,unit='eV'とする
            - 'atoms': list of Atoms
                Atomsのリスト
            - 'name' : list of str
                | Ex) ['EQ0','TS1','EQ3','TS3','EQ2','PT2','EQ5']
                | 'TS','PT'の文字列を含んでいる場合,それぞれグラフ描写時に赤色,緑色で表示される
                
        title: srt
            | Pathを識別するためにタイトルを付ける事が可能.
            | デフォルトはNone
        positions: list of bool
            | solidラインの位置をカスタマイズする
            | 使用方法は   を参照
            | 指定しない場合はすべてTrueリストが作成される
        unit: string
            dataの設定の際に'energy'で設定した場合,与えたエネルギーの単位を記入する
            'kJ/mol','Hartree','eV'のいずれか
        clac_funct: function object
            | dataの設定の際に'atoms'で設定した場合,calculatorを与える関数が必要.
            | デフォルトではPFP.
        """
        self.data = data
        self.title = title
        self.calc_func = calc_func

        if positions is None:
            positions = [True for _ in range(len(self.data))]
        self._positions = positions
        
        if mask is None:
            mask = [True for _ in range(len(self.data)-1)]
        self.mask = mask
        
        self.barrier_less_index = barrier_less_index
        
        if 'atoms' in self.data.columns:
            for atoms in self.data["atoms"]:
                atoms.calc = self.calc_func()
            self.data["energy"] = [atoms.get_potential_energy()*mol/kJ for atoms in self.data["atoms"]]
        else:
            if unit == "eV":
                self.data["energy"] = [i*mol/kJ for i in self.data["energy"]]
            elif unit == "Hartree":
                self.data["energy"] = [i*mol/(kJ*Hartree) for i in self.data["energy"]]
        
    @property
    def data(self):
        return self._data
    
    def get_potential_energy_list(self):
        if 'atoms' in self.data.columns:
            return [atoms.get_potential_energy() for atoms in self.data["atoms"]]
        else:
            raise KeyError("dataにatomsがありません")
        
    def get_atoms_dict(self):
        if 'atoms' in self.data.columns:
            return {name:atoms for atoms,name in zip(self.data["atoms"],self.data["name"])}
        else:
            raise KeyError("dataにatomsがありません")
    
    @data.setter
    def data(self,data):
        if type(data) == dict:
            if not "name" in data.keys():
                raise KeyError("'name'キーがありません")
            self._data = pd.DataFrame(data)
        elif type(data) == pd.DataFrame:
            if not "name" in data.columns:
                raise KeyError("'name'キーがありません")
            self._data = data
        self._positions = [True for _ in range(len(self.data))]
        self._data = self._data.reset_index(drop=True)
    
    def get_energy(self):
        return self.data["energy"].to_list()
    
    def get_name(self):
        return self.data["name"].to_list()
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self, positions):
        n_true = len([i for i in positions if i]) # Trueの個数
        if not n_true == len(self.data):
            raise ValueError("Trueの数がデータの数と一致しません")
        self._positions = positions
        
    def get_solid_df(self):
        solid_xini = [i*2+1 for i,b in enumerate(self.positions) if b]
        solid_xfin = [i+1 for i in solid_xini]
        std_e = self.get_energy()[0]
        solid_yini = [i-std_e for i in self.get_energy()]
        df = pd.DataFrame({"solid_xini":solid_xini,
                           "solid_xfin":solid_xfin,
                           "name":self.get_name(),
                           "solid_yini":solid_yini,
                           "solid_yfin":solid_yini
                           })
        return df
    
    def get_dot_df(self,mask=None):
        solid_df = self.get_solid_df()
        dot_xini = solid_df["solid_xfin"][:-1].to_list()
        dot_xfin = solid_df["solid_xini"][1:].to_list()
        dot_yini = solid_df["solid_yini"][:-1].to_list()
        dot_yfin = solid_df["solid_yini"][1:].to_list()
        dot_name = [f"{self.data['name'][i]}-{self.data['name'][i+1]}" 
                    for i in range(len(self.data)-1)]
        if mask is None:
            mask = self.mask
        dot_xini = [i for m,i in zip(mask,dot_xini) if m]
        dot_xfin = [i for m,i in zip(mask,dot_xfin) if m]
        dot_yini = [i for m,i in zip(mask,dot_yini) if m]
        dot_yfin = [i for m,i in zip(mask,dot_yfin) if m]
        dot_name = [i for m,i in zip(mask,dot_name) if m]
        df = pd.DataFrame({"dot_xini":dot_xini,
                           "dot_xfin":dot_xfin,
                           "dot_name":dot_name,
                           "dot_yini":dot_yini,
                           "dot_yfin":dot_yfin})
        return df
    
    def get_ea(self):
        """律速段階の活性化障壁と反応を出力"""
        ts = self.data['name'].str.contains('TS')
        pt = self.data['name'].str.contains('PT')
        i_ts_pt = [i for i,(t,p) in enumerate(zip(ts,pt)) if any([t,p])] # TS,PTの位置
        e = [self.data['energy'][i]-self.data['energy'][i-1] for i in i_ts_pt]
        ea = max(e)
        ea_idx = i_ts_pt[e.index(ea)]
        reac_name = [self.data['name'][ea_idx-1],self.data['name'][ea_idx],self.data['name'][ea_idx+1]]
        return max(e),reac_name
          
    def write_excel(self,outfile:str,mask=None,**kwargs):
        """Excellにグラフを作成する
        
        その他の引数はreaction_path.functions.to_excell()を参照

        Parameters:
        
        outfile: str
            | Excelファイル名

        """
        table_data = self.data
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df(mask)
        to_excell(outfile,table_data,solid_data,dot_data,**kwargs)
        
    def write_html(self,outfile=None,mask=None,annotation=[],barrier_less_index=None,**kwargs):
        """反応座標グラフをhtmlに保存する

        Parameters:
        
        outfile: str or Path
            | 出力先のhtmlファイル
            | Noneの場合は標準出力する
            
        Examples:
        
            通常
            
            >>> path.write_html("Sample.html")
            
            アノテーション(x座標が14.5の部分に"Sample text"と入力)
            
            >>> path.write_html("Sample.html",annotation=[(14.5,"Sample text")])
        """
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df(mask)
        if barrier_less_index is None:
            barrier_less_index = self.barrier_less_index
        annotation += [(2*i+2.5,"barrier less") for i in self.barrier_less_index]
        if outfile is None:
            return to_plotly(solid_data,dot_data,annotation=annotation)
        html_text = to_html(solid_data,dot_data,annotation=annotation,**kwargs)
        with open(outfile,"w") as f:
            f.write(html_text)
        
    def preview(self,mask=None):
        """Notebook上でグラフを表示(matplotlib)"""
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df(mask)
        to_fig(solid_data,dot_data,xlabel="Energy [kJ/mol]")
        
    def get_fig(self,mask=None):
        """matplotlibのfigを返す"""
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df(mask)
        return to_fig(solid_data,dot_data)
    
    def __str__(self):
        return f"{self.title}: {'->'.join(self.get_name())}"
    
    def topkl(self,outfile):
        """現在のオブジェクトをpickle化する
        
        >>> path.topkl("Path.pickle")
        """
        # with open(outfile, "wb") as f:
        #     pickle.dump(self.todict(),f)
        if 'atoms' in self.data.columns:
            for atoms in self.data["atoms"]:
                del atoms.calc
        with open(outfile, "wb") as f:
            pickle.dump(self,f)
            
            
    @classmethod
    def frompkl(cls,openfile,calc_func=pfp_calculator):
        """topklで作成したpickleファイルからオブジェクト再生する
        
        >>> path = ReactPath.frompkl("Path.pickle")
        """
        # with open(openfile, "rb") as f:
        #     data_dict = pickle.load(f)
        # return self.fromdict(data_dict)
        with open(openfile, "rb") as f:
            react_path = pickle.load(f)
        if 'atoms' in react_path.data.columns:
            for atoms in react_path.data["atoms"]:
                atoms.calc = calc_func()
        return react_path
    
    @classmethod
    def from_reactionstring(cls,obj):
        """Matlantis FeatureのReactionStringFeatureの結果(ReactionStringFeatureResults)からReactPathオブジェクトを作成する"""
        return ReactPath._from_reactionstring(obj,offset=0)
    
    @classmethod
    def _from_reactionstring(cls,obj,offset):
        ts_idxes = obj.transition_state_indices
        atoms_list = []
        name_list = []
        for i,(reaction,ts_idx) in enumerate(zip(reaction_string_images,ts_idxes),start=offset):
            atoms_list.append(reaction[0].ase_atoms.copy())
            name_list.append(f"EQ{i}")
            atoms_list.append(reaction[ts_idx].ase_atoms.copy())
            name_list.append(f"TS{i}")
            fin_name_idx = i
        atoms_list.append(reaction_string_images[-1][-1].ase_atoms.copy())
        name_list.append(f"EQ{fin_name_idx}")
        return self.__class__(data={"atoms":atoms_list,"name":name_list})
    
    def __add__(self,other):
        return ReactPaths([self,other])
    
    def __iadd__(self,other):
        return self + other