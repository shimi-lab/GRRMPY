import pandas as pd
import numpy as np
from ase.units import kJ,mol

from grrmpy.path.functions import to_excell,to_fig, to_html, to_plotly

class ReactPaths():
    """複数のReactPathをマージしてまとめたクラス
    """
    def __init__(self,data):
        """
        
        Parameters:
        
        data:
            ReactPathのリスト
        """
        if type(data)==list:
            self.react_path_list = data
        self.atoms_dict = {}
        for react_path in self.react_path_list:
            self.atoms_dict.update(react_path.get_atoms_dict())
            
    @property
    def data(self):
        data = {
            "name":self.get_atoms_dict().keys(),
            "energy":[atoms.get_potential_energy()*mol/kJ for atoms in self.get_atoms_dict().values()],
        }
        return data
            
    def get_atoms_dict(self):
        return self.atoms_dict
    
    def get_energy_dict(self,std_name=None):
        if std_name is None:
            std_name = self.react_path_list[0].get_name()[0]
        std_energy = self.atoms_dict[std_name].get_potential_energy()
        f = mol/kJ
        energy_dict = {name:(atoms.get_potential_energy()-std_energy)*f
                       for name,atoms in self.get_atoms_dict().items()}
        return energy_dict
    
    def get_title_list_and_solid_data(self,std_name=None):
        df_list = [react_path.get_solid_df() for react_path in self.react_path_list]
        marge_df = pd.concat(df_list)
        duplicated = marge_df[["solid_xini","solid_xfin","name"]].duplicated()
        marge_df = marge_df[~duplicated]
        marge_df = marge_df.reset_index(drop=True)
        energy_list = [self.get_energy_dict(std_name=std_name)[i] for i in marge_df["name"]]
        marge_df["solid_yini"] = energy_list
        marge_df["solid_yfin"] = energy_list
        
        titles = [react_path.title for react_path in self.react_path_list]
        title_list = []
        for df,title in zip(df_list,titles):
            title_list += [title for _ in range(len(df))]
    
        title_list = np.array(title_list)
        boolean = ~duplicated
        title_list = title_list[boolean.tolist()]
        return title_list, marge_df
        
    
    def get_solid_df(self,std_name=None):
        _, marge_df = self.get_title_list_and_solid_data(std_name=std_name)
        return marge_df
    
    def get_dot_df(self,std_name=None):
        df_list = [react_path.get_dot_df() for react_path in self.react_path_list]
        marge_df = pd.concat(df_list)
        marge_df = marge_df[~marge_df[["dot_xini","dot_xfin","dot_name"]].duplicated()]
        marge_df = marge_df.reset_index(drop=True)
        marge_df["dot_yini"] = [self.get_energy_dict(std_name=std_name)[dot_name.split("-")[0]] for dot_name in marge_df["dot_name"]]
        marge_df["dot_yfin"] = [self.get_energy_dict(std_name=std_name)[dot_name.split("-")[1]] for dot_name in marge_df["dot_name"]]
        return marge_df
    
    def write_html(self,outfile=None,std_name=None,**kwargs):
        title_list,solid_data = self.get_title_list_and_solid_data(std_name)
        dot_data = self.get_dot_df(std_name)
        if outfile is None:
            return to_plotly(solid_data,dot_data,title_list=title_list)
        html_text = to_html(solid_data,dot_data,**kwargs,title_list=title_list)
        with open(outfile,"w") as f:
            f.write(html_text)
            
    def write_excel(self,outfile:str,std_name=None,**kwargs):
        table_data = self.data
        solid_data = self.get_solid_df(std_name=std_name)
        dot_data = self.get_dot_df(std_name=std_name)
        to_excell(outfile,table_data,solid_data,dot_data,**kwargs)
        
    def preview(self,std_name=None):
        """Notebook上でグラフを表示(matplotlib)"""
        solid_data = self.get_solid_df(std_name=std_name)
        dot_data = self.get_dot_df(std_name=std_name)
        to_fig(solid_data,dot_data,xlabel="Energy [kJ/mol]")
        
    def get_fig(self,std_name=None):
        """matplotlibのfigを返す"""
        solid_data = self.get_solid_df(std_name=std_name)
        dot_data = self.get_dot_df(std_name=std_name)
        return to_fig(solid_data,dot_data)
                
    def __getitem__(self, item):
        if type(item) == int:
            return self.react_path_list[item]
        elif type(item) == str:
            title_list = np.array([react_path.title for react_path in self.react_path_list])
            react_path_list = np.array(self.react_path_list)
            mutch = title_list==item
            react_path_list = react_path_list[mutch].tolist()
            if len(react_path_list) == 1:
                return react_path_list[0]
            else:
                return react_path_list
            
    def __len__(self):
        return len(self.react_path_list)
    
    def __iter__(self):
        return iter(self.react_path_list)
    
    def __add__(self,other):
        from grrmpy.path.reaction_path import ReactPath
        if type(other) == ReactPaths:
            react_path_list = self.react_path_list + other.react_path_list
            return self.__class__(react_path_list)
        elif type(other) == ReactPath:
            react_path_list = self.react_path_list + [other]
            return self.__class__(react_path_list)
        
    def __iadd__(self,other):
        return self + other
        