from pathlib import PurePath
import networkx as nx
import itertools
import pandas as pd
from matplotlib.colors import rgb2hex
import matplotlib.pyplot as plt
from pyvis.network import Network
from ase.units import kJ,mol,Hartree,eV

#USER
from grrmpy.io.read_listlog import log2atoms,read_connections,read_energies
from grrmpy.calculator import pfp_calculator
from grrmpy.conv.atoms2smiles import atomslist2smileses

def make_color(x,cm="gnuplot"):
    """カラースケール

    Parameters:
    
    x :
        0~1の間の数字
    cm str:
        'jet','viridis','plasma','inferno','magma','cividis','gnuplot','CMRmap','rainbow'
    """
    cmap = plt.get_cmap(cm)
    colorcode = rgb2hex(cmap(x))
    return colorcode

def make_color_scale(val_list,cm="gnuplot"):
    max_val = max(val_list)
    min_val = min(val_list)
    val_scale = [(i-min_val)/(max_val-min_val) for i in val_list]
    color_scale = [make_color(x,cm=cm) for x in val_scale]
    return color_scale

class NodeData():
    def __init__(self,eq_list:list,energies:list,cm="gnuplot",indices=None):
        self.energies = energies
        self.atoms_list = eq_list
        self.all_data = pd.DataFrame({
            "node":[i for i in range(len(self))],
            "label":[f"{i}" for i in range(len(self))],
            "title":[f"EQ{i}\nEnergy:{e:.03f} kJ/mol" for i,e in zip(range(len(self)),self.energies)],
            "color":make_color_scale(self.energies,cm=cm),
            "energy":self.energies,
            "group":self.grouping(indices),
            })
    
    @property
    def df(self):
        return self.all_data[["node","label","title","color"]]
    
    @property
    def data(self):
        """xlwingで動かすプログラムと互換性を保つため"""
        data = {
            "node":[i for i in range(len(self))],
            "name":[f"EQ{i}" for i in range(len(self))],
            "group":self.all_data["group"],
            "E/Hartree":[i*(kJ/mol)/Hartree for i in self.all_data["energy"]],
            "E/kJmol-1":self.all_data["energy"]
        }
        return pd.DataFrame(data)
    
    @property
    def nx_data(self):
        columns = self.df.columns
        return [(node,dict(zip(columns[1:],data))) for i,node,*data in self.df.itertuples()]
    
    def change_color_scale(self,cm):
        self.all_data["color"] = make_color_scale(self.energies,cm)
        
    def grouping(self,indices=None):
        smiles_list = atomslist2smileses(self.atoms_list,target0=indices)
        unique = sorted(set(smiles_list), key=smiles_list.index)
        return [unique.index(smiles) for smiles in smiles_list]
        
    def __len__(self):
        return len(self.atoms_list)
    
class EdgeData():
    def __init__(self,ts_list,connections,energies:list):
        self.energies = energies
        self.atoms_list = ts_list
        self.connections = connections
        self.all_data = pd.DataFrame({
            "edge":[i for i in range(len(self))],
            "label":[f"{i}" for i in range(len(self))],
            "title":[f"TS{i}\n"+
                     f"Energy:{e:.03f} kJ/mol\n"+
                     f"connections:{source}-{target}"
                     for i,e,(source,target)
                     in zip(range(len(self)),self.energies,self.connections)],
            "energy":self.energies,
            "source":self.source,
            "target":self.target,
            })
        
    @property
    def source(self):
        return [source for source,_ in self.connections]
        
    @property
    def target(self):
        return [target for _,target in self.connections]
    
    @property
    def data(self):
        """xlwingで動かすプログラムと互換性を保つため"""
        data = {
            "edge":[i for i in range(len(self))],
            "name":[f"TS{i}" for i in range(len(self))],
            "source":self.source,
            "target":self.target,
            "E/Hartree":[i*(kJ/mol)/Hartree for i in self.all_data["energy"]],
            "E/kJmol-1":self.all_data["energy"]
        }
        return pd.DataFrame(data)
    
    def nx_data(self,self_roop=True):
        """エッジのタプルのリストを返す(??などは消す)"""
        if self_roop:
            edge_list = [(ini_idx,fin_idx) for ini_idx,fin_idx in self.connections 
                        if isinstance(ini_idx, int) and isinstance(fin_idx, int)]
        else:
            edge_list = [(ini_idx,fin_idx) for ini_idx,fin_idx in self.connections 
                        if isinstance(ini_idx, int) and isinstance(fin_idx, int) and ini_idx!=fin_idx]
        return edge_list
    
    def __len__(self):
        return len(self.atoms_list)
    
class NetGraph():
    """
    Parameters:
    
    \*listlog: 
        | 方法1: EQ_list.log,TS_list.log,,PT_list.logを与える.
        | 方法2: eq_list,ts_list,connectionsを与える.
    indices: list of int
        | 原子のindex番号のリスト.
        | 同一構造をグループ化するために用いる.
        | 指定した原子同士の結合状態から同一構造かを判定する.
    comfile:
        comfileパス
    poscar:
        POSCARパス
    calc_func:
        calculatorを返す関数
            
    Examples:

        >>> G = NetGraph("XXX_EQ_list","XXX_TS_list","XXX_PT_list",comfile="XXX.com",poscar="POSCAR",indices=[i for i in range(0:10)])
        or
        >>> G = NetGraph("XXX_EQ_list","XXX_TS_list",comfile="XXX.com",poscar="POSCAR",indices=[i for i in range(0:10)])
        or
        >>> G = NetGraph(eq_list,ts_list,connections,indices=[i for i in range(0:10)])

    """
    def __init__(self,
                 *listlog,
                 indices = None,
                 comfile = None,
                 poscar = None,
                 calc_func=pfp_calculator,
                 ):
        """

        Parameters:
        
        *listlog: 
            | 方法1: EQ_list.log,TS_list.log,,PT_list.logを与える.
            | 方法2: eq_list,ts_list,connectionsを与える.
        indices: list of int
            | 原子のindex番号のリスト.
            | 同一構造をグループ化するために用いる.
            | 指定した原子同士の結合状態から同一構造かを判定する.
        comfile:
            comfileパス
        poscar:
            POSCARパス
        calc_func:
            calculatorを返す関数
        """
        # Node
        if isinstance(listlog[0],str) or isinstance(listlog[0],PurePath):
            # EQ_list.log,TS_list.log,PT_list.logで与えた場合
            # Node
            eq_energies = read_energies(listlog[0])
            eq_energies = [e*Hartree*mol/kJ for e in eq_energies]
            node_atoms = log2atoms(listlog[0],comfile,poscar)
            # EDGE
            ts_energies = list(itertools.chain.from_iterable(
                [read_energies(log) for log in listlog[1:]]
                ))
            ts_energies = [e*Hartree*mol/kJ for e in ts_energies]
            edge_atoms = list(itertools.chain.from_iterable(
                [log2atoms(log,comfile,poscar) for log in listlog[1:]]
                ))
            connections = list(itertools.chain.from_iterable(
                [read_connections(log) for log in listlog[1:]]
                ))
        else:
            # eq_list,ts_list,connectionsで与えた場合
            # NODE
            node_atoms = listlog[0]
            eq_energies = self._set_calc_and_get_energy(node_atoms,calc_func)
            eq_energies = [e*mol/kJ for e in eq_energies]
            # EDGE
            edge_atoms = listlog[1]        
            ts_energies = self._set_calc_and_get_energy(edge_atoms,calc_func)
            ts_energies = [e*mol/kJ for e in ts_energies]
            connections = listlog[2]
        self.node = NodeData(node_atoms,eq_energies,indices=indices)
        self.edge = EdgeData(edge_atoms,connections,ts_energies)
            
    def _set_calc_and_get_energy(self,atoms_list,calc_func):
        try:
            energies = [atoms.get_potential_energy() for atoms in atoms_list]
        except Exception:
            for atoms in atoms_list:
                atoms.calc = calc_func()
            energies = [atoms.get_potential_energy() for atoms in atoms_list]
        return energies
    
    @property
    def forward_energy(self):
        """"単位はkJ/mol,もしマイナスになる場合0とする"""
        node_e = self.node.energies
        edge_e = self.edge.energies
        source_idx = self.edge.source
        return [None if isinstance(i,str) else e-node_e[i] if e-node_e[i]>=0 else 0 
                for e,i in zip(edge_e,source_idx)]
    
    @property
    def reverse_energy(self):
        """"単位はkJ/mol"""
        node_e = self.node.energies
        edge_e = self.edge.energies
        target_idx = self.edge.target
        return [None if isinstance(i,str) else e-node_e[i] if e-node_e[i]>=0 else 0 
                for e,i in zip(edge_e,target_idx)]
    
    @property
    def node_df(self):
        """xlwingで動かすプログラムと互換性を保つため"""
        return self.node.data
    
    @property
    def edge_df(self):
        """xlwingで動かすプログラムと互換性を保つため"""
        df = self.edge.data
        df["forward/kJmol-1"] = self.forward_energy
        df["reverse/kJmol-1"] = self.reverse_energy
        return df
    
    def get_graph(self,self_loop=False,cm="gnuplot",node=None,edge=None):
        """NetWorkXのグラフを返す

        Parameters:
        
        self_loop: bool
            | 自己ループを含める場合True.
        cm : str
            | "jet","hot"などmatplotlibで指定できるカラースケール.

        Returns:
            networkx.Graph: networkxのグラフ
        """
        if node is None:
            node = self.node.df.to_dict(orient='list')["node"]
        if edge is None:
            edge = self.edge.nx_data(self_loop)
        G = nx.Graph()
        self.node.change_color_scale(cm)
        G.add_nodes_from(node)
        G.add_edges_from(edge)
        return G
    
    def write_html(self,html,height="500px",width="100%",cm="gnuplot",notebook=True, show_buttons=False,self_loop=False):
        """htmlにグラフを作成する

        Parameters:
        
        html: str
            保存名.html
        height: str
            縦幅
        width: str
            横幅
        cm: str
            | "jet","hot"などmatplotlibで指定できるカラースケール.
            | デフォルトは"gnuplot".
        notebook: bool
            | Trueした場合,Notebook上に表示される.
            | (現在機能していない)
        show_buttons: bool
            Trueにした場合,全ての機能を搭載したボタンが作成される.
        self_loop: bool
            自己ループを表示する場合はTrue.
        """
        g = Network(height, width,notebook=notebook)
        g.from_nx(self.get_graph(self_loop,cm))
        if show_buttons:
            g.show_buttons(True)   # 全機能使用
        else:
            g.show_buttons(filter_=['physics', 'nodes']) # 一部の機能のみ使用
        g.show(html)
        
    @property
    def _weight_edge(self):
        """最短距離を計算する時に用いる"""
        weight_edge = [[] for _ in range(len(self.node))]
        df = self.edge_df
        s = df.apply( lambda r: set((r['source'], r['target'])), axis=1)
        df = df.loc[s.drop_duplicates().index]
        for i,edge,name,source,target,_,_,forward_e,reverse_e,in df.itertuples():
            if isinstance(source, int) and isinstance(target, int) and not source==target:
                weight_edge[source].append([target,forward_e,name])
                weight_edge[target].append([source,reverse_e,name])
        return weight_edge
    
    @property
    def _group(self):
        """最短距離を計算する時に用いる"""
        df = self.node_df
        group_dict = {str(group_name):df[df["group"]==group_name]["node"].astype('str').tolist() for group_name in df["group"].unique()} 
        return group_dict
    
    @property
    def _node_data_for_graphml(self):
        df = self.node_df
        df = df[df.columns[df.columns != 'E/Hartree']]
        df["color"] = self.node.df["color"]
        key = ["id","name","group","E/kJmol-1","color"]
        node_data = [(node,dict(zip(key,[node,name,str(group)]+vals))) for _,node,name,group,*vals in df.itertuples()]
        return node_data
        
    @property
    def _edge_data_for_graphml(self):
        df = self.edge_df
        df = df[df.columns[df.columns != 'E/Hartree']]
        key = ["edge","name","E/kJmol-1","forward/kJmol-1","reverse/kJmol-1"]
        edge_data = [(source,target,dict(zip(key,[edge,name]+vals))) 
                     for _,edge,name,source,target,*vals in df.itertuples()
                     if not isinstance(source,str) and not isinstance(target,str) and source!=target]
        return edge_data
        
    def write_graphml(self,graphml):
        """
        graphmlファイルを作成する

        Parameters:
        
        graphml: str or Path
            | graphmlパス(.graphml)
        """
        node = self._node_data_for_graphml
        edge = self._edge_data_for_graphml
        nx.write_graphml(self.get_graph(node=node,edge=edge),path)