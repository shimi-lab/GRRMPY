from ase import Atoms
from ase.io.trajectory import TrajectoryReader,SlicedTrajectory
from ase.geometry.analysis import Analysis
from ase.neighborlist import build_neighbor_list,natural_cutoffs
from jinja2 import Environment ,FileSystemLoader
from jinja2.ext import loopcontrols
import json 
import zipfile
import pickle
from pathlib import Path

from grrmpy.blender import default 

def make_py_script(file,Styles):
    """Belnder用のPythonスクリプトを作成する

    Parameters:
    
    file: str(.py or .zip)
        pythonファイル名, Animationが含まれる場合はzip名
        ファイル名がハイフン'-'の場合,標準出力する.(Animationがある場合は無効)
    Styles: BaseStyle object
        | BallAndStick,Stick,SpaceFilling,Animationのオブジェクト
        | 複数のstyleを組み合わせる場合,リストで与える.
    """
    if type(Styles) != list:
        Styles = [Styles]
    for style in Styles:
        into_one_file = True
        if style.style == "animation":
            into_one_file = False
            break
    p = Path(__file__).parent
    env = Environment(loader=FileSystemLoader(p/'template/', encoding='utf8'),extensions=['jinja2.ext.loopcontrols'])
    tmpl = env.get_template("template.py")
    if into_one_file:
        data_list = [style.todict() for style in Styles]
        data = {
            "data_list":data_list
        } 
        pyscript = tmpl.render(data)
        if file == "-":
            return pyscript
        else:
            with open(file,"w") as f:
                f.write(pyscript)
    else:
        data_list = []
        pkl_dict = {}
        for i,style in enumerate(Styles):
            d_dict = style.todict()
            if style.style == "animation":
                filename = f"positions{i}.pkl"
                d_dict["file"] = filename
                pkl_dict[filename] = style
            data_list.append(d_dict) 
        data = {
            "data_list":data_list,
        }
        pyscript = tmpl.render(data)
        write_position_zipfile(file,pyscript,pkl_dict)

def get_unique_bonds(atoms):
    cutoff = natural_cutoffs(atoms, mult=1)
    nl = build_neighbor_list(atoms,cutoff)
    ana = Analysis(atoms, nl=nl)
    bonds_list = ana.unique_bonds[0]
    bonds = []
    for idx_a,idx_list in enumerate(bonds_list):
        for idx_b in idx_list:
            bonds.append((idx_a,idx_b))
    return bonds

class BaseStyle():
    def __init__(self,atoms:Atoms,indices=None):
        """
        
        Parameters:
        
        atoms: Atoms or list of Atoms
            AtomsまたはAtomsのリスト
        indices: list of int
            一部の原子のみを表示する場合,index番号をリストで与える
        
        """
        self.atoms = atoms
        if type(atoms) == Atoms:
            if indices is None:
                indices = [i for i in range(len(atoms))]
            self.indices = indices
            atoms.set_pbc(False)
            self.chemical_symbols = atoms[self.indices].get_chemical_symbols()
            self.positions = atoms[self.indices].get_positions().tolist()
            self.bonds = get_unique_bonds(atoms[self.indices])
        else:
            if indices is None:
                indices = [i for i in range(len(atoms[0]))]
            self.indices = indices
            self.chemical_symbols = atoms[0][self.indices].get_chemical_symbols()
        self.unique_symbols = list(set(self.chemical_symbols)) 
        
    def set_param(self,permited_param,kwargs):
        for attr,default_value in permited_param.items():
            if type(default_value) == dict:
                val = kwargs.get(attr,{})
                val = dict(default_value, **val)
            else:
                val = kwargs.get(attr,default_value)
            setattr(self, attr, val)
            
    def get_parameters(self):
        return {attr:getattr(self, attr) for attr in self.permited_param}
            
    def todict(self,bonds=False):
        attr_list = ["bicolor","colors","radius","scale","sizes","stick_color","subdivision_surface","cartoon"]
        attr_list2 = ["style","chemical_symbols","unique_symbols","positions"]#,"bonds"]
        data_dict = {}
        for attr in attr_list:
            if hasattr(self, attr):
                data_dict[attr] = getattr(self, attr)
        for attr in attr_list2:
            if hasattr(self, attr):
                data_dict[attr] = getattr(self, attr)
        if bonds:
            data_dict["bonds"] = getattr(self, "bonds")
        return data_dict
            
    def write(self,file,bonds=False):
        data_dict = self.todict(bonds=bonds)
        with open(file, "w") as f:
            json.dump(data_dict, f)
        
class BallAndStick(BaseStyle):
    """Ball and Stickのスタイル
    
    bicolor=Trueにすると,Blender上での操作が重くなるので注意(オブジェクト数が多い)
    
    Parameters:
    
    atoms: Atoms
        Atomsオブジェクト
    indices: list of int
        一部の原子のみを表示する場合,index番号をリストで与える
    kwargs:
        bicolor: bool
            bicolorにする場合True.
        cartoon:dict
            apply: bool
                | cartoonを適用する場合Ture.
            IOR: float
                | フレネルのIOR(枠線の太さに相当する)
            color: tuple
                | ミックスのColor2(枠線の色に相当する)
                | 1で規格化されたRGBAで指定する
        colors : dict
            1で規格化したRGBA.
            ex) {'O':(1,0,0,1)}
        radius: float
            stickの半径
        scale: flaot
            | Ballの大きさ.(1の場合.SpaceFillingと同じサイズになる)
            | 元素毎に大きさを変更したい場合はscaleでなくsizesで指定する.
        sizes: dict
            | 元素毎のBallの大きさ.(共有結合半径)
            | {'H':0.46, 'C':0.77}}のように指定
        stick_color: tuple
            | 1で規格化されたRGBA.
            | bicolor=Falseの時のみ有効
        subdivision_surface: dict
            apply : bool
                | 適用する場合True.
            level : int
                | viewポートでのレベル
            render_levels: int
                | Renderレベル
    """
    style = "ball_and_stick"
    def __init__(self,atoms,indices=None,**kwargs):
        """Ball and Stickのスタイル
        
        bicolor=Trueにすると,Blender上での操作が重くなるので注意(オブジェクト数が多い)
        
        Parameters:
        
        atoms: Atoms
            Atomsオブジェクト
        indices: list of int
            一部の原子のみを表示する場合,index番号をリストで与える
        kwargs:
            bicolor: bool
                bicolorにする場合True.
            cartoon:dict
                apply: bool
                    cartoonを適用する場合Ture.
                IOR: float
                    フレネルのIOR(枠線の太さに相当する)
                color: tuple
                    | ミックスのColor2(枠線の色に相当する)
                    | 1で規格化されたRGBAで指定する
            colors : dict
                1で規格化したRGBA.
                ex) {'O':(1,0,0,1)}
            radius: float
                stickの半径
            scale: flaot
                | Ballの大きさ.(1の場合.SpaceFillingと同じサイズになる)
                | 元素毎に大きさを変更したい場合はscaleでなくsizesで指定する.
            sizes: dict
                | 元素毎のBallの大きさ.(共有結合半径)
                | {'H':0.46, 'C':0.77}}のように指定
            stick_color: tuple
                | 1で規格化されたRGBA.
                | bicolor=Falseの時のみ有効
            subdivision_surface: dict
                | - apply : bool
                | - level : int
                | - render_levels: int
        """
        super().__init__(atoms,indices)
        self.check_param()
        self.permited_param = {
            "bicolor":default.bicolor,
            "cartoon":default.cartoon,
            "colors":{symb:colors for symb,colors in default.color.items() if symb in self.unique_symbols},
            "radius":default.radius,
            "scale":default.scale,
            "sizes":{symb:size for symb,size in default.sizes.items() if symb in self.unique_symbols},
            "stick_color":default.bond_color,
            "subdivision_surface":default.subdivision_surface,
            }
        self.set_param(self.permited_param,kwargs)
        
    def check_param(self):
        if not type(self.atoms) == Atoms:
            raise TypeError(f"{self.__class__.__name__}のatomsはAtomsオブジェクトのみをサポートしています.")
        
    def write(self,file):
        super().write(file,bonds=True)
        
    def todict(self):
        return super().todict(bonds=True)
    
class Stick(BaseStyle):
    """Stickのスタイル
    
    bicolor=Trueにすると,Blender上での操作が重くなるので注意(オブジェクト数が多い)
    
    Parameters:
    
    atoms: Atoms
        Atomsオブジェクト
    indices: list of int
        一部の原子のみを表示する場合,index番号をリストで与える
    kwargs:
        bicolor: bool
            bicolorにする場合True.
        cartoon:dict
            apply: bool
                | cartoonを適用する場合Ture.
            IOR: float
                | フレネルのIOR(枠線の太さに相当する)
            color: tuple
                | ミックスのColor2(枠線の色に相当する)
                | 1で規格化されたRGBAで指定する
        colors : dict
            1で規格化したRGBA.
            ex) {'O':(1,0,0,1)} 
        radius: float
            stickの半径
        stick_color: tuple
            | 1で規格化されたRGBA.
            | bicolor=Falseの時のみ有効
        subdivision_surface: dict
            apply : bool
                | 適用する場合True.
            level : int
                | viewポートでのレベル
            render_levels: int
                | Renderレベル
    """
    style = "stick"
    def __init__(self,atoms,indices=None,**kwargs):
        """Stickのスタイル
        
        bicolor=Trueにすると,Blender上での操作が重くなるので注意(オブジェクト数が多い)
        
        Parameters:
        
        atoms: Atoms
            Atomsオブジェクト
        indices: list of int
            一部の原子のみを表示する場合,index番号をリストで与える
        kwargs:
            bicolor: bool
                bicolorにする場合True.
            cartoon:dict
                - apply: bool
                    cartoonを適用する場合Ture.
                - IOR: float
                    フレネルのIOR(枠線の太さに相当する)
                - color: tuple
                    | ミックスのColor2(枠線の色に相当する)
                    | 1で規格化されたRGBAで指定する
            colors : dict
                1で規格化したRGBA.
                ex) {'O':(1,0,0,1)} 
            radius: float
                stickの半径
            stick_color: tuple
                | 1で規格化されたRGBA.
                | bicolor=Falseの時のみ有効
        """
        super().__init__(atoms,indices)
        self.check_param()
        self.permited_param = {
            "bicolor":default.bicolor,
            "cartoon":default.cartoon,
            "colors":{symb:color for symb,color in default.color.items() if symb in self.unique_symbols},
            "radius":default.radius,
            "stick_color":default.bond_color,
            }
        self.set_param(self.permited_param,kwargs)
        
    def check_param(self):
        if not type(self.atoms) == Atoms:
            raise TypeError(f"{self.__class__.__name__}のatomsはAtomsオブジェクトのみをサポートしています.")
        
    def todict(self):
        return super().todict(bonds=True)
        
    def write(self,file):
        super().write(file,bonds=True)

class SpaceFilling(BaseStyle):
    """SpaceFillingのスタイル
    
    bicolor=Trueにすると,Blender上での操作が重くなるので注意(オブジェクト数が多い)
    
    Parameters:
    
    atoms: Atoms
        Atomsオブジェクト
    indices: list of int
        一部の原子のみを表示する場合,index番号をリストで与える
    kwargs:
        cartoon:dict
            apply: bool
                | cartoonを適用する場合Ture.
            IOR: float
                | フレネルのIOR(枠線の太さに相当する)
            color: tuple
                | ミックスのColor2(枠線の色に相当する)
                | 1で規格化されたRGBAで指定する
        colors : dict
            1で規格化したRGBA.
            ex) {'O':(1,0,0,1)}
        scale: flaot
            | Ballの大きさ.デフォルトは1.
            | 元素毎に大きさを変更したい場合はscaleでなくsizesで指定する.
        sizes: dict
            | 元素毎のBallの大きさ.(共有結合半径)
            | {'H':0.46, 'C':0.77}}のように指定
        subdivision_surface: dict
            apply : bool
                | 適用する場合True.
            level : int
                | viewポートでのレベル
            render_levels: int
                | Renderレベル
    """
    style = "space_filling"
    def __init__(self, atoms,indices=None,**kwargs):
        """SpaceFillingのスタイル
                
        Parameters:
        
        atoms: Atoms
            Atomsオブジェクト
        indices: list of int
            一部の原子のみを表示する場合,index番号をリストで与える
        kwargs:
            cartoon:dict
                - apply: bool
                    | cartoonを適用する場合Ture.
                - IOR: float
                    | フレネルのIOR(枠線の太さに相当する)
                - color: tuple
                    | ミックスのColor2(枠線の色に相当する)
                    | 1で規格化されたRGBAで指定する
            colors : dict
                1で規格化したRGBA.
                ex) {'O':(1,0,0,1)}
            scale: flaot
                | Ballの大きさ.デフォルトは1.
                | 元素毎に大きさを変更したい場合はscaleでなくsizesで指定する.
            sizes: dict
                | 元素毎のBallの大きさ.(共有結合半径)
                | {'H':0.46, 'C':0.77}}のように指定
            subdivision_surface: dict
                | - apply : bool
                |   | 適用する場合True.
                | - level : int
                |   | viewポートでのレベル
                | - render_levels: int
                |   Renderレベル
        """
        super().__init__(atoms,indices)
        self.check_param()
        self.permited_param = {
            "cartoon":default.cartoon,
            "colors":{symb:color for symb,color in default.color.items() if symb in self.unique_symbols},
            "scale":default.space_filling_scale,
            "sizes":{symb:size for symb,size in default.sizes.items() if symb in self.unique_symbols},
            "subdivision_surface":default.subdivision_surface,
            }
        self.set_param(self.permited_param,kwargs)
        
    def check_param(self):
        if not type(self.atoms) == Atoms:
            raise TypeError(f"{self.__class__.__name__}のatomsはAtomsオブジェクトのみをサポートしています.")
        
    def todict(self):
        return super().todict(bonds=False)
        
    def write(self,file):
        super().write(file,bonds=False)
    
class Animation(BaseStyle):
    """Animationのスタイル
    
    | 結合の描写が複雑なのでAnimationはSpaceFillingのみしかサポートしていない
    | Animationと他のスタイルを組み合わせることはできるがAnimationで指定した原子のみが動く.
    
    Parameters:
    
    images: Trajectory or list of Atoms
        TrajectoryまたはAtomsのリスト
    indices: list of int
        一部の原子のみを表示する場合,index番号をリストで与える
    kwargs:
        cartoon:dict
            apply: bool
                | cartoonを適用する場合Ture.
            IOR: float
                | フレネルのIOR(枠線の太さに相当する)
            color: tuple
                | ミックスのColor2(枠線の色に相当する)
                | 1で規格化されたRGBAで指定する
        colors : dict
            1で規格化したRGBA.
            ex) {'O':(1,0,0,1)}
        scale: flaot
            | Ballの大きさ.デフォルトは1.
            | 元素毎に大きさを変更したい場合はscaleでなくsizesで指定する.
        sizes: dict
            | 元素毎のBallの大きさ.(共有結合半径)
            | {'H':0.46, 'C':0.77}}のように指定
        start: int
            始めのキーフレームを打つ位置
        step: int
            何フレーム毎にキーを打つか
        subdivision_surface: dict
            apply : bool
                | 適用する場合True.
            level : int
                | viewポートでのレベル
            render_levels: int
                | Renderレベル
    """
    style = "animation"
    def __init__(self, images,indices=None,**kwargs):
        """Animationのスタイル
        
        | 結合の描写が複雑なのでAnimationはSpaceFillingのみしかサポートしていない
        | Animationと他のスタイルを組み合わせることはできるがAnimationで指定した原子のみが動く.
        | bicolor=Trueにすると,Blender上での操作が重くなるので注意(オブジェクト数が多い)
        
        Parameters:
        
        images: Trajectory or list of Atoms
            TrajectoryまたはAtomsのリスト
        indices: list of int
            一部の原子のみを表示する場合,index番号をリストで与える
        kwargs:
            cartoon:dict
                - apply: bool
                    cartoonを適用する場合Ture.
                - IOR: float
                    フレネルのIOR(枠線の太さに相当する)
                - color: tuple
                    | ミックスのColor2(枠線の色に相当する)
                    | 1で規格化されたRGBAで指定する
            colors : dict
                1で規格化したRGBA.
                ex) {'O':(1,0,0,1)}
            scale: flaot
                | Ballの大きさ.デフォルトは1.
                | 元素毎に大きさを変更したい場合はscaleでなくsizesで指定する.
            sizes: dict
                | 元素毎のBallの大きさ.(共有結合半径)
                | {'H':0.46, 'C':0.77}}のように指定
            start: int
                始めのキーフレームを打つ位置
            step: int
                何フレーム毎にキーを打つか
            subdivision_surface: dict
                | - apply : bool
                | - level : int
                | - render_levels: int
        """
        super().__init__(images,indices)
        self.check_param()
        self.permited_param = {
            "cartoon":default.cartoon,
            "colors":{symb:color for symb,color in default.color.items() if symb in self.unique_symbols},
            "scale":default.space_filling_scale,
            "sizes":{symb:size for symb,size in default.sizes.items() if symb in self.unique_symbols},
            "start":default.start,
            "step":default.step,
            "subdivision_surface":default.subdivision_surface
            }
        self.set_param(self.permited_param,kwargs)
        
    def check_param(self):
        if type(self.atoms) == TrajectoryReader or type(self.atoms) == SlicedTrajectory:
            return
        if type(self.atoms) == list:
            if type(self.atoms[0]) == Atoms:
                return
        raise TypeError(f"{self.__class__.__name__}のimagesはTrajectory(TrajectoryReader)またはAtomsのリストです.")
    
    def todict(self):
        # 親クラスを上書き
        attr_list = ["colors","scale","sizes","start","step","subdivision_surface","cartoon"]
        attr_list2 = ["style","chemical_symbols","unique_symbols"]
        data_dict = {}
        for attr in attr_list:
            data_dict[attr] = getattr(self, attr)
        for attr in attr_list2:
            data_dict[attr] = getattr(self, attr)
        return data_dict
    
    def write(self,file):
        super().write(file,bonds=False)
        
        
def write_position_zipfile(zipname,pyscript:str,data:dict):
    """zipファイルにpositions(Animation)を書きこむ
    dataは{"ファイル名(pkl)":Animation}の辞書
    """
    p = Path(zipname)
    if p.suffix != ".zip":
        raise Exception("fileの拡張子は.zipです")
    if p.exists():
        raise FileExistsError(f"{zipname}は既に存在します")
    with zipfile.ZipFile(zipname,"a") as zf:
        for file,animation in data.items():
            with zf.open(file,"w") as f:
                for atoms in animation.atoms:
                    pickle.dump(atoms[animation.indices].get_positions(),f)
        with zf.open(str(p.with_suffix(".py")),"w") as f:
            f.write(pyscript.encode())
        
                    
def write_position_zipfile_for_app(zipname,data:dict,interzip="position"):
    """zipファイルにpositions(Animation)を書きこむ
    dataは{"ファイル名(pkl)":Trajectory}の辞書
    """
    with zipfile.ZipFile(interzip,"a") as zf:
        for file,traj in data.items():
            with zf.open(file,"w") as f:
                for atoms in traj:
                    pickle.dump(atoms.get_positions(),f)
        btn = st.download_button(
            label="Download ZIP",
            data=zf,
            file_name="myfile.zip",
            mime="application/zip"
        )