from ase.optimize import FIRE
from ase.io import write
from ase.neb import interpolate
from ase.units import kJ, mol
from pathlib import Path
import re
from configparser import ConfigParser, ExtendedInterpolation
from tqdm.notebook import tqdm_notebook as tqdm
# USER
from grrmpy.calculator import pfp_calculator
from grrmpy.neb.functions import find_local_minimum_points_and_ea,get_appropriate_nimages,get_imax
from grrmpy.functions import partially_overlapping_groupings
from grrmpy.neb.repetitive_neb import RepetitiveNEB
from grrmpy.opt.repetitive_opt import RepetitiveOpt
from grrmpy.io.html import write_html
from grrmpy.geometry.functions import get_bond_diff
from grrmpy.path.create_reactpath import create_reactpath
from grrmpy.automate.queue import Queue
from grrmpy.io.read_traj import read_traj
try:
    from grrmpy.optimize import FIRELBFGS
    DEFAULT_OPTIMIZER = FIRELBFGS
except:
    from ase.optimize import LBFGS
    DEFAULT_OPTIMIZER = LBFGS

def is_similar(atoms1,atoms2,energy_dif=10,bond_dif=(1,1.0),mult=1,**kwargs):
    """atoms1とatoms2を同一構造とみなすか(NEBでつなぐ必要はないか?)
    
    PinNEBの引数:isomorphism_funcのための関数
    
    基本的には結合状態が同じ且つ,エネルギー差がenergy_dif(kJ/mol)以下であれば同一構造とみなす.
    また,bond_dif[0] 本以下の結合状態の違いであり,その結合の差がbond_dif[1] Å以下の時も同一構造とする.
    
    Parameters:
    
    atoms1: Atoms
        Atomsオブジェクト
    atoms2: Atoms:
        Atomsオブジェクト
    energy_dif: float
        同じ結合状態かつ,エネルギー差がenergy_dif(kJ/mol)の時,同一構造とする
    bond_dif: Tuple[int,flaot]
        | 1要素目は許容できる結合の差(結合本数)
        | 2要素目は結合状態の異なる部分の許容できる結合距離の差(Å)
    """
    energy1 = get_energy(atoms1)
    energy2 = get_energy(atoms2)
    if abs(energy1-energy2)*mol/kJ > energy_dif:
        # energy_dif以下のエネルギーであるか
        return False
    formation, cleavage = get_bond_diff(atoms1,atoms2,mult=mult,**kwargs)
    if len((dif_idx:=formation+cleavage))>bond_dif[0]:
        # bond_dif[0]以下の結合状態の違いでしかないか
        return False
    for a_idx,b_idx in dif_idx:
        d1 = atoms1.get_distance(a_idx,b_idx)
        d2 = atoms2.get_distance(a_idx,b_idx)
        if abs(d1-d2)>bond_dif[1]:
            return False
    return True

def get_energy(atoms,calc_func=pfp_calculator)->float:
    try:
        return atoms.get_potential_energy()
    except RuntimeError:
        atoms.calc = calc_func()
        return atoms.get_potential_energy()
        
class AutoPinNEB():
    def __init__(self,
                 eq_list,
                 connections,
                 constraints=[], # connectionsには??などがあってはならない
                 parallel=True,
                 mic=None,
                 priority=0,
                 maxjob = 20,
                 calc_func=pfp_calculator,
                 eq_folder="EQ_List",
                 ts_folder="TS_List",
                 images_folder="Images",
                 unconv_folder="unconv",
                 logfolder = "log",
                 error_file="ERROR",
                 path_file="AllPath.dat"):
        self.eq_list = eq_list
        self.connections = connections
        self.constraints = constraints
        self.parallel = parallel
        self.priority = priority
        self.maxjob = maxjob
        self.calc_func = calc_func
        self.eq_folder = eq_folder
        self.ts_folder = ts_folder
        self.images_folder = images_folder
        self.unconv_folder = unconv_folder
        self.logfolder = logfolder
        for folder in [self.eq_folder,self.ts_folder,self.images_folder,self.unconv_folder,self.logfolder]:
            self.mkdir(folder)
        self.error_file = error_file
        self.path_file = path_file
        self.eqs = [] # 最適化後のEQのリスト
        self.tss = [] # 最適化後のTSのリスト
        self.path_txt = [] # List[List[str]]
        self.path_title = []
        if mic is None:
            mic = True if any(self.eq_list[0].get_pbc()) else False
        self.mic = mic
        self.param = self.default_params
        self.set_param(**self.param)
        self.preneb: RepetitiveNEB
        self.neb: RepetitiveNEB
        self.opt: DEFAULT_OPTIMIZER
        
    @property
    def default_params(self):
        return {
        "separate_completely":True,
        "minenergy":10,
        "barrier_less":0,
        "nimages":(0.2,32,12),
        "isomorphism_func":is_similar,
        "rename_preneb_file":True,
        "preneb_optimizer":FIRE,
        "preneb_climb_list":[False,False,False,False,False],
        "preneb_maxstep_list":[0.02,0.05,0.1,0.2,0.25],
        "preneb_automaxstep_climb_false":{1:0.05, 0.5:0.1, 0.1:0.2, 0.05:0.3, 0:0.35},
        "preneb_automaxstep_climb_true":{0.1:0.03,0:0.01},
        "preneb_fmax_list":[30,5,1,0.3,0.1],
        "preneb_steps_list":[500,500,500,500,500],
        "preneb_kwargs":{},
        "prenebopt_kwargs":{},
        "neb_optimizer":FIRE,
        "neb_climb_list":[False,False,False,False,True,True],
        "neb_maxstep_list":[0.02,0.05,0.1,None,0.01,None],
        "neb_automaxstep_climb_false":{1:0.05, 0.5:0.1, 0.1:0.2, 0.05:0.3, 0:0.35},
        "neb_automaxstep_climb_true":{0.1:0.03,0:0.01},
        "neb_fmax_list":[0.01,0.05,0.05,0.05,0.05,0.05],
        "neb_steps_list":[100,200,200,2000,500,1500],
        "neb_kwargs":{},
        "nebopt_kwargs":{},
        "opt_optimizer":DEFAULT_OPTIMIZER,
        "opt_maxstep_list":[0.02,0.05,0.1,0.2],
        "opt_automaxstep":{1:0.1, 0.5:0.2, 0.1:0.3, 0.05:0.4, 0:0.5},
        "opt_fmax_list":[0.001,0.001,0.001,0.001],
        "opt_steps_list":[50,100,300,10000],
        "opt_kwargs":{},
        "same_energy_and_fmax_arg":{"times":10,"row":2,"force_digit":4,"energy_digit":6},
        "steep_peak_arg":{"max_ea":800, "nsteps":100},
        }
        
    @property
    def preneb_param(self):
        attrs = ["preneb_optimizer", "preneb_climb_list", "preneb_maxstep_list",
                "preneb_automaxstep_climb_false", "preneb_automaxstep_climb_true",
                "preneb_fmax_list", "preneb_steps_list", "preneb_kwargs",
                "prenebopt_kwargs","same_energy_and_fmax_arg","steep_peak_arg"]
        keys = ["base_optimizer", "climb_list", "maxstep_list",
                "automate_maxstep_false", "automate_maxstep_true",
                "fmax_list", "steps_list", "neb_kwargs",
                "opt_kwargs","same_energy_and_fmax_arg","steep_peak_arg"] 
        return {key:self.param[attr] for attr,key in zip(attrs,keys)}
    
    @property
    def neb_param(self):
        attrs = ["neb_optimizer", "neb_climb_list", "neb_maxstep_list",
                "neb_automaxstep_climb_false", "neb_automaxstep_climb_true",
                "neb_fmax_list", "neb_steps_list", "neb_kwargs",
                "nebopt_kwargs","same_energy_and_fmax_arg","steep_peak_arg"]
        keys = ["base_optimizer", "climb_list", "maxstep_list",
                "automate_maxstep_false", "automate_maxstep_true",
                "fmax_list", "steps_list", "neb_kwargs",
                "opt_kwargs","same_energy_and_fmax_arg","steep_peak_arg"] 
        return {key:self.param[attr] for key,attr in zip(keys,attrs)}
        
    @property
    def opt_param(self):
        attrs = ["opt_optimizer", "opt_maxstep_list", "opt_automaxstep","opt_fmax_list",
                 "opt_steps_list", "opt_kwargs"]
        keys = ["optimizer","maxstep_list","automate_maxstep","fmax_list","steps_list",
                "opt_kwargs"]
        return {key:self.param[attr] for key,attr in zip(keys,attrs)}
    
    def set_param(self,**param):
        default_params = self.param
        unexpected_argument = set(param)-set(default_params)
        if len(unexpected_argument) != 0:
            raise Exception(f"予期しない引数({', '.join(unexpected_argument)})が存在します.\n"+
                            f"'self.default_params??'で引数を確認できます\n"
                           f"Parameters:\n{', '.join(self.default_params.keys())}")
        self.param.update(param)
        for key,val in default_params.items():
            setattr(self, key, val)
            
    def get_param(self):
        """
        
        Parameters:
        
        separate_completely: bool
            | 通常のNEBで収束した時,エネルギーダイアグラム上に極小点があれば,極小点を最適化してNEB分割する.
        minenergy: float
            | minenergy(kJ/mol)以下の活性化エネルギーの場合は極小点があっても無視する.
        barrier_less: float
            | barrier_less(kJ/mol)以下の活性化エネルギーの場合,バリアレスとみなす．
        nimages: int or tuple
            | intので指定した場合,指定した数のイメージ数で新たなNEBイメージを作成する.
            | tupleで指定する場合,3要素(a,b,c)のタプルを与える.
            | aÅ毎にイメージを作成する.最大でもb個,最低でもc個のイメージを作成する.
        isomorphism_func: function object
            | 同一構造判定を行なう関数,デフォルトではis_similarを用いる.
        rename_preneb_file: bool
            | Trueの場合,pre_NEBファイルは消す
        preneb_optimizer: Optimizer
            PreNEBのOptimizer
        preneb_climb_list: list of bool
            | climbのリスト
        preneb_maxstep_list: list of float
            | maxstepのリスト
        preneb_automaxstep_climb_false: 
            | {1:0.05, 0.5:0.1, 0.1:0.2, 0.05:0.3, 0:0.35}
        preneb_automaxstep_climb_true: dict
            {0.1:0.03,0:0.01}
        preneb_fmax_list :list of float
            | 
        preneb_steps_list: list of int
            | 
        preneb_kwargs:
            | 
        prenebopt_kwargs:
            | 
        neb_optimizer: Optimizer
            | 
        "neb_climb_list": list of bool
            | 
        neb_maxstep_list: 
            | [0.02,0.05,0.1,None,0.01,None],
        neb_automaxstep_climb_false:
            | {1:0.05, 0.5:0.1, 0.1:0.2, 0.05:0.3, 0:0.35},
        neb_automaxstep_climb_true":
            | {0.1:0.03,0:0.01},
        neb_fmax_list:
            | [0.01,0.05,0.05,0.05,0.05,0.05],
        neb_steps_list:
            | [100,200,200,2000,500,1500],
        neb_kwargs:
            | {}
        nebopt_kwargs:
            | {}
        opt_optimizer:
            | DEFAULT_OPTIMIZER
        opt_maxstep_list:
            | [0.02,0.05,0.1,0.2]
        opt_automaxstep:
            | {1:0.1, 0.5:0.2, 0.1:0.3, 0.05:0.4, 0:0.5}
        opt_fmax_list:
            | [0.001,0.001,0.001,0.001]
        opt_steps_list:
            | [50,100,300,10000]
        opt_kwargs:
            | {}
        same_energy_and_fmax_arg:
            | {"times":10,"row":2,"force_digit":4,"energy_digit":6},
        steep_peak_arg:
            | {"max_ea":800, "nsteps":100}
        
        """
        return self.param
        
    def mkdir(self,folder):
        p = Path(folder)
        if not p.exists():
            p.mkdir()
            
    def rename(self,old_name,new_name,extention_list=["traj","html","log"]):
            for extention in extention_list:# prenebのReName
                Path(f"{old_name}.{extention}").rename(f"{new_name}.{extention}")
                
    def register_eqs(self,atoms):
        self.eqs.append(atoms)
        eq_num = len(self.eqs)-1
        write(f"{self.eq_folder}/{eq_num}.traj",atoms)
        return eq_num
        
    def register_tss(self,images,name,i,ini_idx,fin_idx): # ts_idx GRRMのTS番号
        imax = get_imax(images,self.barrier_less,self.calc_func)
        if imax is not None: # バリアレスでない場合
            ts = images[imax]
            self.tss.append(ts)
            write(f"{self.ts_folder}/{len(self.tss)-1}.traj",ts)
            new_path_txt = self.make_path_txt(self.path_txt[i],ini_idx,fin_idx,len(self.tss)-1)
        else: # バリアレスの場合
            new_path_txt = self.make_path_txt(self.path_txt[i],ini_idx,fin_idx,None)
        self.path_txt[i] = new_path_txt
        self.write_path(i)
        write(f"{self.images_folder}/{name}.traj",images)            

    def get_energy(self,atoms):
        try:
            return atoms.get_potential_energy()
        except:
            atoms.calc = self.calc_func()
            return atoms.get_potential_energy()
        
    def make_images_and_canc_dH(self,ini_idx,fin_idx,n):
        """imageを作成する
        n: 指定方法は2種類ある
        int: n個のイメージで作成する
        Tuple(dist,max_n,min_n)
        """
        ini = self.eqs[ini_idx]
        fin = self.eqs[fin_idx]
        if type(n)!=int:
            n = get_appropriate_nimages(ini,fin,n[0],n[1],n[2],self.mic)
        images = [ini.copy() for _ in range(n+1)]+[fin.copy()]
        dH = abs(self.get_energy(ini)-self.get_energy(fin))
        interpolate(images,mic=self.mic)
        return images, dH # この地点ではCalculatorもConstraintsも設定していない
            
    def _first_push(self,i,ini_idx,fin_idx):
        ini_atoms = self.eqs[ini_idx]
        fin_atoms = self.eqs[fin_idx]
        if self.isomorphism_func(ini_atoms,fin_atoms):
            return [] # ini,finが同じ構造だった時
        self.mkdir(f"{self.logfolder}/{i}")
        images,dH = self.make_images_and_canc_dH(ini_idx,fin_idx,self.nimages)
        self.write_path(i)
        return [[-dH,[ini_idx,fin_idx,0],images,i]]
    
    def make_first_push(self,i,ini_idx,fin_idx):
        return lambda ini_idx=ini_idx,fin_idx=fin_idx: self._first_push(i,ini_idx,fin_idx)
    
    def run_and_push(self,data,layer):
        _,[ini_idx,fin_idx,width],images,i = data
        folder = f"{self.logfolder}/{i}"
        preneb_name = f"pre_EQ{ini_idx}-EQ{fin_idx}({layer}-{width})"
        neb_name = f"EQ{ini_idx}-EQ{fin_idx}({layer}-{width})"
        _, images,minimum_point = self.run_neb(images,f"{folder}/{preneb_name}",preneb=True) # PreNEB
        if len(minimum_point) > 0:
            # 極小値がある場合
            eq_num_list, atoms_list = self.run_opts(images,minimum_point,f"{folder}/opt_{preneb_name}",nabname=f"{folder}/{preneb_name}") # 構造最適化とEQの保存
            eq_num_list = [ini_idx] + eq_num_list + [fin_idx]
            atoms_list = [images[0]] + atoms_list + [images[-1]]
            push_data = self.make_push_data(eq_num_list,atoms_list,i) # 新たなImagesを作成
            return push_data
        else:
            if self.rename_preneb_file: # pre_NEBのrename
                self.rename(f"{folder}/{preneb_name}",f"{folder}/{neb_name}")
            converged,images,minimum_point = self.run_neb(images,f"{folder}/{neb_name}") # PreNEB
            if len(minimum_point) > 0:    
                if not converged or self.separate_completely:
                    eq_num_list, atoms_list = self.run_opts(images,minimum_point,f"{folder}/opt_{neb_name}",nabname=f"{folder}/{neb_name}") # 構造最適化とEQの保存
                    eq_num_list = [ini_idx] + eq_num_list + [fin_idx]
                    atoms_list = [images[0]] + atoms_list + [images[-1]]
                    push_data = self.make_push_data(eq_num_list,atoms_list,i) # 新たなImagesを作成
                    return push_data
                else:
                    self.register_tss(images,f"{i}_EQ{ini_idx}-EQ{fin_idx}({layer}-{width})",i,ini_idx,fin_idx)
                    return []
            else:
                if converged:
                    self.register_tss(images,f"{i}_EQ{ini_idx}-EQ{fin_idx}({layer}-{width})",i,ini_idx,fin_idx)
                    return []
                else:
                    write(f"{self.unconv_folder}/({i}){neb_name}.traj",images) # 未収束の保存
                    return []
    
    def write_path(self,i):
        path_path = f"{self.logfolder}/{i}/Path.dat"
        with open(path_path,"w") as f:
            f.write("[Dir]\n")
            f.write(f"EQ: ./{self.eq_folder}/\n")
            f.write(f"TS: ./{self.ts_folder}/\n")
            f.write(f"[{self.path_title[i]}]\n")
            f.write(f"path: {self.path_txt[i]}\n")
        path_obj = create_reactpath(path_path)
        path_obj.write_html(f"{self.logfolder}/{i}/Path.html")
        
    def write_all_path(self):
        p = Path(self.path_file)
        with open(p,"w") as f:
            f.write("[Dir]\n")
            f.write(f"EQ: ./{self.eq_folder}/\n")
            f.write(f"TS: ./{self.ts_folder}/\n")
            for title,path in zip(self.path_title,self.path_txt):
                f.write(f"[{title}]\n")
                f.write(f"path: {path}\n")
        paths_obj = create_reactpath(self.path_file)
        paths_obj.write_html(outfile=p.with_suffix(".html"))
        
    def make_path_txt(self,path_txt,ini_name,fin_name,add_data):
        """
        path_txt: 現在のpath_txt
        add_data: str
            EQを追加する
        add_data: int:
            TSを追加する
        add_data: None
            バリアレスを追加
        """
        if add_data is None:
            new_path_txt = re.sub(f"(EQ{ini_name}).(EQ{fin_name})",f"\\1~\\2",path_txt)
        elif type(add_data) == int:
            new_path_txt = re.sub(f"(EQ{ini_name}).(EQ{fin_name})",f"\\1;TS{add_data};\\2",path_txt)
        elif type(add_data)==str:
            new_path_txt = re.sub(f"EQ{ini_name}.EQ{fin_name}",f"{add_data}",path_txt)
        return new_path_txt
                
    def before_run(self,restart:int):
        if restart != 0:
            config  = ConfigParser(interpolation= ExtendedInterpolation())
            config.optionxform = str # キーが小文字になることを防ぐ
            config.read(self.path_file, encoding='utf-8')
            title_list = config.sections()
            if "Dir" in title_list:
                title_list.remove("Dir")
            self.path_title = title_list
            self.path_txt = [config[title]["path"] for title in title_list]
            self.eqs,_ = read_traj(self.eq_folder)
            self.tss,_ = read_traj(self.ts_folder)
        else:
            for i in range(len(self.eq_list)):
                write(f"{self.eq_folder}/{i}.traj",self.eq_list[i])
                self.eqs.append(self.eq_list[i])   
            self.path_title = [f"{ts_idx}(EQ{ini_idx}-EQ{fin_idx})" for ts_idx,(ini_idx,fin_idx) 
                        in enumerate(self.connections[restart:],start=restart)]
            self.path_txt = [f"EQ{ini_idx};EQ{fin_idx}" 
                        if self.isomorphism_func(self.eqs[ini_idx],self.eqs[fin_idx])
                        else f"EQ{ini_idx}|EQ{fin_idx}" 
                        for ini_idx,fin_idx in self.connections[restart:]]
        self.write_all_path()
                    
    def run_neb(self,images,name,preneb=False,):
        self.neb = RepetitiveNEB(images,self.constraints,self.parallel,self.mic,self.calc_func,name)
        self.neb.attach(lambda:write_html(f"{name}.html",self.neb.images))
        self.neb.attach(lambda:write(f"{name}.traj",self.neb.images))
        param = self.preneb_param if preneb else self.neb_param
        converged = self.neb.run(**param)
        minimum_point,_ = find_local_minimum_points_and_ea(self.neb.images,minenergy=self.minenergy)
        write_html(f"{name}.html",self.neb.images,highlight=minimum_point)
        return converged,self.neb.images,minimum_point
    
    def run_opt(self,atoms,name):
        self.opt = RepetitiveOpt(atoms,self.constraints,name,self.calc_func)
        param = self.opt_param
        converged = self.opt.run(**param)
        eq_num = self.register_eqs(self.opt.atoms) if converged else None
        return [eq_num,self.opt.atoms] # 収束した場合eq_num=EQ番号,未収束の場合None
    
    def run_opts(self,images,minimum_point,name,nabname):
        # 最適化を行なう. [EQ番号] 未収束の場合EQ番号はNoneが入る
        eq_num_list,atoms_list,minimum_point = zip(*[self.run_opt(images[idx],f"{name}-({idx}).log")+[idx] for idx in minimum_point])
        write_html(f"{nabname}.html",images,highlight=minimum_point,     
            annotation=[(i,f"Succeed-->EQ{eq_num}" if eq_num is not None else "Failed") 
                         for i,eq_num in zip(minimum_point,eq_num_list)])
        return list(eq_num_list), list(atoms_list)

    def make_push_data(self,eq_num_list, atoms_list,ts_idx): # ts_idxはGRRM上のTS番号
        push_data = []
        path_txt = f"EQ{eq_num_list[0]}"
        for i,((ini_eq_num,ini_atoms),(fin_eq_num,fin_atoms)) in enumerate(partially_overlapping_groupings(zip(eq_num_list,atoms_list),2,1)):
            if not self.isomorphism_func(ini_atoms,fin_atoms): # 違う構造の時
                images,dH = self.make_images_and_canc_dH(ini_eq_num,fin_eq_num,self.nimages)
                push_data.append([-dH,[ini_eq_num,fin_eq_num,i],images,ts_idx])
                path_txt += f"|EQ{fin_eq_num}"
            else:
                path_txt += f";EQ{fin_eq_num}"
        new_path_txt = self.make_path_txt(self.path_txt[ts_idx],eq_num_list[0],eq_num_list[-1],path_txt)
        self.path_txt[ts_idx] = new_path_txt
        self.write_path(ts_idx)
        return push_data
    
    def run(self,restart:int=0):
        self.before_run(restart)
        for i,(ini_idx,fin_idx) in enumerate(tqdm(self.connections[restart:]),start=restart):
            if isinstance(ini_idx,str) or isinstance(fin_idx,str):
                continue
            self.first_push = self.make_first_push(i,ini_idx,fin_idx)
            que = Queue(executor=self,priority=self.priority,maxjob=self.maxjob)
            que.run()
            self.write_all_path()
            if Path("STOP.txt").exists():
                Path("STOP.txt").unlink()
                break

class PinNEB(AutoPinNEB):
    def __init__(self, 
                 images, 
                 constraints=[], 
                 parallel=True,
                 mic=None, 
                 calc_func=pfp_calculator, 
                 eq_folder="EQ_List",
                 logfolder = "log",
                 images_folder="Images_List",
                 unconv_folder="unconv",
                 ts_folder="TS_List",
                 errorfile="ERROR",
                 path_html="Path.html",
                 path_dat = "Path.dat",
                 priority = 0,
                 maxjob = 20,):
        self.images = images
        self.constraints = constraints
        self.parallel = parallel
        if mic is None:
            mic = True if any(images[0].get_pbc()) else False
        self.mic = mic
        self.calc_func = calc_func
        self.eq_folder = eq_folder
        self.images_folder = images_folder
        self.unconv_folder = unconv_folder
        self.ts_folder = ts_folder
        self.logfolder = logfolder
        self.errorfile = errorfile
        self.path_html = path_html
        self.path_dat = path_dat
        self.priority = priority
        self.maxjob = maxjob
        self.mkdir(f"{self.logfolder}")
        self.mkdir(f"{self.ts_folder}")
        self.mkdir(f"{self.eq_folder}")
        self.mkdir(f"{self.images_folder}")
        self.mkdir(f"{self.unconv_folder}")
        
        self.path_txt = ""
        self.eqs = []
        self.tss = []
        self.param = self.default_params
        self.set_param(**self.param)
        
    @property
    def default_params(self):
        return super().default_params
    
    @property
    def preneb_param(self):
        return super().preneb_param
    
    @property
    def neb_param(self):
        return super().neb_param
    
    @property
    def opt_param(self):
        return super().opt_param
    
    def set_param(self,**param):
        super().set_param(**param)
            
    def get_param(self):
        return super().get_param()
        
    def register_eqs(self,atoms):
        self.eqs.append(atoms)
        eq_num = len(self.eqs)-1
        write(f"{self.eq_folder}/{eq_num}.traj",atoms)
        return eq_num
        
    def register_tss(self,images,name,ini_idx,fin_idx): # ts_idx GRRMのTS番号
        imax = get_imax(images,self.barrier_less,self.calc_func)
        if imax is not None: # バリアレスでない場合
            ts = images[imax]
            self.tss.append(ts)
            write(f"{self.ts_folder}/{len(self.tss)-1}.traj",ts)
            new_path_txt = self.make_path_txt(self.path_txt,ini_idx,fin_idx,len(self.tss)-1)
        else: # バリアレスの場合
            new_path_txt = self.make_path_txt(self.path_txt,ini_idx,fin_idx,None)
        self.path_txt = new_path_txt
        self.write_path()
        write(f"{self.images_folder}/{name}.traj",images)            

    def make_path_txt(self,path_txt,ini_name,fin_name,add_data):
        """
        path_txt: 現在のpath_txt
        add_data: str
            EQを追加する
        add_data: int:
            TSを追加する
        add_data: None
            バリアレスを追加
        """
        if add_data is None:
            new_path_txt = re.sub(f"(EQ{ini_name}).(EQ{fin_name})",f"\\1~\\2",path_txt)
        elif type(add_data) == int:
            new_path_txt = re.sub(f"(EQ{ini_name}).(EQ{fin_name})",f"\\1;TS{add_data};\\2",path_txt)
        elif type(add_data)==str:
            new_path_txt = re.sub(f"EQ{ini_name}.EQ{fin_name}",f"{add_data}",path_txt)
        return new_path_txt
    
    def write_path(self):
        with open(self.path_dat,"w") as f:
            f.write("[Dir]\n")
            f.write(f"EQ: ./{self.eq_folder}/\n")
            f.write(f"TS: ./{self.ts_folder}/\n")
            f.write(f"[reaction]\n")
            f.write(f"path: {self.path_txt}\n")
        path_obj = create_reactpath(self.path_dat)
        path_obj.write_html(self.path_html)
    
    def first_push(self):
        """必須メソッド"""
        self.path_txt = "EQ0|EQ1"
        ini = self.images[0]
        fin = self.images[-1]
        if self.isomorphism_func(ini,fin):
            return [] # ini,finが同じ構造だった時
        dH = abs(self.get_energy(ini)-self.get_energy(fin))
        self.write_path()
        return [[-dH,[0,1,0],self.images]]
    
    def make_push_data(self,eq_num_list,atoms_list):
        push_data = []
        path_txt = f"EQ{eq_num_list[0]}"
        for i,((ini_eq_num,ini_atoms),(fin_eq_num,fin_atoms)) in enumerate(partially_overlapping_groupings(zip(eq_num_list,atoms_list),2,1)):
            if not self.isomorphism_func(ini_atoms,fin_atoms): # 違う構造の時
                images,dH = self.make_images_and_canc_dH(ini_eq_num,fin_eq_num,self.nimages)
                push_data.append([-dH,[ini_eq_num,fin_eq_num,i],images])
                path_txt += f"|EQ{fin_eq_num}"
            else:
                path_txt += f";EQ{fin_eq_num}"
        new_path_txt = self.make_path_txt(self.path_txt,eq_num_list[0],eq_num_list[-1],path_txt)
        self.path_txt = new_path_txt
        self.write_path()
        return push_data 

    def run_and_push(self,data,layer):
        _,[ini_idx,fin_idx,width],images = data
        preneb_name = f"pre_EQ{ini_idx}-EQ{fin_idx}({layer}-{width})"
        neb_name = f"EQ{ini_idx}-EQ{fin_idx}({layer}-{width})"
        _, images,minimum_point = self.run_neb(images,f"{self.logfolder}/{preneb_name}",preneb=True) # PreNEB
        if len(minimum_point) > 0:
            # 極小値がある場合
            eq_num_list, atoms_list = self.run_opts(images,minimum_point,f"{self.logfolder}/opt_{preneb_name}",nabname=f"{self.logfolder}/{preneb_name}") # 構造最適化とEQの保存
            eq_num_list = [ini_idx] + eq_num_list + [fin_idx]
            atoms_list = [images[0]] + atoms_list + [images[-1]]
            push_data = self.make_push_data(eq_num_list,atoms_list) # 新たなImagesを作成
            return push_data
        else:
            if self.rename_preneb_file: # pre_NEBのrename
                self.rename(f"{self.logfolder}/{preneb_name}",f"{self.logfolder}/{neb_name}")
            converged,images,minimum_point = self.run_neb(images,f"{self.logfolder}/{neb_name}") # PreNEB
            if len(minimum_point) > 0:    
                if not converged or self.separate_completely:
                    eq_num_list, atoms_list = self.run_opts(images,minimum_point,f"{self.logfolder}/opt_{neb_name}",nabname=f"{self.logfolder}/{neb_name}") # 構造最適化とEQの保存
                    eq_num_list = [ini_idx] + eq_num_list + [fin_idx]
                    atoms_list = [images[0]] + atoms_list + [images[-1]]
                    push_data = self.make_push_data(eq_num_list,atoms_list) # 新たなImagesを作成
                    return push_data
                else:
                    self.register_tss(images,f"EQ{ini_idx}-EQ{fin_idx}({layer}-{width})",ini_idx,fin_idx)
                    return []
            else:
                if converged:
                    self.register_tss(images,f"EQ{ini_idx}-EQ{fin_idx}({layer}-{width})",ini_idx,fin_idx)
                    return []
                else:
                    write(f"{self.unconv_folder}/{neb_name}.traj",images) # 未収束の保存
                    return []

    def run(self):
        self.register_eqs(self.images[0])
        self.register_eqs(self.images[-1])
        que = Queue(executor=self,priority=self.priority,maxjob=self.maxjob)
        que.run()
        
    
