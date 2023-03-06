from ase.optimize import FIRE,LBFGS
from ase.io import write
from ase.neb import NEB,interpolate
from ase.units import kJ, mol
import numpy as np
from scipy.signal import argrelmax,argrelmin
# USER
from grrmpy.calculator import pfp_calculator
from grrmpy.neb.functions import to_html_nebgraph
from grrmpy.functions import copy_images, set_calc2images, set_const2images, get_energy_list, n_sampling
from grrmpy.optimize.attach import automate_maxstep,same_energy_and_fmax,steep_peak
from grrmpy.optimize.optimizer import OptWithCond

NEB_MAXSTEP_CLIMB_FALSE = {10:0.01,5:0.05,2:0.07,1:0.1,0.1:0.2,0.07:0.3,0.06:0.4,0:0.5,}
NEB_MAXSTEP_CLIMB_TRUE = {1:0.05,0:0.01,}
OPT_MAXSTEP = {10:0.1,5:0.2,2:0.3,0:0.35,}


class RepetitiveNEB():
    """
    | 一つのNEBイメージに対して,異なる条件で複数回NEB計算を行なうためのクラス.
    | デフォルトではclimb,maxsteps,fmax,stepsなど,NEBやOptimizerの引数を変えられる.
    
    Parameters:
    
    images: list of Atoms
        | Atomsのリスト
    constraints: ASE Constraint
        | constrainオブジェクト,複数制約をかける場合リストで与える.
    parallel: bool
        | 並列計算を行なう場合True.
    mic: bool
        | 最小イメージ規則にする場合True. デフォルトでは周期境界条件の場合自動でTrueにする.
    calc_func: functions object
        | claculatorを返す関数
    name:
        | 指定した名前でlogファイルが作成される.
            
    Properties:
    
    | RepetitiveNEBのプロパティ(の一部).
    | このクラスを継承するし, 新たなクラスを作成することでカスタマイズされたNEBを作成することができる.
    | 例えば,updateにイメージを切断,挿入する関数を与えることでSNEBは作成された.
        
    images: list of Atoms
        | 現在実行中のNEBのimages
    neb: NEB
        | 現在実行中のneb
    opt: Optimizer
        | 現在実行中のOptimizer
    converged: bool
        | 直近のNEBの収束状況
    """
    _default_param = {
            "base_optimizer":FIRE,
            "climb_list":[False,False,True],
            "maxstep_list":[0.03,0.1,None],
            "fmax_list":[0.05,0.05,0.05],
            "steps_list":[200,2000,2000],
            "neb_kwargs":{},
            "opt_kwargs":{},
            "automate_maxstep_true":NEB_MAXSTEP_CLIMB_TRUE,
            "automate_maxstep_false":NEB_MAXSTEP_CLIMB_FALSE,
            "same_energy_and_fmax_arg":{"times":10,"row":2,"force_digit":4,"energy_digit":6},
            "steep_peak_arg":{"max_ea":300, "nsteps":100},
        }
    
    def __init__(self,images, constraints=[], parallel=True, mic=None, calc_func=pfp_calculator,name=None,):
        """異なる条件で繰り返しNEB計算を行なう

        Parameters:
        
        images: list of Atoms
            Atomsのリスト
        constraints: ASE Constraint
            constrainオブジェクト,複数制約をかける場合リストで与える.
        parallel: bool
            並列計算を行なう場合True.
        mic: bool
            最小イメージ規則にする場合True. デフォルトでは周期境界条件の場合自動でTrueにする.
        calc_func: functions object
            claculatorを返す関数
        name:
            指定した名前でlogファイルが作成される.
        """
        self.images = copy_images(images)
        self.calc_func = calc_func
        self.constraints = constraints
        set_calc2images(self.images,self.calc_func)
        set_const2images(self.images, self.constraints)
        
        self.name=self.__class__.__name__ if name is None else name
        self.parallel = parallel
        self.constraints = constraints
        if mic is None:
            mic = True if any(self.images[0].get_pbc()) else False
        self.mic = mic
        self.attach_list = []
        self.update_list = []
        self.attach_stop_list = []
    
    @property
    def default_param(self):
        """run()実行時のデフォルトの引数"""
        return self._default_param
    
    def check_param(self):
        neb = ["climb_list","maxstep_list","fmax_list","steps_list"]
        neb_len = [len(getattr(self, i)) for i in neb]
        if not all([i==neb_len[0] for i in neb_len]):
            raise KeyError(f"{', '.join(neb)}は同じ要素数である必要があります")

    def set_param(self,**param):
        default_param = self.default_param
        unexpected_argument = set(param)-set(default_param)
        if len(unexpected_argument) != 0:
            raise Exception(f"予期しない引数({', '.join(unexpected_argument)})が存在します.\n"+
                            f"'self.default_param??'で引数を確認できます\n"
                           f"Parameters:\n{', '.join(self.default_param.keys())}")
        default_param.update(param)
        for key,val in default_param.items():
            setattr(self, key, val)
                
    def run_neb(self, climb_list, maxstep_list, fmax_list, steps_list):
        for climb,maxstep,fmax,steps in zip(climb_list,maxstep_list,fmax_list,steps_list):
            self.neb = NEB(self.images,climb=climb,parallel=self.parallel,**self.neb_kwargs)
            self.opt = OptWithCond(self.base_optimizer, self.neb,
                                   logfile=f"{self.name}.log",
                                   maxstep=maxstep if maxstep else 0.2)
            self.opt.attach_stop(same_energy_and_fmax(self.opt,**self.same_energy_and_fmax_arg))
            self.opt.attach_stop(lambda:steep_peak(self.opt,**self.steep_peak_arg))
            if maxstep is None:
                self.opt.attach(
                    lambda:automate_maxstep(
                        self.opt,
                        self.automate_maxstep_true if climb else self.automate_maxstep_false))
            for args,kwargs in self.attach_list:
                self.opt.attach(*args,**kwargs)
            for args,kwargs in self.attach_stop_list:
                self.opt.attach_stop(*args,**kwargs)
            self.converged = self.opt.run(fmax=fmax,steps=steps)
            self.do_update()
        return self.converged
    
    def attach(self,*args,**kwargs):
        """ASEのOptimizerと同じ"""
        self.attach_list.append((args,kwargs))
        
    def attach_stop(self,*args,**kwargs):
        """
        | True or Falseを与える関数を設定する
        | ASEのoptimizerのようにイタレーション毎に与えた関数が実行される.
        | 関数の戻り値がTrueになった場合, 計算が停止する.
        """
        self.attach_stop_list.append((args,kwargs))
    
    def update(self,function):
        """1回のNEB計算が終了する毎に実行される関数を指定できる"""
        self.update_list.append(function)
        
    def do_update(self):
        for function in self.update_list:
            function()
        
    def run(self,**params):
        """
        | 計算を実行する
        | 設定可能な引数がdefault_paramで確認できる.
        
        Parameters:
        
        base_optimizer: Optimizer
            | NEB計算に使用するOptimizer.
            | ここで指定したOptimizerを継承して作成されるOptWithCondオプティマイザーが実際には使用される.
        climb_list: list of bool
            毎回のNEB計算のclib
        maxstep_list: list of float
            | 毎回のNEB計算のmaxstep_list
            | 要素にNoneを指定した場合,forceに合わせて自動でmaxstepを変化させる.
            | どのように変化させるかはautomate_maxstep_true,automate_maxstep_falseで指定できる.
        fmax_list: list of float
            | 毎回のNEB計算のfmax
        steps_list: list of int
            | 毎回のNEB計算のsteps
        neb_kwargs: dict
            | NEBのインスタンス時の引数を指定したい場合に用いる
            | 例えばばね定数kを変えたい場合{"k"=0.2}のようにする
            | ここで指定した値は全てのNEB計算に適用される
        opt_kwargs: dict
            | Optimizerのインスタンス引数を指定したい場合に用いる.
        automate_maxstep_true: dict
            | climb=True時にautomaxstepを使用した際のmaxstep値
        automate_maxstep_false: dict
            | climb=False時にautomaxstepを使用した際のmaxstep値
        same_energy_and_fmax_arg: dict
            | 停止基準
        steep_peak_arg:
            | 停止基準
        """
        self.set_param(**params)
        self.check_param()
        return self.run_neb(self.climb_list, self.maxstep_list, self.fmax_list, self.steps_list)
    
    
class CutNEB(RepetitiveNEB):
    def __init__(self,images,constraints=[],parallel=True,calc_func=pfp_calculator,**kargs):
        from grrmpy.io.html import write_html #循環インポートになるため
        super().__init__(images,constraints,parallel,calc_func,**kargs)
        self.barrierless = False
        self.attach(lambda:write(f"{self.name}.traj",self.images))
        self.attach(lambda:write_html(f"{self.name}.html",self.images))
            
    def adjust_images(self,images,final_nimages,reverse=False):
        """
        imagesに新らたにimageを追加したり,削除したりしてfinal_nimages+1個のイメージになるように調節する
        """
        if final_nimages+1 == len(images):
            return copy_images(images)
        if final_nimages+1 < len(images):
            # イメージを削除する
            new_images = copy_images(n_sampling(images,final_nimages+1,True))
        else:
            add_n = final_nimages-len(images)+1 # 新たに追加するimages数
            box = len(images)-1 # imageとimegeの間の数(この間に新たなimageを挿入する)
            basic_n = add_n // box
            extra_n = add_n % box
            insert_list = [basic_n for _ in range(box)]
            for i in range(extra_n):
                i = i if reverse else -(i+1)
                insert_list[i] += 1
            new_images = []
            for i,n in enumerate(insert_list):
                p_images = [images[i].copy() for _ in range(n+1)]+[images[i+1].copy()]
                interpolate(p_images,self.mic)
                new_images+=p_images[:-1]
            new_images+=[images[-1]]
        return new_images

    def cut_images(self,threshold,nimages):
        energy_list = get_energy_list(self.images,self.calc_func)
        energy_list = np.array([energy*mol/kJ for energy in energy_list]) # 単位換算&np.array化
        if len(argrelmax(energy_list)[0])==0 and len(argrelmin(energy_list)[0])==0:
            self.barrierless = True
            return # バリアレスの場合
        
        # 削除するindex番号を調べる
        ts_i = energy_list[1:-1].argmax()+1 # 両端以外で最大値を探す
        ini_e = energy_list[:ts_i].min()
        ts_e = energy_list[ts_i]
        fin_e = energy_list[ts_i:].min()
        ini_threshold_val = (ts_e-ini_e)*threshold*0.01 + ini_e
        fin_threshold_val = (ts_e-fin_e)*threshold*0.01 + fin_e
        
        ini_del_idxes = [i for i,b in enumerate(energy_list<ini_threshold_val) if b and i<ts_i] # 削るindex番号
        fin_del_idxes = [i for i,b in enumerate(energy_list<fin_threshold_val) if b and i>ts_i] # 削るindex番号

        # いき過ぎたイメージの削除を阻止する
        if ts_i-self.neighborhood <= 0:
            ini_del_idxes = [] # 削除するindex番号
        elif ts_i-self.neighborhood in ini_del_idxes:
            ini_del_idxes = [i for i in range(0,ts_i-self.neighborhood)]
            
        if ts_i+self.neighborhood >= len(energy_list)-1:
            fin_del_idxes = []
        elif ts_i+self.neighborhood in fin_del_idxes:
            fin_del_idxes = [i for i in range(ts_i+self.neighborhood+1,len(energy_list))]
   
        self.new_ini_idx = ini_del_idxes[-1]+1 if ini_del_idxes != [] else 0 # write_progress_htmlに伝えるため
        self.new_fin_idx = fin_del_idxes[0]-1 if ini_del_idxes != [] else len(self.images)-1
        # 両端を削る
        ini_deleted_images = [image for i,image in enumerate(self.images[:ts_i+1]) if not i in ini_del_idxes]
        fin_deleted_images = [image for i,image in enumerate(self.images[ts_i:],start=ts_i) if not i in fin_del_idxes]
        # nimagesに合わせ新たにimagesを追加または削除
        ini_images = self.adjust_images(ini_deleted_images,int((nimages-1)/2),reverse=True)
        fin_images = self.adjust_images(fin_deleted_images,int((nimages-1)/2))
        self.images = ini_images + fin_images[1:]
        set_calc2images(self.images,self.calc_func)
        set_const2images(self.images,self.constraints)
                
    @property
    def default_param(self):
        param = super().default_param
        this_param = {
            "nimages_list":[13 for _ in range(len(param["climb_list"])-1)], #親クラスのデフォルト値が変更されてもいいように
            "threshold_list":[10 for _ in range(len(param["climb_list"])-1)],
            "neighborhood":1,
            }
        param.update(this_param)
        return param
    
    def write_progress_html(self,text=True):
        with open(f"{self.name}_progress.html","a") as f:
            html_text = to_html_nebgraph(self.neb,self.calc_func,False,include_plotlyjs="cdn")
            if text:
                f.write(f"<p>{self.new_ini_idx},{self.new_fin_idx}</p>")
            f.write(html_text)
    
    def check_param(self):
        super().check_param()
        if len(self.climb_list) != len(self.nimages_list)+1 or len(self.climb_list) != len(self.threshold_list)+1:
            raise Exception(f"len(climb_list)==len(nimages_list)+1==len(threshold_list)+1を満たす必要があります")
        if any([i%2==0 for i in self.nimages_list]):
            raise Exception(f"nimages_listは奇数のリストです")
        
    def set_param(self,**params):
        super().set_param(**params)
        image_change_iter = iter(
            [lambda:self.cut_images(threshold,nimages)
             for threshold,nimages in zip(self.threshold_list,self.nimages_list)]+
            [lambda:None])
        self.update(lambda:next(image_change_iter)())
        write_html_iter = iter(
            [lambda:self.write_progress_html(True) for _ in self.threshold_list]
            + [lambda:self.write_progress_html(False)]
            )
        self.update(lambda:next(write_html_iter)())

class SNEB(RepetitiveNEB):
    def __init__(self,**kargs):
        super().__init__(**kargs)
        
    def updata_images(self,images):
        return images

    def check_param(self):
        super().check_param()
        if len(self.climb_list) != len(self.nimages_list)+1:
            raise KeyError(f"len(climb_list)==len(nimages_list)+1を満たす必要があります")
        
    def set_param(self,**params):
        super().set_param(**params)
        attach_list = [{"function":lambda:write(f"{self.name}.traj",self.images)},
                       {"function":lambda:write_html(f"{self.name}.html",self.images)}]
        self.attach_list += attach_list