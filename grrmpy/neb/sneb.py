from ase.optimize import FIRE
from ase.neb import NEB,interpolate
# USER
from grrmpy.neb.functions import is_barrier_less,_separate_images
from grrmpy.optimize.optimizer import OptWithCond
from grrmpy.optimize.attach import same_energy_and_fmax,steep_peak
from grrmpy.calculator import pfp_calculator
from grrmpy.optimize.attach import automate_maxstep
from grrmpy.neb.functions import to_html_nebgraph


class SNEB():
    def __init__(
        self, 
        *data,
        constraints=[],
        html = "progress.html",
        logfile = "-",
        trajectory = None,
        mic = None,
        parallel = True,
        calc_func = pfp_calculator,
        ):
        """

        Parameters:
        
        data:
            3つの方法で指定できる
            - 方法1(2つ引数を与える)
                ini,finのイメージを指定(イメージ数は13個になる)
            - 方法2(3つ引数を与える)
                ini,fin,イメージ数  
            - 方法3(1つの引数を与える)
                images
        constraints: list of constraint
            constraintsを設定する場合constraintのオブジェクトを指定する.
            複数指定する場合はリストで与える
        html: str
            htmlファイル名
        logfile:
            logファイルの出力ファイル名
        trajectory:
            trajファイルの出力ファイル名
        mic: bool
            Noneの場合,dataが周期境界の場合自動でTrueにする.
        parallel: bool
            NEBの際に並列計算をする場合True.
        calc_func: functions object
            calculatorを返す関数
        """
        # 引数の整理
        self.constraints = constraints
        self.calc_func = calc_func
        self.html = html
        self.logfile = logfile
        self.trajectory = trajectory
        self.barrier_less = None
        self.parallel = parallel
        
        # imagesの作成
        if len(data) == 1:
            self.images = data[0]
            if mic is None:
                self.mic = True if all(self.images[0].get_pbc()) else False
            else:
                self.mic = mic
        else:
            ini = data[0]
            fin = data[1]
            n_images = 13 if len(data) == 2 else data[2] # デフォルトのイメージ数は13
            self.images = [ini.copy() for _ in range(n_images+1)]+[fin.copy()]  
            if mic is None:
                self.mic = True if all(ini.get_pbc()) else False
            else:
                self.mic = mic
            interpolate(self.images,mic=self.mic)
        for image in self.images:
            del image.constraints
            image.set_constraint(self.constraints)
            image.calc = self.calc_func()
                
    def write_html(self,comment,graph=True,**kargs):
        if self.html:
            with open(self.html,"a") as f:
                f.write(comment)
                if graph:
                    html_text = to_html_nebgraph(
                        self.neb,
                        self.calc_func,
                        False,
                        include_plotlyjs="cdn",
                        **kargs)
                    f.write(html_text)
                        
    def update_images(self,nimages,tolerance,threshold):
        """
        | バリアレスだった場合,更新はしないでNoneを返す
        | イメージを削った場合,削った部分のリスト(新たにini,finとするidx)を返す
        """
        if is_barrier_less(self.images, self.barrier_less_energy, calc_func=self.calc_func):
            return None
        # imageを削る
        self.images, ini_idx, fin_idx = _separate_images(self.images,tolerance,threshold)
        # 削ったimagesの間にimageを追加
        
        
        return [ini_idx,fin_idx]
            
    def irun(self,fmax,steps,maxstep,climb):
        self.neb = NEB(self.images,climb=climb,parallel=self.parallel,**self.neb_kwargs)
        self.opt = OptWithCond(
            self.optimizer,
            self.neb,
            logfile=self.logfile,
            trajectory=self.trajectory,
            maxstep = maxstep if maxstep else 0.2,
            **self.opt_kwargs,
            )
        self.opt.attach_stop(steep_peak(self.opt,self.max_ea,self.nsteps))
        self.opt.attach_stop(same_energy_and_fmax(self.opt,self.times,self.row))
        if maxstep:
            ms = self.maxstep_true if climb else self.maxstep_false
            self.opt.attach(lambda:automate_maxstep(self.opt,maxstep=ms))
        for attach in self.attach_list:
            self.opt.attach(**attach)
        converged = self.opt.run(fmax=fmax,steps=steps)
        return converged
    
    def run_neb(self,climb,maxstep):
        """recalcで再度NEB計算を行なう"""
        self.neb = NEB(self.images,climb=climb,parallel=self.parallel,**self.neb_kwargs)
        self.opt = OptWithCond(
            self.optimizer,
            self.neb,
            logfile=self.logfile,
            trajectory=self.trajectory,
            maxstep = 0.2,
            **self.opt_kwargs)
        self.opt.attach_stop(steep_peak(self.opt,self.max_ea,self.nsteps))
        self.opt.attach_stop(same_energy_and_fmax(self.opt,self.times,self.row))
        self.opt.attach(lambda:automate_maxstep(self.opt,maxstep=maxstep))
        for attach in self.attach_list:
            self.opt.attach(**attach)
        converged = self.opt.run(fmax=self.fmax_list[-1],steps=self.recalc_steps)
        return converged
    
    @property
    def imax(self):
        self.neb.imax
    
    @property
    def default_param(self):
        """SEB.run()の引数
        
        基本的な引数:
        
        optimizer: obj
            NEB計算を行なう際のOptimizer
        nimages_list: list of int
            各NEB計算をいくつのイメージ数で行なうか
        maxstep_list: list of float
            | 各NEB計算のmaxstep. Noneを要素とした場合,forceに合わせて自動で変更する.
            | どのように自動化するのはauto_maxstep_climb_false,auto_maxstep_climb_tureで指定できる
        climb: bool
            | Falseにした場合,最後のNEBのみclimb=Trueとする.
            | Trueにした場合,各NEBでFalseもTrueも計算する
        
        収束基準:
        
        fmax_list: list of float
            各NEB計算のfmax値
        steps_list: list of int
            各NEB計算のstep数
        max_ea: float
            max_ea(kJ/mol)以上のTSになった瞬間,計算を停止する.
        nsteps: int
            | max_ea(kJ/mol)以上のTSになった瞬間,計算を停止するが,
            | 始めのnsteps回のイタレーションの間はこの限りではない.
        times: int
            times回連続で同じfmax,energy繰り返す事態がrow回あれば停止する.
        row: int
            times回連続で同じfmax,energy繰り返す事態がrow回あれば停止する.
            
        NEBイメージの切断:
        
        tolerance : float
            tolerance(kJ/mol)以下の谷間は無視してTSの山を分離する
        threshold: list of float
            TSの山を分離した後, さらに両端のイメージをTSの高さに対してthreshold%削る
        barrier_less_energy: float
            barrier_less_energy(kJ/mol)以下のエネルギー変化はバリアレスとする
        min_images: int
            TSを削る際に最低でも残すイメージ数
            
        再計算:
            
        recalc: boolean
            | NEBイメージを分離したことで,それまでTSが存在していたのにバリアレスに判定になってしまった場合に    
            | 初期イメージをSNEBでなくただのNEB計算を行なう場合True.
            | Falseにした場合はそこで計算を中止するため,計算は失敗した扱いとなる
        recalc_steps: int
            | recalを行なう場合のsteps数.
        recalc_fmax: float
            recalcを行なう場合のfmax
            
        カスタマイズ:
        
        attach_list: list of dictionary
            | optにattachを取り付ける場合に用いる  
            | [{function:func_obj,interval:int},]の形で指定する
        neb_kwargs: dict
            | NEBのインスタンス引数を指定できる
        opt_kwargs: dict
            | OPTのインスタンス引数を指定できる
        auto_maxstep_climb_false: dictionary
            | maxstep_listでNoneを指定した際にこの値に従ってmaxstepが決定される(climb=False条件で)
        auto_maxstep_climb_true: dictionary
            | maxstep_listでNoneを指定した際にこの値に従ってmaxstepが決定される(climb=True条件で)   
    
        """
        return {
            "optimizer": FIRE,
            "nimages_list":[13, 14, 15],
            "maxstep_list": [0.02, 0.05, None],
            "climb": False,
            "fmax_list" : [0.1, 0.05, 0.03],
            "steps_list" : [1000, 1000, 2000],
            "max_ea":300.0,
            "nsteps":30,
            "times":12,
            "row":2,
            "tolerance":5.0,
            "threshold":[10.0, 10.0],
            "barrier_less_energy":5.0,
            "min_images":3,
            "recalc":True,
            "recalc_steps":3000,
            "recalc_fmax":0.05,
            "attach_list": [],
            "neb_kwargs":{},
            "opt_kwargs":{},
            "auto_maxstep_climb_false":{10.0 : 0.01,
                                        5.0 : 0.05,
                                        2.0 : 0.07, 
                                        1.0 : 0.1,
                                        0.1 : 0.2,
                                        0.07 : 0.3,
                                        0.06 : 0.4,
                                        0.0 : 0.5,},
            "auto_maxstep_climb_true":{1.0 : 0.05,
                                       0.0 : 0.01,},
        }
                
    def set_param(self,**kwargs):
        param = self.default_param
        param.update(kwargs)
        self.optimizer = param["optimizer"]
        self.nimages_list = param["nimages_list"]
        self.maxstep_list = param["maxstep_list"]
        self.climb = param["climb"]
        self.fmax_list = param["fmax_list"]
        self.steps_list = param["steps_list"]
        self.max_ea = param["max_ea"]
        self.nsteps = param["nsteps"]
        self.times = param["times"]
        self.row = param["row"]
        self.tolerance = param["tolerance"]
        self.threshold = param["threshold"]
        self.barrier_less_energy = param["barrier_less_energy"]
        self.min_images = param["min_images"]
        self.recalc = param["recalc"]
        self.recalc_steps = param["recalc_steps"]
        self.recalc_fmax = param["recalc_fmax"]
        self.attach_list = param["attach_list"]
        self.neb_kwargs = param["neb_kwargs"]
        self.opt_kwargs = param["opt_kwargs"]
        self.maxstep_false = param["auto_maxstep_climb_false"]
        self.maxstep_true = param["auto_maxstep_climb_true"]
        self.check_param()
        self.run_param = self.gen_run_param()  
        
    def check_param(self):
        param_len = [len(self.fmax_list),
                     len(self.steps_list),
                     len(self.maxstep_list),
                     len(self.nimages_list),
                     len(self.threshold)+1,]
        if all(val == param_len[0] for val in param_len):
            return
        else:
            raise ValueError(f"len(fmax_list), len(steps_list),"+
                             " len(maxstep_list), len(nimages_list)"+
                             "len(threshold)+1 を満たす必要があります")
    
    def gen_run_param(self):
        """[[fmax,steps,maxstep,climb,updata]]に変換する"""
        run_param = []
        for fmax,steps,maxstep in zip(self.fmax_list,self.steps_list,self.maxstep_list):
            if self.climb:
                run_param.append([fmax,steps,maxstep,False,False])
                run_param.append([fmax,steps,maxstep,True,True])
            else:
                run_param.append([fmax,steps,maxstep,False,True])
        if self.climb:
            run_param[-1][-1] = False # 最後だけはimagesのupdateをしない
        else:
            # 最後はclimb=Trueで計算する
            run_param.append([self.fmax_list[-1],self.steps_list[-1],self.maxstep_list[-1],True,False])
        return run_param
    
    def run(self,**kwargs):
        """SNEB計算を実行する
        
        kwargsはSNEB.default_paramを参照する
        """
        self.set_param(**kwargs)
        # 計算
        check_has_ts = [] # バリアレスならFalseが入る
        for fmax,steps,maxstep,climb,update in self.run_param:
            converged = self.irun(fmax,steps,maxstep,climb)
            if update:
                if is_barrier_less(self.images, self.barrier_less_energy, calc_func=self.calc_func):
                    check_has_ts.append(False)
                    if any(check_has_ts):
                        # 今までTSがあったのに削ったことでバリアレスになった場合
                        self.write_html(f"{climb=},{converged=}\n{comment}")
                        break
                else:
                    ini_idx,fin_idx = self.update_images()
                    comment = f"{ini_idx}番目,{fin_idx}番目を新たなini,finとした."
                    check_has_ts.append(True)
            else:
                comment = ""
            self.write_html(f"{climb=},{converged=}\n{comment}")
        
        # 戻り値の出力and再計算(うまくいかなかった場合)
        if all(check_has_ts):
            # 全ての過程でTSが存在した場合
            self.barrier_less = False      
        elif not any(check_has_ts):
            # ずっとバリアレスだった場合
            self.barrier_less = True
        else:
            # 上のfor文でbreakしたものがここに入る
            # TSがあったのに途中からバリアレスになってしまったもの
            if self.recalc:
                converged = self.run_neb(climb=False,maxstep=self.maxstep_false)
                self.write_html(f"再計算:climb=False,{converged=}\n")
                if converged:
                    converged = self.run_neb(climb=True,maxstep=self.maxstep_true)
                    self.write_html(f"再計算:climb=True,{converged=}\n")
                if is_barrier_less(self.images, self.barrier_less_energy, calc_func=self.calc_func):
                    self.barrier_less = True
                else:
                    self.barrier_less = False
            else:
                self.barrier_less = False
                converged = False
        self.write_html(f"{self.barrier_less=}",graph=False)
        return converged