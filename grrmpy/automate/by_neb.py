from ase.neb import interpolate
from ase.optimize import FIRE, LBFGS
from ase.io import write,iread, Trajectory
from ase.vibrations import Vibrations
from ase.units import kJ,mol
from ase.geometry import find_mic
from pathlib import Path
import csv
from pprint import pprint
import numpy as np
from ase.build.rotate import minimize_rotation_and_translation

#USER
from grrmpy.calculator import pfp_calculator
from grrmpy.neb.auto_neb import SNEB
from grrmpy.io.html import write_html
from grrmpy.vibrations.functions import to_html_table_and_imode,find_ts_idx,to_html_graph
from grrmpy.functions import (to_html_energy_diagram,
                              minimize_rotation_and_translation_for_specified_indices_only,
                              connected_components)
from grrmpy.path import ReactPath
try:
    from grrmpy.optimize import FIRELBFGS
    defaultoptimizer = FIRELBFGS
except:
    defaultoptimizer = FIRE


class SinglePath():
    """
        Parameters:
        
        \*data: 
            | 3種類の方法でinput構造を与えることができる.

            ini,fin構造と初めの中間イメージの数を与える(3つの引数を与える)
                >>> sp = SinglePath(ini,fin,21,indices=[1,2,3,4]) # 21個内部にイメージを作成
            ini,fin構造を与える(2つの引数を与える.初期イメージの数は13になる)
                >>> sp = SinglePath(ini,fin,indices=[1,2,3,4])
            images構造を与える(1つの引数を与える.calculatorを付ける必要はない)
                >>> sp = SinglePath(images,indices=[1,2,3,4])
                
            Constraintsがある場合
                
                >>> f = FixAtoms(indices=[0,1,2,3])
                >>> sp = SinglePath(ini,fin,indices=[1,2,3,4],constraints=f)
                
                >>> f = FixAtoms(indices=[0,1,2,3])
                >>> c = FixBondLength(0, 1)
                >>> sp = SinglePath(ini,fin,indices=[1,2,3,4],constraints=[f,c])
                
            imageの数は奇数個を推奨する.
        indices: list
            振動数計算を行なう時に動かす原子のindex番号のリスト
        parallel: bool
            並列計算を行なう場合True.
        mic: bool
            | 最小イメージ規則を適用する場合True.
            | Noneの場合,周期境界条件で計算する場合自動でTrueにする.
            | デフォルトはTrue.
        calc_func: function object
            calculatorを返す関数. デフォルトは ``pfp_calculator``
        debug: bool
            debug.logを出力する(デバック用)
            
        Note:
            | EQ_list.traj, TS_list.trajが既にディレクトリ中にある場合,上書きされてしまうので
            | 1つフォルダ内で複数の計算を行なわないようにする!!!
    """
    def __init__(self,
                 *data,
                 indices = None,
                 parallel:bool=True,
                 mic:bool=None,
                 calc_func=pfp_calculator,
                 debug=False,
                 constraints=[]):
        """
        
        EQ構造,TS構造,PT構造はEQ_list.traj,TS_list.traj,PT_list.trajに保存される.
        CONNECTIONSの情報はTS_CONNECTIONS.csv,PT_CONNECTIONS.csvに保存される

        Parameters:
        
        *data: 
            3種類の方法でinput構造を与えることができる.
     
            ・ini,fin構造と初めの中間イメージの数を与える(3つの引数を与える)
            ・ini,fin構造を与える(2つの引数を与える.初期イメージの数は13になる)
            ・images構造を与える(1つの引数を与える.calculatorを付ける必要はない)
                   imageの数は奇数個を推奨する.
        indices: list
            振動数計算を行なう時に動かす原子のindex番号のリスト
        parallel: bool
            並列計算を行なう場合True.
        mic: bool
            | 最小イメージ規則を適用する場合True.
            | Noneの場合,周期境界条件で計算する場合自動でTrueにする.
            | デフォルトはTrue.
        calc_func: function object
            calculatorを返す関数. デフォルトは ``pfp_calculator``
        constraints:constraint object or list
            | constraintsまたはconstraintsのリスト
        debug: bool
            debug.logを出力する(デバック用)
            
        Note:
            EQ_list.traj, TS_list.trajが既にディレクトリ中にある場合,上書きされてしまうので
            1つフォルダ内で複数の計算を行なわないようにする!!!
        """
        self.parallel = parallel
        self.calc_func = calc_func
        self.indices = indices # vibrations用
        self.debug = debug
        self.constraints = constraints
        ### 初期イメージの作成 ###
        if len(data)==1:
            data_type="images"
            self.images = data[0]
        else:
            data_type="ini_fin"
            ini = data[0]
            fin = data[1]
            nimages = 13 # Deffaultのイメージ数
            if len(data) == 3:
                nimages = data[2]
            self.images = [ini.copy()]+[ini.copy() for _ in range(nimages)]+[fin.copy()]
        
        if mic is None:
            self.mic = True if any(self.images[0].get_pbc()) else False
        else:
            self.mic = mic
        
        if data_type == "ini_fin":
            interpolate(self.images,mic=self.mic,apply_constraint=False)
        for image in self.images:
            self.set_calculator(image)
            image.set_constraint(self.constraints)
        
            
    def write_connections(self,file,data_list,mode="a"):
        with open(file,mode) as f:
            writer = csv.writer(f)
            writer.writerow(data_list)
        
    def write_eq(self,atoms):
        self.eq_traj.write(atoms)
        self.eq_count += 1
        
    def write_ts(self,atoms):
        self.ts_traj.write(atoms)
        self.ts_count += 1 
         
    def write_pt(self,atoms):
        self.pt_traj.write(atoms)
        self.pt_count += 1
    
    def debug_log(self,text):
        if self.debug:
            with open("debug.log","a") as f:
                f.write(str(text)+"\n") 
    
    @property
    def default_param(self):
        """run()のデフォルトのパラメータを返す.
        
        General:
        
        stopping_criterion: int
            停止基準. SNEBの計算をstopping_criterion回行なうと停止する.
        struct_check_threshold: list of float
            | 2つの構造が同一構造であるとみなす基準.
            | a,b,c,d,e,fの6つの要素を指定する.基準1,基準2のどちらかを満たしている時,同一構造とする.
            | 基準1:
            | 電子エネルギーの差, 原子間距離のRMS誤差,原子間距離の最大距離がそれぞれ
            | a kJ/mol, b%, cÅ 以下の時, 同一構造であるとみなす
            | 基準2: 
            | エネルギーの差がd kJ/mol以下であり,indicesで指定した原子を重ね合わせた時,
            | (最小二乗法で重ねる)原子間距離のRMS誤差,原子間距離の最大距離がそれぞれ,
            | e%, fÅ 以下の時, 同一構造であるとみなす. 
        minimize_rotation_and_translation: boolean
            | 基準1で2つの構造の座標が近くなるように再配置して解析する場合Ture.
            | (最小二乗値が小さくなるよう再配置を行なう)

        SNEB:
        
        optimizer: Optimizer class
            NEB計算に使用するOptimizerの種類
        nimages: list of int
            各NEB計算の中間イメージの数
        maxstep: list of float
            | 各NEB計算でのmaxstep値のリスト
            | Noneの場合,forceに合わせて自動的にmaxstepsを変化させる.
        fmax: list of float
            各NEB計算でのfmax値のリスト
        steps: list of int
            各NEB計算でのsteps値のリスト
        tolerance: float
            NEBバンド中のtolerance(kJ/mol)以下の谷間はEQとは考えないで計算する.
        threshold: list of float
            活性化エネルギーに対してthreshold%以下の点は削り新たなNEBイメージを作成して計算する
        dist: float
            | 新たなNEBイメージを作成する際,TSの両隣に,TSからdistÅ離れた位置にイメージを作成する.
            | TS付近にイメージを密集させるために行なう.
        min_nimages: int
            | 新たなイメージを作成するにあたって,それまでのNEBを削ることになるが,
            | 少なくともmin_nimages個は確保するようにする.
        with_stop: bool
            | optimizerにstop機能を付ける場合True.
            | max_n回連続でEnergyとForceが同じ状況がtimes回あった時計算を終了する.
            | grrmpy.optimize.functions.add_stop_to()を参照
        max_n: int:
            with_stop=Trueの場合に指定.
            with_stop(add_stop_to)の引数. デフォルト12
        times: int
            with_stop=Trueの場合に指定.
            with_stop(add_stop_to)の引数. デフォルト2
        climb_steps: integer
            | climb=Falseで収束しなかった場合,その後のTrueでのNEB計算のStep数をclimb_steps回にする.
            | Noneの場合,climb=Falseの時と同じstepsで計算する.
            | 一番最後のNEB計算ではclimb_stepsに関係なくclimb=Falseの時と同じstepsで計算する.
        climb: bool
            | Falseの場合,climbはFalse,False,False,False...Trueで行なう.
            | Trueの場合はclimbはFalse,True,Fase,True,False....Trueで行なう
            | Falseにした場合,climb_stepsを設定する意味はなくなる.
            
        IRC:
        
        optimizer1: Optimizer calss
            振動数計算後のエネルギーダイアグラム解析の結果,TSがはっきり分かる場合に使用するOptimizer
        optimizer2: Optimizer calss
            振動数計算後のエネルギーダイアグラム解析の結果,TSがはっきり見えない場合に使用するOptimizer
        maxsteps:float or None
            | maxsteps値. maxstepsの値を大きく設定したとしても,始めの200回はmaxstep=0.03で計算する.
            | Noneを指定した場合,forceに合わせて自動的にmaxstepを調整する.
        dif: float
            | 振動数計算後のエネルギーダイアグラム解析の結果,
            | TS(極大値)と谷間(極小値)の差が,dif(kJ/mol)以下の時,optimizer2を使って計算する.
        calc_notop_ts: boolean
            | 振動数計算後のエネルギーダイアグラム解析の結果,TSがない(極大値がない)時,
            | IRC計算を行なわないならFalse.
        """
        return {
            "General":{
                "stopping_criterion":20,
                "struct_check_threshold":[10.0, 0.25, 2.5, 10.0, 0.25, 2.5],
                "minimize_rotation_and_translation" : True,
            },
            "SNEB":{
                "optimizer":FIRE,
                "nimages":[17, 11, 13],
                "maxstep":[0.2, 0.3, None],
                "fmax":[0.09, 0.06, 0.03],
                "steps":[500, 1000, 3000],
                "tolerance":5.0,
                "threshold":[0, 0],
                "dist":0.1,
                "min_nimages":3,
                "with_stop":True,
                "max_n":12,
                "times":4,
                "climb_steps":15,
                "climb":False,
                },
            "IRC":{
                "optimizer1":defaultoptimizer,
                "optimizer2":defaultoptimizer,
                "maxstep":None,
                "fmax":0.001,
                "steps":40000,
                "dif":2.0,
                "calc_notop_ts":True,
            },
            }
        
    def set_calculator(self,atoms):
        if not atoms.get_calculator():
            atoms.calc = self.calc_func()
            atoms.set_constraint(self.constraints) # 多分なくても良い
        
    def run_sneb(self,*images):
        self.sneb = SNEB(*images,
                            logfile=f"SNEB{self.iter_count}.log",
                            html=f"SNEB{self.iter_count}_progress.html",
                            mic=self.mic,
                            parallel=self.parallel,
                            calc_func=self.calc_func,
                            optimizer=self.neb_optimizer,
                            with_stop = self.with_stop,
                            max_n = self.max_n,
                            times = self.times,
                            constraints=self.constraints)
        self.sneb.attach(lambda:write_html(f"SNEB{self.iter_count}.html", self.sneb.images))
        self.sneb.attach(lambda:write(f"SNEB{self.iter_count}.traj", self.sneb.images))
        sneb_converged = self.sneb.run(
            nimages=self.nimages,
            maxstep=self.neb_maxstep,
            fmax=self.neb_fmax,
            steps=self.neb_steps,
            tolerance=self.tolerance,
            threshold=self.threshold,
            dist=self.neb_dist,
            min_nimages=self.min_nimages,
            climb_steps = self.climb_steps,
            climb = self.climb
        )     
        return sneb_converged
    
    def run_vib(self):
        self.vib = Vibrations(self.ts,self.indices,name=f"vib{self.iter_count}")
        self.vib.run()
        tb_text,imode = to_html_table_and_imode(self.vib,full_html=False,include_plotlyjs="cdn")
        if type(imode) == int:
            self.vib.write_mode(n=imode)
            self.vimages = [i for i in iread(f"vib{self.iter_count}.{imode}.traj")]
            fig_text = to_html_energy_diagram(
                self.vimages,
                calc_func=self.calc_func,
                full_html=False,
                unit="kJ/mol",
                title="Vibration Energy Diagram",
                xaxis_title="",
                yaxis_title=None,
                include_plotlyjs="cdn")
            success = True
        else:
            success = False
            fig_text = ""
        with open(f"vib{self.iter_count}.html","w") as f:
            f.write(tb_text+fig_text)
        return success
    
    def run_irc(self,ts_idx:int,r_use_newton:bool, f_use_newton:bool):
        self.ini = self.vimages[ts_idx-1].copy()
        self.ini.calc = self.calc_func()
        self.ini.set_constraints(self.constraints)
        self.fin = self.vimages[ts_idx+1].copy()
        self.fin.calc = self.calc_func()
        self.fin.set_constraints(self.constraints)
        r_converged = self.irun_irc(self.ini, r_use_newton, "reverse")
        f_converged = self.irun_irc(self.fin, f_use_newton, "forward")
        return [r_converged, f_converged]
        
    def irun_irc(self,atoms,use_newton,name):
        optimizer = self.optimizer2 if use_newton else self.optimizer1
        def get_args(maxsteps):
            try:
                if optimizer==FIRELBFGS:
                    arg = {"switch":0.04,"maxstep_fire":maxsteps,"maxstep_lbfgs":maxsteps}
                else:
                    arg = {"maxstep":maxsteps}
            except:
                arg = {"maxstep":maxsteps}
            return arg
        ### 始めにmaxstep=0.03で200回計算しておく
        self.irc_opt = optimizer(atoms,
                                 logfile = f"IRC{self.iter_count}_{name}.log",
                                 **get_args(0.03))
        self.irc_opt.attach(lambda:write(f"IRC{self.iter_count}_{name}.traj", atoms))
        self.irc_opt.run(fmax=self.irc_fmax, steps=200)
        ### 計算
        self.irc_opt = optimizer(atoms,
                                 logfile = f"IRC{self.iter_count}_{name}.log",
                                 **get_args(self.irc_maxstep))
        self.irc_opt.attach(lambda:write(f"IRC{self.iter_count}_{name}.traj", atoms))
        converged = self.irc_opt.run(fmax=self.irc_fmax, steps=self.irc_steps)
        return converged
    
    def _get_dif_energy(self,atoms1, atoms2):
        """エネルギー差の判定"""
        dif_e = abs(atoms1.get_potential_energy()-atoms2.get_potential_energy()) # eV単位
        dif_e *= mol/kJ # kJ/mol単位
        return dif_e
    
    def _get_check_rmse_and_maxdist(self,atoms1,atoms2):
        """RMS誤差と最大距離の判定"""
        pos1 = atoms1.get_positions()
        pos2 = atoms2.get_positions()
        d_mat = pos2 - pos1
        if self.mic:
            d_mat = find_mic(d_mat, atoms2.get_cell(), atoms2.pbc)[0]
        d_line = np.linalg.norm(d_mat,axis=1)
        rms_error = np.sqrt(np.sum(d_line**2)/len(atoms1))
        return rms_error, np.max(d_line)
    
    def check_structures(self,atoms1,atoms2):
        """2つの構造が同一構造であるか判断する.
        
            | エネルギーの差, 原子間距離のRMS誤差,原子間距離の最大距離がそれぞれ
            | a kJ/mol, b%, cÅ 以下の時, 同一構造であるとみなす. 
            | もしくはエネルギーの差がd kJ/mol以下であり,indicesの原子を重ね合わせた時,
            | 原子間距離のRMS誤差,原子間距離の最大距離がそれぞれ,
            | e%, fÅ 以下の時, 同一構造であるとみなす. 
        """
        a,b,c,d,e,f = self.struct_check_threshold
        
        dif_energy = self._get_dif_energy(atoms1, atoms2) # エネルギー差を取得
        
        atoms1_cp = atoms1.copy()
        atoms2_cp = atoms2.copy()
        mols_idx1 = [i for i in connected_components(atoms1_cp,self.indices)] # 分子を抽出
        mols_idx2 = [i for i in connected_components(atoms2_cp,self.indices)] # 分子を抽出
        if mols_idx1 == mols_idx2: # 同じ分子で構成されていた場合
            for idxs in mols_idx1:
                indices = list(idxs)
                minimize_rotation_and_translation_for_specified_indices_only(
                    atoms1_cp,
                    atoms2_cp,
                    indices
                    )
            rmse1, max_dist1 = self._get_check_rmse_and_maxdist(atoms1_cp,atoms2_cp)
            same_geo = True
        else:
            rmse1, max_dist1 = False,False
            same_geo = False
        
        atoms1_cp = atoms1.copy()
        atoms2_cp = atoms2.copy()   
        if self.mrt:
            minimize_rotation_and_translation(atoms1_cp,atoms2_cp)  
        rmse2, max_dist2 = self._get_check_rmse_and_maxdist(atoms1_cp,atoms2_cp)
        
        # 同じ分子で構成されていた場合,rmse1を返す
        rmse = rmse1 if rmse1 else rmse2
        max_dist = max_dist1 if max_dist1 else max_dist2
        
        self.debug_log(f"{dif_energy},{rmse1},{max_dist1},{rmse2},{max_dist2}")

        if all([dif_energy<a, rmse2<b, max_dist2<c]) or all([dif_energy<d, b<rmse<e, max_dist<f]):
            return True, same_geo, rmse
        else:
            return False, same_geo, rmse
    
    def analyze_connection(self,ts_n,ini_n,fin_n,ts,ini,fin):
        """CSVファイルに書き込む情報を作成する
        
        ["TS/PT","ini","fin","ini_energy","fin_energy","forward_Ea","reverce_Ea"]

        Parameters:
        ts_n (int): TS番号
        ini_n (int): iniのEQ番号
        fin_n (int): finのEQ番号
        ts (Atoms): TSのAtoms
        ini (Atoms): iniのEQのAtoms
        fin (Atoms): finのEQのAtoms
        """
        self.set_calculator(ini)
        self.set_calculator(fin)
        self.set_calculator(ts)
        ini_e = ini.get_potential_energy()*mol/kJ
        fin_e = fin.get_potential_energy()*mol/kJ
        ts_e = ts.get_potential_energy()*mol/kJ
        forward_ea = ts_e-ini_e
        reverse_ea = ts_e-fin_e
        text = [ts_n,ini_n,fin_n,ini_e,fin_e,forward_ea,reverse_ea]
        return text
    
    def create_path(self):
        """Pathオブジェクトを作成する"""
        name = self.order
        atoms = [self.atoms_dict[i] for i in self.order]
        self.path = ReactPath({"name":name,"atoms":atoms})
        self.path.write_html("Path.html")
        self.path.topkl("Path.pickle")
        
    def run(self, param=None):
        """計算を実行する

        Parameters:
        
        param: dict
            | 種々のパラメータは辞書で与える
            | デフォルトの値はget_param()で取得できる
            
        Example:
        
        IRCのoptimizer1をBFGS変えたい時
        
        >>> sp = SinglePath(ini,fin)
        >>> param = sp.default_param
        >>> param["IRC"]["optimizer1"] = BFGS
        >>> sp.run(param)
        """
        ### 番号の初期化 ###
        self.eq_count = 0
        self.ts_count = 0
        self.pt_count = 0
        
        # connection #
        self.ts_connections = []
        self.pt_connections = []
        
        ### ファイル関係 ###
        self.eq_list_file = "EQ_list.traj"
        self.ts_list_file = "TS_list.traj"
        self.pt_list_file = "PT_list.traj"
        self.ts_info_file = "TS_CONNECTIONS.csv"
        self.pt_info_file = "PT_CONNECTIONS.csv"
        
        if Path(self.eq_list_file).exists():
            raise Exception("既にEQ_list.trajファイルがあります")
        if Path(self.ts_list_file).exists():
            raise Exception("既にTS_list.trajファイルがあります")
        
        self.eq_traj = Trajectory(self.eq_list_file,mode="a")
        self.ts_traj = Trajectory(self.ts_list_file,mode="a")
        self.pt_traj = Trajectory(self.pt_list_file,mode="a")
        self.write_eq(self.images[0])
        self.write_eq(self.images[-1])
        
        ts_coulmn = ["TS","ini","fin","ini_energy","fin_energy","forward_Ea","reverce_Ea"]
        self.write_connections(self.ts_info_file, ts_coulmn, "w")
        pt_coulmn = ["PT","ini","fin","ini_energy","fin_energy","forward_Ea","reverce_Ea"]
        self.write_connections(self.pt_info_file, pt_coulmn, "w")
        
        ###パラメータ関係
        if param is None:
            param = self.default_param
        with open("PARAM.txt","w") as f:
            pprint(param, stream=f)
        self.stopping_criterion = param["General"]["stopping_criterion"]
        self.struct_check_threshold = param["General"]["struct_check_threshold"]
        self.mrt = param["General"]["minimize_rotation_and_translation"]
        self.neb_optimizer = param["SNEB"]["optimizer"]
        self.first_nimages = param["SNEB"]["nimages"][0]
        self.nimages = param["SNEB"]["nimages"][1:]
        self.neb_maxstep = param["SNEB"]["maxstep"]
        self.neb_fmax = param["SNEB"]["fmax"]
        self.neb_steps = param["SNEB"]["steps"]
        self.tolerance = param["SNEB"]["tolerance"]
        self.threshold = param["SNEB"]["threshold"]
        self.neb_dist = param["SNEB"]["dist"]
        self.min_nimages = param["SNEB"]["min_nimages"]
        self.with_stop = param["SNEB"]["with_stop"]
        self.climb_steps = param["SNEB"]["climb_steps"]
        self.max_n = param["SNEB"]["max_n"]
        self.times = param["SNEB"]["times"]
        self.climb = param["SNEB"]["climb"]
        self.optimizer1 = param["IRC"]["optimizer1"]
        self.optimizer2 = param["IRC"]["optimizer2"]
        self.irc_maxstep = param["IRC"]["maxstep"]
        self.irc_fmax = param["IRC"]["fmax"]
        self.irc_steps = param["IRC"]["steps"]
        self.irc_dif = param["IRC"]["dif"]
        self.calc_notop_ts = param["IRC"]["calc_notop_ts"]
        #QUE
        self.que = [[0,1]] # 実際には一回目はinitで作成したself.imagesを用いる
        
        self.sneb_ini = self.images[0].copy()
        self.sneb_fin = self.images[-1].copy()
        self.order = ["EQ0","EQ1"]
        self.atoms_dict = {"EQ0":self.sneb_ini,"EQ1":self.sneb_fin}

        self.first_calculation = True
        
        self.iter_count = -1
        while len(self.que)> 0 and self.stopping_criterion >= self.iter_count+1:
            self.iter_count += 1
            self.debug_log(f"{self.order}")
            self.debug_log(f"QUE:{self.que}")
            self.create_path()
            sneb_ini_idx,sneb_fin_idx = self.que.pop(0)
            self.debug_log(f"SNEB:{sneb_ini_idx}-{sneb_fin_idx}")
            
            if self.first_calculation:
                """始めのSNEB計算(始めのみinitで作成したimagesで計算する)"""
                self.first_calculation = False
                if not self.run_sneb(self.images):
                    self.debug_log(f"SNEB収束せず終了")
                    continue
            else:
                self.sneb_ini = self.atoms_dict[f"EQ{sneb_ini_idx}"]
                self.sneb_fin = self.atoms_dict[f"EQ{sneb_fin_idx}"]       
                if not self.run_sneb(self.sneb_ini,self.sneb_fin,self.first_nimages,self.constraints):
                    self.debug_log(f"SNEB収束せず終了")
                    continue
            
            imax = self.sneb.imax
            self.ts = self.sneb.images[imax].copy()
            self.ts.calc = self.calc_func()
            
            ### VIB計算と虚振動の確認 ###
            if not self.run_vib():
                self.debug_log(f"虚振動がない")
                self.write_pt(self.ts)
                continue # 虚振動がない時はSellaを行なうように今後したい

            ### VIB エネルギーダイアグラムの分析 ###
            ts_idx,(r_use_newton,f_use_newton,n) = find_ts_idx(self.vimages,
                                         dif=self.irc_dif,
                                         calc_func=self.calc_func)
            if n == 0 and not self.calc_notop_ts: #TSが見えない(極大値がない)時,
                self.debug_log(f"ピークが見えない．終了")
                self.write_pt(self.ts)
                continue # TSを確認できない場合,Sellaを行なうように今後したい
            else:
                self.debug_log(f"TSがほぼ見えないが計算を続行")
                self.write_ts(self.ts)
                self.atoms_dict[f"TS{self.iter_count}"] = self.ts.copy()

            ### IRC計算 ###
            self.debug_log(f"IRC計算")
            r_converged,f_converged = self.run_irc(ts_idx, r_use_newton, f_use_newton)
            if r_converged:
                self.write_eq(self.ini)
                self.irc_ini_eq_n = self.eq_count-1 # IRC計算の結果ini構造となったEQ番号 
                self.atoms_dict[f"EQ{self.irc_ini_eq_n}"] = self.ini.copy()    
                self.debug_log(f"Revers収束,EQ{self.irc_ini_eq_n}として保存")
            if f_converged:
                self.write_eq(self.fin)
                self.irc_fin_eq_n = self.eq_count-1 # IRC計算の結果fin構造となったEQ番号
                self.atoms_dict[f"EQ{self.irc_fin_eq_n}"] = self.fin.copy() 
                self.debug_log(f"Forward収束,EQ{self.irc_fin_eq_n}として保存")
            if not all([r_converged,f_converged]):
                self.debug_log(f"IRCの一方が収束しなかったので終了")
                continue
            
            ### CONNECTION情報の分析と書き込み ###
            text = self.analyze_connection(self.ts_count-1, self.eq_count-2, self.eq_count-1, self.ts, self.ini, self.fin)
            self.write_connections(self.ts_info_file,text)

            ### 同一構造か判定とQUEへのイメージの追加 ###
            # self.sneb.images[0]         # SNEBのiniのAtoms
            # sneb_ini_idx                # SNEBのini番号
            # self.sneb.images[-1]        # SNEBのiniのAtoms
            # sneb_fin_idx                # SNEBのfin番号
            # self.ini                    # IRCのiniのAtoms
            irc_ini_idx = self.eq_count-2 # IRCのini番号
            # self.fin                    # IRCのiniのAtoms
            irc_fin_idx = self.eq_count-1 # IRCのfii番号
            
            self.debug_log(f"SNEB:{sneb_ini_idx}-{sneb_fin_idx}"+
                           f"-->IRC:{irc_ini_idx}-{irc_fin_idx}")
            self.sneb_ini.calc = self.calc_func()
            self.sneb_fin.calc = self.calc_func()
            check1,val1_1,val1_2 = self.check_structures(self.ini, self.fin) #IRC_ini,IRC_finが同じか
            check2,val2_1,val2_2 = self.check_structures(self.sneb_ini, self.ini) #SNEB_ini,IRC_iniが同じか
            check3,val3_1,val3_2 = self.check_structures(self.sneb_fin, self.fin) #SNEB_fin,IRC_finが同じか
            check4,val4_1,val4_2 = self.check_structures(self.sneb_ini, self.fin) #SNEB_ini,IRC_finが同じか
            check5,val5_1,val5_2 = self.check_structures(self.ini, self.sneb_fin) #IRC_ini,NEB_finが同じか
            
            insert_idx = self.order.index(f"EQ{sneb_ini_idx}")
            if check1 and not check2 and not check3:
                self.que.append([sneb_ini_idx,irc_ini_idx])
                self.que.append([irc_fin_idx,sneb_fin_idx])
                self.order[insert_idx+1:0] = [f"EQ{irc_ini_idx}",f"TS{self.iter_count}",f"EQ{irc_fin_idx}"]
                self.debug_log(0)
                continue
            elif not check1 and check2 and not check3:
                self.que.append([irc_fin_idx,sneb_fin_idx])
                self.order[insert_idx+1:0] = [f"EQ{irc_ini_idx}",f"TS{self.iter_count}",f"EQ{irc_fin_idx}"]
                self.debug_log(1)
                continue
            elif not check1 and not check2 and check3:
                self.que.append([sneb_ini_idx,irc_ini_idx])
                self.order[insert_idx+1:0] = [f"EQ{irc_ini_idx}",f"TS{self.iter_count}",f"EQ{irc_fin_idx}"]
                self.debug_log(2)
                continue
            elif not check1 and check4 and not check5:
                self.que.append([irc_ini_idx,sneb_fin_idx])
                self.order[insert_idx+1:0] = [f"EQ{irc_fin_idx}",f"TS{self.iter_count}",f"EQ{irc_ini_idx}"]
                self.debug_log(3)
                continue
            elif not check1 and not check4 and check5:
                self.que.append([sneb_ini_idx,irc_fin_idx])
                self.order[insert_idx+1:0] = [f"EQ{irc_fin_idx}",f"TS{self.iter_count}",f"EQ{irc_ini_idx}"]
                self.debug_log(4)
                continue
            elif not check1 and check2 and check3:
                self.debug_log(5)
                self.order[insert_idx+1:0] = [f"EQ{irc_ini_idx}",f"TS{self.iter_count}",f"EQ{irc_fin_idx}"]
                continue
            elif not check1 and check4 and check5:
                self.debug_log(6)
                self.order[insert_idx+1:0] = [f"EQ{irc_fin_idx}",f"TS{self.iter_count}",f"EQ{irc_ini_idx}"]
                continue
            else:
                if [val2_1,val3_1].count(True) > [val4_1,val5_1].count(True):
                    self.que.append([sneb_ini_idx,irc_ini_idx])
                    self.que.append([irc_fin_idx,sneb_fin_idx])
                    self.order[insert_idx+1:0] = [f"EQ{irc_ini_idx}",f"TS{self.iter_count}",f"EQ{irc_fin_idx}"]
                    self.debug_log(7)
                    continue
                elif [val2_1,val3_1].count(True) < [val4_1,val5_1].count(True):
                    self.que.append([sneb_ini_idx,irc_fin_idx])
                    self.que.append([irc_ini_idx,sneb_fin_idx])
                    self.order[insert_idx+1:0] = [f"EQ{irc_fin_idx}",f"TS{self.iter_count}",f"EQ{irc_ini_idx}"] 
                    self.debug_log(8)
                    continue
                else:
                    rmse_list = [val2_2,val3_2,val4_2,val5_2]
                    idx = rmse_list.index(min(rmse_list))
                    if idx==0 or idx==1:
                        self.que.append([sneb_ini_idx,irc_ini_idx])
                        self.que.append([irc_fin_idx,sneb_fin_idx])
                        self.order[insert_idx+1:0] = [f"EQ{irc_ini_idx}",f"TS{self.iter_count}",f"EQ{irc_fin_idx}"] 
                        self.debug_log(9)
                        continue
                    else:
                        self.que.append([sneb_ini_idx,irc_fin_idx])
                        self.que.append([irc_ini_idx,sneb_fin_idx])
                        self.order[insert_idx+1:0] = [f"EQ{irc_fin_idx}",f"TS{self.iter_count}",f"EQ{irc_ini_idx}"] 
                        self.debug_log(10)
                        continue  
        self.create_path()       
        return True if len(self.que) == 0 else False      