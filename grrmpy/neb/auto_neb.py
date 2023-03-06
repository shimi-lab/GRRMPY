from ase.optimize import FIRE,LBFGS
from ase.neb import NEB, interpolate
from ase.geometry.analysis import Analysis
from ase.geometry import find_mic
from ase.units import kJ,mol
from ase import Atoms
import numpy as np
from scipy import signal

from math import sqrt
import warnings
from functools import reduce
from operator import mul

# User
from grrmpy.calculator import pfp_calculator
from grrmpy.optimize import automate_maxstep
from grrmpy.optimize.functions import add_stop_to
from grrmpy.neb.functions import to_html_nebgraph

class ANEB():
    def __init__(self,
                 *data,
                 logfile="-",
                 html=None,
                 mic=None,
                 parallel=True,
                 calc_func=pfp_calculator,
                 optimizer=FIRE,
                 with_stop=True,
                 max_n=12,
                 times=2,
                 constraints=[]):
        """Adaptive NEBを行なう.

        Parameters:
        
        *data: 
            | 3種類の方法でinput構造を与えることができる.
            | 
            | - ini,fin構造と初期イメージの数を与える(3つの引数を与える)
            | - ini,fin構造を与える(2つの引数を与える.初期イメージの数は13にする)
            | - images構造を与える(1つの引数を与える.calculatorを付ける必要はない)
            |       imageの数は奇数個にすることが推奨される.
        logfile: str or Path
            logファイルを出力する場合,ファイル名.
        html: str or Path
            途中のNEBイメージを保存する場合にファイル名を設定.
        mic : bool
            周期境界条件の際, 最小距離でイメージを作成する場合.True.
            設定しない場合(None). ini構造が周期境界条件である場合.自動的にTrueにする.
        parallel: bool
            並行処理を行なう場合.True. デフォルトはTrue.
        calc_func: object
            calculatorを返す関数.デフォルトはpfpのcalculator
        optimizer:
            NEB計算を行なう際のOptimizer. デフォルトはFIRE.
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
        constraints: constrain or list
            適用するconstraint
            
        Note:
            arg以外の引数は全てキーワード引数なので注意!
            
        Examples:
            >>> aneb = ANEB(ini,fin,15) # 15イメージで作成される
            >>> # または
            >>> aneb = ANEB(ini,fin) # 13イメージで作成される
            >>> # または
            >>> aneb = ANEB(images)  
        
        """
        self.attach_dict = {} # {func:interval}
        self.logfile = logfile
        self.parallel=parallel
        self.calc_func = calc_func
        self.mic = mic
        self.html = html
        self.constraints = constraints
        #: 全てのNEBイメージが2次元リストで保存されている
        self.archive = []
        if with_stop:
            optimizer = add_stop_to(optimizer,max_n,times)
        self.optimizer = optimizer
        
        if len(data) == 1 and type(data[0])==list:
            data_type = "ini_fin"
            #: 現在のNEBイメージ
            self.images = data[0]
        elif type(data[0])==Atoms and type(data[1])==Atoms:
            data_type = "images"
            ini = data[0]
            fin = data[1]
            nimages = 13 if len(data)==2 else data[2]
            self.images = [ini.copy()]+[ini.copy() for _ in range(nimages)]+[fin.copy()]
        else:
            raise Exception("Dataはini,fin構造またはNEB imagesである必要があります")
            
        if mic is None:
            if any(self.images[0].get_pbc()):
                self.mic = True
            else:
                self.mic = False
        
        if data_type == "images":
            interpolate(self.images,mic=self.mic,apply_constraint=False)

        for image in self.images:
            image.calc = self.calc_func()
            image.set_constraint(self.constraints)
        
        if self.html:
            with open(self.html,"w") as f:
                f.write("<meta charset='UTF-8' />")
                
    @property
    def imax(self):
        """NEB計算上でTSと判断されたindex番号"""
        return self.neb.imax
    
    def get_near_ts(self,ts,atoms,dist):
        """
        ts:ts構造
        atoms:ts構造の隣の構造
        """
        pos1 = atoms.get_positions()
        pos2 = ts.get_positions()
        d = pos2 - pos1
        if self.mic:
            d = find_mic(d, ts.get_cell(), ts.pbc)[0]
        n = np.linalg.norm(d)/dist
        d /= (n - 1.0)
        new_pos = pos1 + (n-2) * d # n-2(ts付近)
        unconstrained_image = ts.copy()
        unconstrained_image.set_positions(new_pos,apply_constraint=True)
        unconstrained_image.calc = self.calc_func()
        return unconstrained_image
        
    def updata_images(self,nsampling,nimages,i=None,dist=0.1):
        """
        | TS構造から,forward側,reverse側にそれぞれnsampling個隣の構造を
        | 新たなini,fin構造としてイメージを作成する.
        | ini,finの間にnimages個のイメージを作成する.
        | iが指定されていない場合,neb.imaxをTS構造とする.
        | TS構造からdist[Å]離れた位置にイメージを作成する
        
        Note:
            | 
        """
        if i is None:
            i = self.neb.imax
        ini_idx = i - nsampling
        fin_idx = i + nsampling
        if i == 0:
            near_ts = self.get_near_ts(self.images[i],self.images[i+1],dist)
            fin_atoms = self.images[nsampling*2]
            imgs = [near_ts.copy()]+[near_ts.copy() for _ in range(nimages-2)]+[fin_atoms.copy()]
            interpolate(imgs,mic=self.mic,apply_constraint=False)
            self.images = [self.images[0].copy()] + imgs
        elif i ==  len(self.images)-1:
            near_ts = self.get_near_ts(self.images[i],self.images[i-1],dist)
            ini_atoms = self.images[-(nsampling*2+1)]
            imgs = [ini_atoms.copy()]+[ini_atoms.copy() for _ in range(nimages-2)]+[near_ts.copy()]
            interpolate(imgs,mic=self.mic,apply_constraint=False)
            self.images = imgs + [self.images[-1].copy()]
        else:
            if ini_idx < 0:
                ini_idx = 0
                fin_idx = nsampling*2
            elif fin_idx >  len(self.images)-1:
                fin_idx = -1
                ini_idx = -(nsampling*2+1)
            ini_near_ts = self.get_near_ts(self.images[i],self.images[i-1],dist)
            fin_near_ts = self.get_near_ts(self.images[i],self.images[i+1],dist)
            if nimages%2 == 0:
                n = int((nimages-2)/2)
                n = (n,n-1)
            else:
                n = int((nimages-3)/2)
                n = (n,n)
            imgs1 = [self.images[ini_idx].copy()]+[self.images[0].copy() for _ in range(n[0])]+[ini_near_ts.copy()]
            interpolate(imgs1,mic=self.mic,apply_constraint=False)
            img2 = [fin_near_ts.copy()]+[fin_near_ts.copy() for _ in range(n[1])]+[self.images[fin_idx].copy()]
            interpolate(img2,mic=self.mic,apply_constraint=False)
            self.images = imgs1 + [self.images[i].copy()]+img2
        for image in self.images:
            image.calc = self.calc_func()
            image.set_constraint(self.constraints)
        return ini_idx, fin_idx
            
    def make_neb(self,climb:bool):
        self.neb = NEB(self.images,climb=climb,parallel=self.parallel)
        
    def make_opt(self,maxstep:float=None):
        if maxstep is None:
            self.opt = self.optimizer(self.neb, logfile=self.logfile)
            self.opt.attach(lambda:automate_maxstep(self.opt),interval=10)
        else:
            self.opt = self.optimizer(self.neb,maxstep=maxstep, logfile=self.logfile)
        
    def attach(self,func,interval:int=1):
        self.attach_dict[func]=interval
        
    def write_html(self,comment,**kargs):
        if self.html is None:
            return
        with open(self.html,"a") as f:
            html_text = to_html_nebgraph(self.neb,self.calc_func,False,include_plotlyjs="cdn",**kargs)
            f.write(comment)
            f.write(html_text)
        
    def iter_run(self,fmax:float,steps:int,climb=True,maxstep:float=None):
        self.make_neb(climb)
        self.make_opt(maxstep)
        for func,i in self.attach_dict.items():
            self.opt.attach(func,interval=i)
        converged = self.opt.run(fmax=fmax,steps=steps)    
        self.archive.append([image.copy() for image in self.images])
        return converged
    
    def run(self,
            fmax=[0.12, 0.1, 0.07, 0.05],
            steps=[500, 1000, 1000,2000],
            maxstep=[0.2, 0.2, 0.1, None],
            nimages=[13, 9, 12],
            nsampling=[3, 3, 3],
            dist=0.1):
        """
        
        | ini,fin構造の間に3つのイメージを作成し,緩い収束条件でNEB計算を行なう.
        | 計算後のTS構造の両隣の構造を新たなini,fin構造としてNEB計算を行なう.
        | この操作を合計len(nimages)回繰り返す(NEB計算をlen(nimages)回行なう).
        | 最後のNEB計算のみ収束条件を厳しく(=fmax)する.
        | デフォルトではn=4回のNEB計算を行なう.
        
        Parameters:
        
        fmax: float or list of float
            | floatで与えた場合,[fmax*2,fmax*2,fmax*2,...fmax]の収束条件で計算する.
            | (最後だけ厳しい収束条件で計算する)
            | listで与えた場合, [fmax_0, fmax_1, fmax_2,...fmax_n]の収束条件で計算する.
        steps: integers or list of integers
            | integersで与えた場合,[steps,steps,steps,...steps*5]のステップ数で計算する.
            | listで与えた場合, [steps_0, steps_1, steps_2,...steps_n]のステップ数で計算する.
        maxstep: float or list of float
            | maxstep. Noneの場合,Forceに合わせて自動的に調整する.
            | 詳しくはgrrmpy.optimize.attach.automate_maxsteps()を参照
            | floatで与えた場合,maxstepの値で計算する.
            | listで与えた場合, [maxstep_0, maxstep_1, maxstep_2,...maxstep_n]の収束条件で計算する.
            | None,またはNoneの要素を与えた場合,forceに合わせて自動的に調整する.
        nimages: list of integers
            | 毎回のNEB計算の中間イメージ数.
            | [nimages_0, nimages_1, nimages_2,...nimages_n]イメージ数で作成する.
        nsampling: integers or list of integers
            | TS構造からforward側,reverse側にそれぞれnsampling個隣の構造を
            | 新たなini,fin構造としてイメージを作成する.
            | integersで与えた場合,nsamplingのイメージ数で作成する.
            | listで与えた場合, [nsampling_0, nsampling_1, nsampling_2,...nsampling_n]イメージ数で作成する.
            | Noneの場合,自動的に調整する.
        dist: float
            | 新たなNEBイメージを作成する際に,TS付近に2つイメージを作成する.
            | TS構造からの変位をdistで設定する
            
        Retrun:
            bool: 収束した場合True.
        """
        ###引数の検証####
        if type(nimages) != list:
            raise TypeError("nimagesはリストで与えて下さい.")
        
        if type(fmax) == list:
            if len(fmax) != len(nimages)+1:
                raise IndexError("fmaxはflaotまたはnimagesと同じ要素数のflaotのリストで与えください.")
        else:
            fmax = [fmax*2 for _ in range(len(fmax)-1)]+[fmax]
            
        if type(steps) == list:
            if len(steps) != len(nimages)+1:
                raise IndexError("stepsはintまたはnimagesと同じ要素数のintのリストで与えください.")
        else:
            steps = [steps for _ in range(len(fmax)-1)]+[steps*5]
            
        if type(maxstep) == list:
            if len(maxstep) != len(nimages)+1:
                raise IndexError("maxstepはflaotまたはnimagesと同じ要素数のflaotのリストで与えください.")
        else:
            maxstep = [maxstep for _ in range(len(fmax))]
            
        if nsampling is None:
            nsampling = [int((n-2)/2) if n%2 == 0 else int((n-1)/2) for n in nsampling[:-1]]
        elif type(nsampling) == list: 
            if len(nsampling) != len(nimages):
                raise IndexError("nsamplingはintまたは(nimagesの要素数-1)の要素数のintのリストで与えください.")
        else:
            nsampling = [nsampling for _ in range(len(fmax)-1)]
        for s,n in zip(nsampling,[len(self.images)-2]+nimages[:-1]):
            if s*2+1 > n+2:
                raise("nsamplingに不適切な値が含まれています,nsamplingが大きすぎる可能性があります")
        ################
        ###何パーセントzoomするか出力###
        zooms = [(2*j)/(i+1) for i,j,in zip([len(self.images)-2]+nimages[:-1],nsampling)] 
        zoom = reduce(mul, zooms)*100
        if self.logfile == "-":
            print(f"zoom:{zoom}%")
        else:
            with open(self.logfile,"a") as f:
                f.write(f"zoom:{zoom}%\n")
            
        try:
            for f,s,m,n,sp in zip(fmax[:-1],steps[:-1],maxstep[:-1],nimages,nsampling):
                c = self.iter_run(fmax=f,steps=s,climb=False,maxstep=m)
                self.write_html(f"<p>climb=False, converged={c}</p>")
                s = s if c else 15
                c = self.iter_run(fmax=f,steps=s,climb=True,maxstep=m)
                ini_idx,fin_idx = self.updata_images(sp,n,dist=dist)
                self.write_html(f"<p>climb=True, converged={c}</p><p>{ini_idx}番,{fin_idx}番をini,finに選択</p>")
            # 最後のNEB計算
            c = self.iter_run(fmax=fmax[-1],steps=steps[-1],climb=False,maxstep=maxstep[-1])
            self.write_html(f"<p>climb=False, converged={c}</p>")
            c = self.iter_run(fmax=fmax[-1],steps=steps[-1],climb=True,maxstep=maxstep[-1])
            self.write_html(f"<p>climb=True, converged={c}</p>")
            return c
        except Exception as e:
            raise Exception(e)

class SNEB(ANEB):
    """
        Parameters:
        
        \*data: 
            | 3種類の方法でinput構造を与えることができる.
            | 
            | - ini,fin構造と初めの中間イメージの数を与える(3つの引数を与える)
            | - ini,fin構造を与える(2つの引数を与える.初期イメージの数は13になる)
            | - images構造を与える(1つの引数を与える.calculatorを付ける必要はない)
            |       imageの数は奇数個を推奨する.
            
            Examples:
            
                インスタンス化方法
            
                >>> sneb = SNEB(ini,fin,html="NEB_progress.html") # 13個で中間イメージが作られる
                >>> # または
                >>> sneb = SNEB(ini,fin,15,html="NEB_progress.html") # 15個で中間イメージが作られる
                >>> # または
                >>> sneb = SNEB(images,html="NEB_progress.html")
                
                Constraintsがある場合
                
                >>> f = FixAtoms(indices=[0,1,2,3])
                >>> sneb = SNEB(ini,fin,html="NEB_progress.html",constraints=f)
                
                >>> f = FixAtoms(indices=[0,1,2,3])
                >>> c = FixBondLength(0, 1)
                >>> sneb = SNEB(ini,fin,html="NEB_progress.html",constraints=[f,c])
                
        logfile: str or Path
            logファイルを出力する場合,ファイル名.
        html: str or Path
            途中のNEBイメージを保存する場合にファイル名を設定.
        mic : bool
            周期境界条件の際, 最小距離でイメージを作成する場合.True.
            設定しない場合(None). ini構造が周期境界条件である場合.自動的にTrueにする.
        parallel: bool
            並行処理を行なう場合.True. デフォルトはTrue.
        calc_func: object
            calculatorを返す関数.デフォルトはpfpのcalculator
        optimizer: class
            NEB計算を行なう際のOptimizer. デフォルトはFIRE.
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
            
        Note:
            \*data以外の引数は全てキーワード引数になるので注意!!
    """
    def __init__(self,
                 *data,
                 logfile="-",
                 html=None,
                 mic=None,
                 parallel=True,
                 calc_func=pfp_calculator,
                 optimizer=FIRE,
                 with_stop=True,
                 max_n=12,
                 times=2,
                 constraints=[]):
        """Separative NEB
        
        | 緩い収束条件でNEB計算を行ない,最もエネルギーの高い点をTSとする.
        | TS構造の両隣の構造のエネルギーを順次確認し,エネルギーの谷間があった場合,
        | その点を新たなini,fin構造としてNEB計算を行なう.
        | これらの一連の計算を何度か行ない,最後に厳しい収束条件でNEB計算を行なう.
        | 
        | 浅い谷間は無視するため, toleranceの引数によりtolerance kJ/mol以下の谷間は無視する.
        | NEB計算ではTS近くに構造を密に配置する事で無駄な計算量を減らし計算時間を短縮できる.
        | NEBバンドの裾野が広い場合,計算時間が長くなるため,thresholdの引数により,
        | 山の高さのthreshold%以内までの点は谷間でなくても切断する.
        
        Parameters:
        
        *data: 
            | 3種類の方法でinput構造を与えることができる.
            | 
            | - ini,fin構造と初めの中間イメージの数を与える(3つの引数を与える)
            | - ini,fin構造を与える(2つの引数を与える.初期イメージの数は13になる)
            | - images構造を与える(1つの引数を与える.calculatorを付ける必要はない)
            |       imageの数は奇数個を推奨する.
        logfile: str or Path
            logファイルを出力する場合,ファイル名.
        html: str or Path
            途中のNEBイメージを保存する場合にファイル名を設定.
        mic : bool
            周期境界条件の際, 最小距離でイメージを作成する場合.True.
            設定しない場合(None). ini構造が周期境界条件である場合.自動的にTrueにする.
        parallel: bool
            並行処理を行なう場合.True. デフォルトはTrue.
        calc_func: object
            calculatorを返す関数.デフォルトはpfpのcalculator
        optimizer: class
            NEB計算を行なう際のOptimizer. デフォルトはFIRE.
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
        constraints: constraint obj or list
            適用するconstraintのオブジェクトまたはそのリスト
            
        Note:
            *data以外の引数は全てキーワード引数になるので注意!!
            
        Examples:
        
            >>> sneb = SNEB(ini,fin,html="NEB_progress.html") # 13個で中間イメージが作られる
            >>> # または
            >>> sneb = SNEB(ini,fin,15,html="NEB_progress.html") # 15個で中間イメージが作られる
            >>> # または
            >>> sneb = SNEB(images,html="NEB_progress.html")

        """
        super().__init__(*data,
                         logfile=logfile,
                         html=html,
                         mic=mic,
                         parallel=parallel,
                         calc_func=calc_func,
                         optimizer=optimizer,
                         with_stop=with_stop,
                         max_n=max_n,
                         times=times,
                         constraints=constraints)
        
    def updata_images(self,nimages,tolerance,threshold,dist,min_nimages,i=None):
        if i is None:
            imax = self.neb.imax
        else:
            imax = i
        energies = np.array([i.get_potential_energy()*mol/kJ for i in self.images])
        ts_e = energies[imax]

        max_idx = list(signal.argrelmax(np.array(energies))[0]) # 極大値のindexを取得
        min_idx = list(signal.argrelmin(np.array(energies))[0]) # 極小値のindexを取得
        
        m = int((min_nimages-1)/2) # tsの±m隣は必ず含める
        forward_max_idx = list(reversed(list(filter(lambda x: x<imax-m, max_idx))))        
        reverse_max_idx = list(filter(lambda x: x>imax+m, max_idx))
        forward_min_idx = list(reversed(list(filter(lambda x: x<imax-m, min_idx))))
        reverse_min_idx = list(filter(lambda x: x>imax+m, min_idx))
        
        # iniのindex番号を探す(toleranceを考慮しつつ)
        ini_idx = 0
        if len(forward_max_idx) >1 and len(forward_min_idx)>0:
            if forward_max_idx[0] > forward_min_idx[0]:
                """iがneb.maxでなく,任意の値にしている時,iが極大値でない可能性があるため"""
                forward_max_idx = forward_max_idx[1:]
        for ma,mi in zip(forward_max_idx,forward_min_idx):
            """!!!zip中のforward_max_idx,forward_min_idxは要素数が異なる場合もある!!!"""
            if energies[ma]-energies[mi] <= tolerance:
                continue
            else:
                ini_idx = mi
                break
            
        # finのindex番号を探す(toleranceを考慮しつつ)
        fin_idx = len(self.images)-1
        if len(reverse_max_idx) >1 and len(reverse_min_idx)>0:
            if reverse_max_idx[0] < reverse_min_idx[0]:
                """iがneb.maxでなく,任意の値にしている時,iが極大値でない可能性があるため"""
                reverse_max_idx = reverse_max_idx[1:]
        for ma,mi in zip(reverse_max_idx,reverse_min_idx):
            if energies[ma]-energies[mi] <= tolerance:
                continue
            else:
                fin_idx = mi
                break
            
        ini_H = ts_e - energies[ini_idx] # 正反応の活性化エネルギー(もちろん正確ではないが)
        fin_H = ts_e - energies[fin_idx] # 逆反応の活性化エネルギー(もちろん正確ではないが)
        
        # 再度thresholdを考慮し,iniのindex番号を探す
        ini_idx_copy = ini_idx
        for i in range(ini_idx_copy,imax):
            if energies[i]-energies[ini_idx_copy] < ini_H*threshold*0.01:
                """正反応の活性化エネルギーのthreshold%以下の点があった場合それを新たなini_idxとする"""
                ini_idx = i
                
        # 再度thresholdを考慮し,finのindex番号を探す
        fin_idx_copy = fin_idx
        for i in reversed(range(imax, fin_idx_copy+1)):
            if energies[i]-energies[fin_idx_copy] < fin_H*threshold*0.01:
                """正反応の活性化エネルギーのthreshold%以下の点があった場合それを新たなini_idxとする"""
                fin_idx = i 
                
        if fin_idx-ini_idx+1 < min_nimages:
            n = int((min_nimages-1)/2)
            ini_idx = imax-n
            fin_idx = imax+n
            if ini_idx < 0:
                ini_idx = 0
                fin_idx = min_nimages-1
            if fin_idx > len(self.images)-1:
                fin_idx = len(self.images)-1
                ini_idx = len(self.images)-min_nimages        
        
        ts_images = self.images[imax]
        if ini_idx==0 and fin_idx==len(self.images)-1:
            pass
        elif ini_idx == imax:
            fin_near_ts = self.get_near_ts(self.images[imax],self.images[imax+1],dist)
            fin_images = [fin_near_ts.copy()]+[fin_near_ts.copy() for _ in range(nimages-1)]+[self.images[fin_idx].copy()]
            interpolate(fin_images,mic=self.mic,apply_constraint=False)
            self.images = [ts_images.copy()]+fin_images
        elif fin_idx == imax:
            ini_near_ts = self.get_near_ts(self.images[imax],self.images[imax-1],dist)
            ini_images = [self.images[ini_idx].copy()]+[self.images[ini_idx].copy() for _ in range(nimages-1)]+[ini_near_ts.copy()]
            interpolate(ini_images,mic=self.mic,apply_constraint=False)
            self.images = ini_images + [ts_images.copy()]
        else:
            ini_near_ts = self.get_near_ts(self.images[imax],self.images[imax-1],dist)
            fin_near_ts = self.get_near_ts(self.images[imax],self.images[imax+1],dist)
            ini_images = [self.images[ini_idx].copy()]+[self.images[ini_idx].copy() for _ in range(int((nimages-3)/2))]+[ini_near_ts.copy()]
            interpolate(ini_images,mic=self.mic,apply_constraint=False)
            fin_images = [fin_near_ts.copy()]+[fin_near_ts.copy() for _ in range(int((nimages-3)/2))]+[self.images[fin_idx].copy()]
            interpolate(fin_images,mic=self.mic,apply_constraint=False)
            self.images = ini_images + [ts_images.copy()] + fin_images
 
        for image in self.images:
            image.calc = self.calc_func()
            image.set_constraint(self.constraints)
        return ini_idx, fin_idx
            
    def run(self,
            nimages=[13, 9, 12],
            maxstep=[0.2, 0.2, 0.1, None],
            fmax=[0.12, 0.1, 0.07, 0.05],
            steps=[500, 1000, 1000,2000],
            tolerance=5,
            threshold=[30,20,20],
            dist=0.1,
            min_nimages=3,
            climb_steps=15,
            climb = False,):
        """
        Parameters:
        
        nimages: list of integers
            毎回のNEB計算のイメージの数, 5以上の奇数.
        maxstep: list of float
            毎回のNEB計算のmaxstep値,Noneを要素にとった場合,fmaxに応じて自動調整する.
        fmax: list of float
            毎回のNEB計算のfmax値
        steps: list of integers
            毎回のNEB計算のsteps値
        tolerance: float
            tolerance(kJ/mol)以下の谷間は無視する.デフォルトは5 kJ/mol
        threshold: list of float
            山の高さのthreshold%以下の点は切り捨てる.デフォルトは30%
        dist: float
            | 新たなNEBイメージを作成する際に,TS付近に2つイメージを作成する.
            | TS構造からの変位をdistで設定する
        min_nimages: integer
            | NEBバンドを切断する際,最低min_nimages個の点を残す.デフォルトは3
            | (一個だけ突出しているTS(Wみたいな形)のNEBバンドがあった時のために設定する)
        climb_steps: integer
            | climb=Falseで収束しなかった場合,その後のTrueでのNEB計算のStep数をclimb_steps回にする.
            | Noneの場合,climb=Falseの時と同じstepsで計算する.
            | 一番最後のNEB計算ではclimb_stepsに関係なくclimb=Falseの時と同じstepsで計算する.
        climb: bool
            | Falseの場合,climbはFalse,False,False,False...Trueで行なう.
            | Trueの場合はclimbはFalse,True,Fase,True,False....Trueで行なう
            | Falseにした場合,climb_stepsを設定する意味はなくなる.
        """
        ####引数の検証######
        if any(list(map(lambda x:x<=5,nimages))):
            raise ValueError(f"nimagesはmin_nimages+2={min_nimages+2}以上を設定してください.")
        if min_nimages%2 == 0:
            raise ValueError("min_nimagesは奇数且つ,nimages以下の数を設定してください.")
        arg_check1 = len(nimages)+1!=len(maxstep)
        arg_check2 = len(nimages)+1!=len(fmax)
        arg_check3 = len(nimages)+1!=len(steps)
        arg_check4 = len(nimages)!=len(threshold)
        if any([arg_check1,arg_check2,arg_check3,arg_check4]):
            raise IndexError("maxstep,fmax,steps,threshold+1はthresholdと同じ要素数である必要があります.")
        ####################
        self.climb_steps = climb_steps
        self.clibm = climb
        
        try:
            for f,s,m,n,t in zip(fmax[:-1],steps[:-1],maxstep[:-1],nimages,threshold):
                converged = self.iter_run(fmax=f,steps=s,climb=False,maxstep=m)
                if self.clibm: # self.climb=Trueの場合climb=Trueでも計算する
                    self.write_html(f"<p>climb={self.neb.climb}, converged={converged}</p>")
                    if not converged:
                        s = s if self.climb_steps is None else self.climb_steps
                    converged = self.iter_run(fmax=f,steps=s,climb=True,maxstep=m)
                ini_idx,fin_idx = self.updata_images(n,tolerance,t,dist,min_nimages)
                self.write_html(f"<p>climb={self.neb.climb}, converged={converged}</p><p>{ini_idx}番,{fin_idx}番をini,finに選択</p>")
            # 最後のNEB計算
            converged = self.iter_run(fmax=fmax[-1],steps=steps[-1],climb=False,maxstep=maxstep[-1])
            self.write_html(f"<p>climb={self.neb.climb}, converged={converged}</p>")
            converged =  self.iter_run(fmax=fmax[-1],steps=steps[-1],climb=True,maxstep=maxstep[-1])
            self.write_html(f"<p>climb={self.neb.climb}, converged={converged}</p>")
            return converged
        except Exception as e:
            raise Exception(e)

