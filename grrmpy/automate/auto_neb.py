from ase.optimize import FIRE
from ase.io import write
from ase import Atoms
from ase.neb import NEB,interpolate
from pathlib import Path
import traceback
from tqdm.notebook import tqdm_notebook as tqdm

# USER
from grrmpy.calculator import pfp_calculator
from grrmpy.functions import set_calc2images,set_const2images,copy_images,get_unique_connections
from grrmpy.neb.repetitive_neb import RepetitiveNEB
from grrmpy.io import write_html

class AutoNEB():
    """連続的にNEB計算を行なう
    
    | attach(),attach_stop(),update(),next_i()に関数を与える事ができる.
    | attach()とattach_stop()はイタレーション毎に実行される.
    | update()はNEB計算が終了する毎に実行される.
    | next_i()は1サンプルの計算を終える毎に実行される.
    
    | STOP.txtファイルをディレクトリ中に設置すると計算を一時停止,またはスキップできる.
    | STOP.txtの中に
    | STOP = TUREとすると,計算が一時停止される.
    | SKIP = TRUEとすると現在実行中サンプルの計算が中止される.
    
    | STOP.txtで停止を行なうとRESTART.txtが作成される.
    | ディレクトリ中にRESTART.txtが存在する場合,途中から再開する.
    
    Parameters:
    
    eq_list: list of Atoms
        | EQのリスト(Atoms)
    connections:
        | CONNECTiONS情報
        | Ex) [[0,1],[2,4],[3,5],...]
    constraints: Constraintオブジェクト 
        | FixAtomsなどの制約.
    parallel: bool
        | NEBで並列計算を行なう場合True.
    mic: bool
        | 最小イメージ規則を適用する場合,True.
        | デフォルトでは周期境界条件で計算する場合自動でTrueにする.
    calc_func: function object
        | calculatorを返す関数
    savefolder: str
        | 収束した構造を保存するためのフォルダ
    unconvergedfolder: str
        | 収束しなかった構造を保存するためのフォルダ
    logfolder: str
        | logファイルを格納するフォルダ
    htmlfolder:
        | エネルギーダイアグラムのhtmlを保存するためのフォルダ
    errorfile:
        | ERRORファイルのパス
        | 計算中にエラーが発生した場合,ERRORファイルに書き込まれる.
    """
    def __init__(self, eq_list, connections, constraints=[], parallel=True, mic = None, 
                 calc_func=pfp_calculator, savefolder="Structures",unconvergedfolder="unconv",
                 logfolder="log", htmlfolder="html", errorfile="ERROR"):
        """連続的にNEB計算を行なう
        
        attach(),attach_stop(),update(),next_i()に関数を与える事ができる.
        attach()とattach_stop()はイタレーション毎に実行される.
        update()はNEB計算が終了する毎に実行される.
        next_i()は1サンプルの計算を終える毎に実行される.
        
        STOP.txtファイルをディレクトリ中に設置すると計算を一時停止,またはスキップできる.
        STOP.txtの中に
        STOP = TUREとすると,計算が一時停止される.
        SKIP = TRUEとすると現在実行中サンプルの計算が中止される.
        
        STOP.txtで停止を行なうとRESTART.txtが作成される.
        ディレクトリ中にRESTART.txtが存在する場合,途中から再開する.
        
        Parameters:
        
        eq_list: list of Atoms
            EQのリスト(Atoms)
        connections:
            CONNECTiONS情報
            Ex) [[0,1],[2,4],[3,5],...]
        constraints: Constraintオブジェクト 
            FixAtomsなどの制約.
        parallel: bool
            NEBで並列計算を行なう場合True.
        mic: bool
            最小イメージ規則を適用する場合,True.
            デフォルトでは周期境界条件で計算する場合自動でTrueにする.
        calc_func: function object
            calculatorを返す関数
        savefolder: str
            収束した構造を保存するためのフォルダ
        unconvergedfolder: str
            収束しなかった構造を保存するためのフォルダ
        logfolder: str
            logファイルを格納するフォルダ
        htmlfolder:
            エネルギーダイアグラムのhtmlを保存するためのフォルダ
        errorfile:
            ERRORファイルのパス
            計算中にエラーが発生した場合,ERRORファイルに書き込まれる.
        """
        self.calc_func= calc_func
        self.constraints = constraints
        self.connections = self.gen_connections(connections)
        self.parallel = parallel
        self.eq_list = copy_images(eq_list)
        self.logfolder = logfolder
        self.errorfile = errorfile
        self.savefolder = savefolder
        self.htmlfolder = htmlfolder
        self.unconvergedfolder = unconvergedfolder
        
        if mic is None:
            mic = True if any(self.eq_list[0].get_pbc()) else False
        self.mic = mic
        self.attach_list = []
        self.attach_stop_list = []
        self.update_list = []
        self.next_i_list = []
        
    def make_images(self,ini,fin,nimage):
        images = [ini.copy() for _ in range(nimage+1)]+[fin.copy()]
        interpolate(images,mic=self.mic)
        set_calc2images(images,self.calc_func)
        set_const2images(images, self.constraints)
        return images
    
    def gen_connections(self,connections):
        new_connections,idx = get_unique_connections(connections)
        new_connections = [c if i in idx else None for i, c in enumerate(connections)]
        return new_connections
    
    def make_folder(self,foldername):
        p = Path(foldername)
        if not p.exists():
            # フォルダが存在しなければ作成
            p.mkdir()
            
    def check_stop_file(self):
        p = Path("STOP.txt")
        if p.exists():
            pass
                
    def attach(self,*args,**kwargs):
        """
        | ASEのOptimizerのattachと同じ.
        | イタレーション毎にattachで与えた関数が実行される.
        """
        self.attach_list.append((args,kwargs))
        
    def attach_stop(self,*args,**kwargs):
        """
        | grrmpy.optimize.optimizer.OptWithCondのattach_stopと同じ
        | True,Falseを返す関数を設定する.
        | イタレーション毎に関数が実行されTrueを返した時に計算が停止する
        """
        self.attach_stop_list.append((args,kwargs))
        
    def update(self,function):
        """
        | grrmpy.neb.repetitive_neb.RepetitiveNEBのuodateと同じ
        | updateに関数を与えると,NEBの計算が終わるごとに実行される
        """
        self.update_list.append(function)
        
    def next_i(self,function):
        """
        1サンプル計算が終える度に実行される関数を与えられる.
        """
        self.next_i_list.append(function)
        
    def write_errorlog(self,text):
        with open(f"{self.errorfile}_{id(self)}","a") as f:
            f.write(text)    
        
    def make_nimages_list(self,nimages):
        nimages_list = [nimages for _ in range(len(self.connections))] if type(nimages)==int else nimages
        if len(nimages_list) != len(self.connections):
            raise("nimagesの要素数がconnectionの要素数と一致しません")
        return nimages_list
    
    @property
    def default_param(self):
        """run()の際に使用されるデフォルトの引数.
        
        climb_list,maxstep_list,fmax_list,steps_listの要素数は一致している必要がある.
        
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
        return RepetitiveNEB._default_param
    
    @property
    def images(self):
        return self.rneb.images
    
    def run(self,nimages,restart=0,**param):
        """
        
        Parameters:
        
        nimages: int or list of int
            NEBのイメージの数. listで指定した場合はサンプル毎に違うイメージ数でNEB計算できる.
            
        その他の引数はdefault_paramで確認する
        
        Examples:
            
            >>> auto_neb.run(16,steps_list=[100,1000,1500])
            
            以下のコードでも同じ処理を行なえる
            
            >>> param = auto_neb.default_param
            >>> param["steps_list"] = [100,1000,1500]
            >>> auto_neb.run(16,**param)
            
            途中から計算を再開するには
            
            >>> auto_neb.run(16,restart=8) # NEB8から計算を開始する
            
        """
        self.make_folder(self.logfolder)
        self.make_folder(self.savefolder)
        self.make_folder(self.htmlfolder)
        self.make_folder(self.unconvergedfolder)
        nimages_list = self.make_nimages_list(nimages)
        nimages_list = nimages_list[restart:]
        connections = self.connections[restart:]
        for i,(nimage,connection) in enumerate(zip(tqdm(nimages_list),connections),start=restart):
            self.i = i
            self.connection = connection
            try:
                if not connection is None:
                    ini = self.eq_list[connection[0]]
                    fin = self.eq_list[connection[1]]
                    images = self.make_images(ini,fin,nimage)
                    self.rneb = RepetitiveNEB(
                        images,
                        constraints=self.constraints, 
                        parallel=self.parallel, 
                        mic=self.mic, 
                        calc_func=self.calc_func, 
                        name=f"{self.logfolder}/{i}")
                    self.rneb.attach(
                        lambda:write_html(
                            f"{self.htmlfolder}/{i}(EQ{connection[0]}-EQ{connection[1]}).html",
                            self.rneb.images))
                    for args,kwargs in self.attach_list:
                        self.rneb.attach(*args,**kwargs)
                    for args,kwargs in self.attach_stop_list:
                        self.rneb.attach_stop(*args,**kwargs)     
                    for function in self.update_list:
                        self.rneb.update(function)
                    self.converged = self.rneb.run(**param)
                    if self.converged:
                        write(f"{self.savefolder}/{i}.traj",self.rneb.images)
                    else:
                        write(f"{self.unconvergedfolder}/{i}.traj",self.rneb.images)
                    for functions in self.next_i_list:
                        functions()
            except Exception as e:
                self.write_errorlog(f"{i}番目の計算:\n{e}\n{traceback.format_exc()}\n")
