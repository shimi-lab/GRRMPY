from pathlib import Path
from ase.io import write
from ase.optimize import LBFGS
from tqdm.notebook import tqdm_notebook as tqdm

# USER
from grrmpy.io import log2atoms
from grrmpy.optimize.attach import automate_maxstep
from grrmpy.calculator import pfp_calculator
try:
    from grrmpy.optimize import FIRELBFGS
except:
    pass

class AutoOpt():
    """最適化後の構造は'Structure'フォルダ内にtrajファイルで保存される.
    
        計算後の構造を一括で読み込むには
        
        >>> import grrmpy.io import read_traj
        >>> atoms_list = read_traj('Structure')

        Parameters:
        
        atomslist: list of Atoms
            Atomsのリスト
        optimizer: object
            使用するOptimizer.デフォルトはLBFGS.
        constraints: ASE constraint
            | FixAtoms等の制約.
            | 複数設定する場合はリストで与える.
            | eq_list中のAtomsに既にconstraintがある場合,改めて設定する必要はない.
        trajectory: bool
            | Trueの場合,最適化途中の構造をtrajに保存する.
            | 'trajectory'フォルダー内に保存される.
        logfile: bool
            | Trueの場合, logファイルを保存する.
            | 'log'フォルダー内に保存される.
        calc_func: object
            calculatorを返す関数
    """
    def __init__(self, 
                 atomslist,
                 optimizer = LBFGS,
                 constraints = [],
                 trajectory = False,
                 logfile = True,
                 calc_func = pfp_calculator,
                 errorfile = "ERROR",
                 traj_foldername = "trajectory",
                 log_foldername = "log",
                 save_foldername = "Structures"):
        """
        
        最適化後の構造は'Structure'フォルダ内にtrajファイルで保存される.

        Parameters:
        
        atomslist: list of Atoms
            Atomsのリスト
        optimizer: object
            使用するOptimizer.デフォルトはLBFGS.
        constraints: ASE constraint
            FixAtoms等の制約.
            複数設定する場合はリストで与える.
            eq_list中のAtomsに既にconstraintがある場合,改めて設定する必要はない.
        trajectory: bool
            Trueの場合,最適化途中の構造をtrajに保存する.
            'trajectory'フォルダー内に保存される.
        logfile: bool
            Trueの場合, logファイルを保存する.
            'log'フォルダー内に保存される.
        calc_func: object
            calculatorを返す関数
        """
        self.optimizer = optimizer
        self.trajectory = trajectory
        self.logfile = logfile
        self.maxstep_dict = None
        self.attach_list = []
        
        self.atomslist = atomslist
        for atoms in self.atomslist:
            atoms.set_constraint(constraints)
            atoms.calc = calc_func()
            
        # フォルダ名,ファイル名
        self.errorfile = f"{errorfile}_{id(self)}"
        self.log_foldername = log_foldername
        self.traj_foldername = traj_foldername
        self.save_foldername = save_foldername
        
        # フォルダの作成
        self.make_folder(self.save_foldername)
        if self.trajectory:
            self.make_folder(self.traj_foldername)
        if self.logfile:
            self.make_folder(self.log_foldername)

    def make_folder(self,foldername):
        p = Path(foldername)
        if not p.exists():
            # フォルダが存在しなければ作成
            p.mkdir()
        
    def set_maxstep(self,maxstep):
        if type(maxstep) == list:
            self.maxstep = maxstep
        else:
            self.maxstep = [maxstep]
            
    def set_steps(self,steps):
        if type(steps) == list:
            self.steps = steps
        else:
            self.steps = [steps]
            
    def set_automaxstep(self,maxstep_dict):
        """auto_maxstepsを用いる場合のパラメータを変更する
        
        Examples:
        
            >>> obj.set_automaxstep({10:0.1, 5:0.2, 2:0.3, 0:0.35})
            
        必ず0のキーを含める必要があるので注意する.
        """
        self.maxstep_dict = maxstep_dict
            
    def check_param(self):
        if len(self.maxstep) != len(self.steps):
            raise Exception("maxstepとstepsの要素数が一致しません")
        if self.restart is not None and self.indexes is not None:
            raise Exception("restartとindexesは同時に指定できません")
    
        
    def errorlog(self,massage):
        with open(self.errorfile,"a") as f:
            f.write(massage)
            f.write("\n")
            
    def attach(self,*args,**kwargs):
        self.attach_list.append((args,kwargs))
        
    def irun(self,atoms,name:int,optimizer,maxstep_list,steps_list,fmax,**kwargs):
        logfile = f"{self.log_foldername}/{name}.log" if self.logfile else None
        trajectory = f"{self.traj_foldername}/{name}.traj" if self.trajectory else None
        savefile = f"{self.save_foldername}/{name}.traj"
        try:
            for maxstep,steps in zip(maxstep_list,steps_list):
                ms = 0.2 if maxstep is None else maxstep
                if optimizer == FIRELBFGS:
                    opt = FIRELBFGS(atoms,logfile=logfile,maxstep_fire=ms,maxstep_lbfgs=ms,**kwargs)
                else:
                    opt = optimizer(atoms,maxstep=ms,logfile=logfile,trajectory=trajectory,**kwargs)
                if maxstep is None:
                    opt.attach(lambda:automate_maxstep(opt,self.maxstep_dict))
                for args,kwargs in self.attach_list:
                    opt.attach(*args,**kwargs)
                opt.run(fmax=fmax,steps=steps)
            if opt.converged:
                write(savefile,atoms)
                return True
            else:
                self.errorlog(f"{name}の計算:未収束")
        except Exception as e:
            self.errorlog(f"{name}の計算:\n{e}")
        return False
           
    def run(self,
            maxstep_list=[0.05,0.2],
            steps_list=[200,10000],
            fmax=0.001,
            indexes = None,
            restart = None,
            **kwargs):
        """
        
        Parameters:
        
        maxstep_list: float or list of float
            | maxstep.
            | optimizeをFIRELBFGSにした場合,maxstep_fire,maxstep_lbfgsどちらもmaxstepで指定した値になる.
        steps_list: int or list of int
            steps
        fmax: float
            収束条件
        indexes: list of int
            | index番号で指定した構造のみを計算する
            | indexesとrestartは同時に指定できない.
        restart: int
            | 指定した番号から計算を開始する.
            | indexesとrestartは同時に指定できない.
        kwargs:
            その他,optimizer(LBFGSなど)ののインスタンス引数を指定できる.
        """
        self.set_maxstep(maxstep_list)
        self.set_steps(steps_list)
        self.indexes = indexes
        self.restart = restart
        self.check_param()
        
        if self.indexes is None:
            self.indexes = [i for i in range(len(self.atomslist))]
        if self.restart is not None:
            self.indexes = [i for i in range(self.restart,len(self.atomslist))]
        for i in tqdm(self.indexes):
            atoms = self.atomslist[i]
            self.irun(atoms,i,self.optimizer,maxstep_list,steps_list,fmax,**kwargs)
            