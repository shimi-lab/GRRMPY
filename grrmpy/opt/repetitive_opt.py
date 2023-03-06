from ase.optimize import LBFGS
# USER
from grrmpy.calculator import pfp_calculator
from grrmpy.optimize.attach import automate_maxstep
try:
    from grrmpy.optimize import FIRELBFGS
    OPTIMIZER = FIRELBFGS
except:
    OPTIMIZER = LBFGS

class RepetitiveOpt():
    """異なる条件で繰り返し構造最適化を行なう
    
    Parameters:
    
    atoms: Atoms
        | Atoms
    constraints: constraint object
        | ASEのconstraint.複数の制約をかける場合はリストで与える
    logfile (str, optional):
        | 作成するlogファイルのパス."-"の場合標準出力
    calc_func: function object
        | calculatorを返す関数
        
    Note:
        | ASEのOptimzizerとは異なり,引数のatoms自体は変化しない
        | 最適化後の構造はは,RepetitiveOpt.atoms に保存されている.
        
        >>> opt = RepetitiveOpt(atoms)
        >>> opt.run()
        >>> opt.atoms # <--これで最適化後のAtomsを取り出せる   
    """
    def __init__(self,atoms,constraints=[],logfile="-",calc_func=pfp_calculator):
        """異なる条件で繰り返し構造最適化を行なう

        Parameters:
        
        atoms: Atoms
            Atoms
        constraints: constraint object
            ASEのconstraint.複数の制約をかける場合はリストで与える
        logfile (str, optional):
            作成するlogファイルのパス."-"の場合標準出力
        calc_func: function object
            calculatorを返す関数
            
        Note:
            ASEのOptimzizerとは異なり,引数のatoms自体は変化しない
            最適化後の構造はは,RepetitiveOpt.atoms に保存されている.
        """
        self.atoms = atoms.copy()
        self.constraint = constraints
        self.calc_func = calc_func
        self.atoms.calc = self.calc_func()
        self.atoms.set_constraint(self.constraint)
        self.logfile = logfile
        self.attach_list = []
        self.update_list = []
        
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
        self.check_param()
            
    def check_param(self):
        attr = ["maxstep_list","fmax_list","steps_list"]
        attr_len = [len(getattr(self, i)) for i in attr]
        if not all([i==attr_len[0] for i in attr_len]):
            raise Exception(f"{', '.join(neb)}は同じ要素数である必要があります")
                
    def attach(self,*args,**kwargs):
        """ASEのOptimizerのattachと同じ"""
        self.attach_list.append((args,kwargs))
        
    def update(self,functions):
        """条件を変える前にここに設定した関数が走る"""
        self.update_list.append(functions)
    
    @property
    def default_param(self):
        """runのデフォルトの引数
        
        optimizer: Optimizer
            使用するOptimizer
        maxstep_list: list of float
            | maxstepのリスト
            | Noneの要素を指定した場合はfmaxに応じて自動でmaxstep値を変化させる.
            | どのように変化させるかはautomate_maxstepで指定できる.
        fmax_list: list of float
            famxのリスト
        steps_list: list of int
            stepsのリスト
        opt_kwargs: dict
            optimizerのインスタンス時に指定したい引数があれば設定する.
        automate_maxstep: dict
            grrmpy.optimizer.optimizer.automate_maxstepを参照
        """
        return {
            "optimizer": OPTIMIZER,
            "maxstep_list": [0.02,None],
            "fmax_list": [0.001,0.001],
            "steps_list": [200,10000],
            "opt_kwargs":{},
            "automate_maxstep":{}
        }
        
    def run(self,**param):
        """最適化計算を行なう.引数はdefault_paramを参照"""
        self.set_param(**param)
        for maxstep,fmax,steps in zip(self.maxstep_list,self.fmax_list,self.steps_list):
            ms = 0.2 if maxstep is None else maxstep
            if self.optimizer == FIRELBFGS:
                opt_param = {"maxstep_fire":ms,"maxstep_lbfgs":ms}
            else:
                opt_param = {"maxstep":ms}
            self.opt = self.optimizer(self.atoms,
                                      logfile=self.logfile,
                                      **opt_param,
                                      **self.opt_kwargs)
            for args,kwargs in self.attach_list:
                self.opt.attach(*args,**kwargs)
            if maxstep is None:
                self.opt.attach(lambda:automate_maxstep(self.opt,))
            converged = self.opt.run(fmax=fmax,steps=steps)
            for functions in self.update_list:
                functions()
        return converged
        