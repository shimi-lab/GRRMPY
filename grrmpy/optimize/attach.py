from ase.neb import  NEB
from ase.io import Trajectory
from ase.units import mol,kJ,Hartree
from math import sqrt
from functools import partial 

# User
from grrmpy.calculator import pfp_calculator
from grrmpy.functions import get_diff,get_fmax

def optimize_eq(neb,calc_func=pfp_calculator):
    """iniとfinを最適化しながらNEB計算を行なう

    Parameters:
    
    neb: neb object
        NEBオブジェクト
    calc_func: function object
        calculatorを返す関数
    """    
    if get_fmax(neb.images[0]) > get_fmax(neb.images[1]):
        ini = neb.images[1].copy()
        ini.calc = calc_func()
        neb.images[0] = ini
    if get_fmax(neb.images[-1]) > get_fmax(neb.images[-2]):
        fin = neb.images[-2].copy()
        fin.calc = calc_func()
        neb.images[-1] = fin

     
neb_maxstep_climb_false = {
    10:0.01,
    5:0.05,
    2:0.07,
    1:0.1,
    0.1:0.2,
    0.07:0.3,
    0.06:0.4,
    0:0.5,
}

neb_maxstep_climb_true = {
    1:0.05,
    0:0.01,
}

opt_maxstep = {
    10:0.1,
    5:0.2,
    2:0.3,
    0:0.35,
}

def automate_maxstep(opt,maxstep=None):
    """maxstepの値を変更しながら最適化を行なう

    Parameters:
    
    opt: optimizer object
        Optimizerオブジェクト
    maxstep (_type_, optional)
        | {fmax:maxstep}の辞書で与える.
        | Noneの場合はNEB(climb=True,False),opt等に合わせて自動で設定する
    """
    if maxstep is None:
        if type(opt.atoms) == NEB:
            maxstep = neb_maxstep_climb_true if opt.atoms.climb else neb_maxstep_climb_false
        else:
            maxstep = opt_maxstep
            
    maxstep = sorted(maxstep.items(),key=lambda x:x[0],reverse=True)
    now_fmax = get_fmax(opt.atoms)
    for fmax,step in maxstep:
        if now_fmax > fmax:
            opt.maxstep = step
            break


def write_traj(outfile,opt,dist=0.25,mic=None,method=0):
    """最適化中の構造をある条件に従って保存する.
    
    | IRC Pathのようなものを作成することができる
    | 
    | method=0 : TS構造からの距離がdistずつ離れるたびに記録
    | method=1 : 緩和中の移動距離がdistずつ移動するたびに記録

    Parameters:
    
    outfile: str
        保存先のtrajファイルパス.既に存在するファイルを設定するとエラー
    opt: Optimizer object
        Optimizer
    dist: float
        dist Å構造が変化する度に構造を保存する.
    mic:
        | 周期境界条件で最小移動規則を適用し距離を算出する場合,True
        | Noneの場合,計算する構造が周期境界を持っている場合は自動でTrueにする.
    method: Int
        メソッド
        
    Note:
        | 他のattachとは異なりlambda文でのattachを行なわないことに注意する.
        | Exampleを参照する.
        
    Example:
    
        >>> opt = LBFGS(atoms)
        >>> opt.attach(write_traj('forward_path.traj',opt))
    """
    return partial(_write_traj, outfile, HoldAtoms(opt,dist,mic), dist, method)

def _write_traj(outfile,holdatoms,dist,method):
    """
    maxstepを変えて再計算した時にも追記できるような発展的な使用にはこれを使った方がいいかも
    """
    if holdatoms.first_write:
        mode = "w"
        holdatoms.first_write = False
    else:
        mode = "a"
    atoms1 = holdatoms.atoms
    atoms2 = holdatoms.opt.atoms
    if method == 0:
        method0(atoms1,atoms2,outfile,holdatoms,dist,mode)
    elif method == 1:
        method1(atoms1,atoms2,outfile,holdatoms,dist,mode)
    
def method0(atoms1,atoms2,outfile,holdatoms,dist,mode):
    if get_diff(atoms1,atoms2,mic=holdatoms.mic) > holdatoms.dist:
        Trajectory(outfile, mode=mode, atoms=atoms2).write()
        holdatoms.dist += dist
        
def method1(atoms1,atoms2,outfile,holdatoms,dist,mode):
    if get_diff(atoms1,atoms2,mic=holdatoms.mic) > dist:
        Trajectory(outfile, mode=mode, atoms=atoms2).write()
        holdatoms.atoms = atoms2.copy()

class HoldAtoms():
    """最適化時に構造を保持するためのクラス
    """
    def __init__(self,opt,dist,mic):
        self.first_write = True
        self.opt = opt
        self.atoms = opt.atoms.copy()
        self.dist = dist
        if mic is None:
            self.mic= True if any(self.atoms.get_pbc()) else False
        else:
            self.mic = mic
            

def same_energy_and_fmax(opt,times,row,force_digit=4,energy_digit=6):
    """同じEnergyとfmaxを繰り返す場合に停止する.

    OptWithCondと上で使用する.

    Parameters:
    
    opt: optimizer
        Optimizer(OptWithCond)
    times: int
        times回連続で同じfmax,energy繰り返す事態がrow回あれば停止する.
    row: int
        times回連続で同じfmax,energy繰り返す事態がrow回あれば停止する.
    force_digit: int
        小数点以下force_digit桁が一致するか
    energy_digit: integer
        小数点以下energy_digit桁が一致するか

    Examples:
    
        >>> from grrmpy.optimize import OptWithCond
        >>> opt = OptWithCond(FIRE,atoms,maxstep=0.2)
        >>> opt.attach_stop(same_energy_and_fmax(opt,10,2))
    """
    class data_class():
        def __init__(self,opt,times,row,force_digit=4,energy_digit=6):
            self.fmax = 0
            self.energy = 0
            self.opt = opt
            self.times = times
            self.times_cont = 0
            self.row = row
            self.row_count = 0
            self.force_digit = force_digit
            self.energy_digit = energy_digit

        def compare(self):
            pre_f = self.fmax
            now_f = get_fmax(self.opt.atoms)
            pre_e = self.energy
            now_e = self.opt.atoms.get_potential_energy()
            fd = self.force_digit
            ed = self.energy_digit
            if round(pre_f,fd) == round(now_f,fd) and round(pre_e,ed) == round(now_e,ed):
                if self.times_cont == self.times:
                    self.times_cont = 0
                    if self.row_count == self.row:       
                        return True
                    self.row_count += 1
                self.times_cont+=1  
            else:
                self.fmax = now_f
                self.energy = now_e
                return False

    def _same_energy_and_fmax(opt,data_class):
        return data_class.compare()
    
    return lambda c=data_class(opt,times,row,force_digit,energy_digit):_same_energy_and_fmax(opt,data_class=c)

def steep_peak(opt,max_ea,nsteps):
    """
    
    | NEB計算中にあるイメージと隣のイメージとのエネルギー差がmax_ea(kJ/mol)以上になれば計算を停止する.
    | 計算初期はイメージが高く出るのでnsteps(回)イタレーションが終わった後から計測を始める.

    Parameters:
    
    opt: optimizer
        Optimizerオブジェクト   
    max_ea: float
        隣のイメージとのエネルギー差がmax_ea(kJ/mol)以上になれば計算を停止する
    nsteps (_type_):
        始めのnsteps回のあとから計測を始める

    Examples:
    
        >>> from grrmpy.optimize import OptWithCond
        >>> neb = NEB(images,climb=False,parallel=True)
        >>> opt = OptWithCond(FIRE,neb,maxstep=0.2)
        >>> opt.attach_stop(lambda:steep_peak(opt,300,30))
    """
    energy_list = [i.get_potential_energy()*mol/kJ for i in opt.atoms.images]
    if opt.nsteps > nsteps:
        for i, _ in enumerate(energy_list[:-1]):
            if abs(energy_list[i+1]-energy_list[i]) > max_ea:
                return True
    return False