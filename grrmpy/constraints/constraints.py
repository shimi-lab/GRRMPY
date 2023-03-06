from ase.constraints import FixConstraint
from ase.geometry import find_mic
import numpy as np

class ForProduct(FixConstraint):
    """
    
    | 指定した終構造に向かって最適化される.
    | NEBイメージを作成する際にinterpolateでうまく作れない時等に使用するとよいかも.

    """
    def __init__(self, ini_atoms, fin_atoms, f_ext=2, indices=None, dif=0.3, mic=None):
        """
        Parameters:
        
        ini_atoms: Atoms
            初期構造(最適化計算を行なう構造)
        fin_atoms:
            終構造,予め最適化を行なっておく必要がある.
        f_ext:flaot
            終構造に向かって行く際にかける力
        indices:
            | 反応に関わる原子のindex番号.
            | 指定しなくても計算可能，
            | 指定しないと収束しない場合があるので指定することを勧める.
        dif: float
            終構造との原子の位置の差がdif(Å)以下の場合,力を加えない
        mic: bool
            周期境界の場合は自動でTrueにする(mic=Noneの場合)
            
        Examples:
        
            >>> from grrmpy.optimize.attach import write_traj
            >>> f = ForProduct(ini,fin,indices=[i for i in range(20,28)])
            >>> ini.set_constraint(f)
            >>> opt =LBFGS(ini,maxstep=0.05)
            >>> opt.attach(write_traj("NEB_image.traj",opt,dist=0.25,method=1)) #NEBイメージ作成が目的の場合
            >>> opt.run()
        """
        if mic is None:
            self.mic = True if any(ini_atoms.get_pbc()) else False
        else:
            self.mic = mic
        self.cell = ini_atoms.get_cell()
        self.pbc = ini_atoms.pbc
        
        if indices is None:
            self.indices = [False for _ in range(len(ini_atoms))]
        else:
            self.indices = [False if i in indices else True for i in range(len(ini_atoms))]
            
        self.external_force = f_ext
        self.objective_pos = fin_atoms.get_positions()
        self.dif = dif
        
        self.first_dif = fin_atoms.get_positions()-ini_atoms.get_positions()
        if self.mic:  
            self.first_dif = find_mic(self.first_dif, self.cell, self.pbc)[0]
        
    def adjust_forces(self, atoms, forces):
        """self.difより変化量が小さい時はforceを加えない"""
        dist = self.objective_pos - atoms.positions
        if self.mic:
            dist = find_mic(dist, self.cell, self.pbc)[0]
        ## 元々fin構造に近かったのに離れてしまったものに大きなforceをかけるため
        coef0 = (np.abs(np.linalg.norm(dist,axis=1)/np.linalg.norm(self.first_dif,axis=1))).reshape(-1, 1)
        ## fin構造との差分の絶対値を取得 ##
        coef1 = np.linalg.norm(dist,axis=1).reshape(-1, 1)
        ## fin構造との差分の単位ベクトルを算出 ##
        norm_dif_f = dist/coef1
        ## forceの単位ベクトルを取得 ##
        norm_opt_f = forces/np.linalg.norm(forces,axis=1).reshape(-1, 1)
        ## 係数2(optの示すforceのベクトルとdifとの差分のベクトルのcos類似度を算出) ##
        coef2 = np.array([self.cos_similarity(i,j)
                                for i,j in zip(norm_dif_f,norm_opt_f)]).reshape(-1, 1)
        ## 外力を算出 ##
        mask = np.linalg.norm(dist,axis=1) < self.dif
        mask[self.indices] = True
        coef0[mask],coef1[mask],coef2[mask] = 0,0,0
        ex_force = self.external_force * norm_dif_f * (coef0*8*coef1**2 + coef1 + coef2)
        forces += ex_force
    
    def adjust_positions(self, atoms, new):
        pass
    
    def cos_similarity(a,b):
        """コサイン類似度(から1引き,符号変え,2割ったもの)
        0の時完全に一致,1の時180異なる.
        """
        cos_sim = np.inner(a,b)/(np.linalg.norm(a)* np.linalg.norm(b))
        cos_sim = -(cos_sim-1)/2
        return cos_sim