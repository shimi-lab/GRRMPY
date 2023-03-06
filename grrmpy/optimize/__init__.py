from grrmpy.optimize.attach import optimize_eq,automate_maxstep,write_traj,same_energy_and_fmax,steep_peak
from grrmpy.optimize.optimizer import OptWithCond
import warnings

try:
    from matlantis_features.ase_ext.optimize import FIRELBFGS #今後matlantis_featuresのアップデートの際に場所が変更される恐れがあるため.
    __all__ = ["optimize_eq","automate_maxstep","write_traj","same_energy_and_fmax","steep_peak",
               "FIRELBFGS",
               "OptWithCond"]
    
except:
    warnings.warn('matlantis_featuresのFIRELBFGSのディレクトリの位置が変更されたためインポートできませんでした')
    __all__ = ["optimize_eq","automate_maxstep","write_traj","same_energy_and_fmax","steep_peak",
               "OptWithCond"]