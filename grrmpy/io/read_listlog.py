import numpy as np
from ase import Atoms
from ase.io import read

# USER
from grrmpy.io.read_com import frozen2atoms

def log2atoms(logfile,com=None,poscar=None,constraints=[]):
    """logファイルからAtomsオブジェクトに変換する

    Parameter:
    
    logfile: str or Path
        _list.logファイルパス.
    com:  str or Path
        | FrozenAtomsを追加する場合にcomファイルを指定.
        | その他comファイルの情報を確認したい時はgrrmpy.structure.comfile.COMを参照
    poscar: str or Path
        POSCARパスを設定した場合.セル情報を読み取りAtomsオブジェクトに適用する.
    constraints: ase.constraints
        | 与えたconstraintsは全てのAtomsに適用される.
        | 複数のconstraintsを与えたい場合はリストで与える
        
    Returns:
        Atoms: Atomsオブジェクト
    """
    logtext = _open_file(logfile)
    hash_idx, energy_idx = _read_mark_of_list(logtext)
    natoms = _get_natoms(logtext, energy_idx)
    positions_list = (_int2float(logtext[i+1:i+natoms+1]) for i in hash_idx) # 座標を取得
    elements = _get_chemical_symbols(logtext, energy_idx)
    if com:
        frozen = frozen2atoms(com)
        atoms_list = [Atoms(elements,positions)+frozen for positions in positions_list]
    else:
        atoms_list = [Atoms(elements,positions) for positions in positions_list]
    
    if poscar:
        pos_atoms = read(poscar,format="vasp")
        cell = pos_atoms.get_cell()
        for atoms in atoms_list:
            atoms.set_cell(cell)
            atoms.set_pbc(True)
    
    for atoms in atoms_list:
        atoms.set_constraint(constraints)
    return atoms_list

def read_positions(logfile):
    """logファイルからpositionsを読み取る
    
    Parameters:
    
    logfile: string
        *_list.logファイルパス

    Return:
        list of np.array: 座標のリスト
    """
    logtext = _open_file(logfile)
    hash_idx, energy_idx = _read_mark_of_list(logtext)
    natoms = _get_natoms(logtext, energy_idx)
    positions_list = (_int2float(logtext[i+1:i+natoms+1]) for i in hash_idx) # 座標を取得
    return positions_list

def read_energies(logfile):
    """logファイルからenergyのリストを返す
    
    Parameters:
    
    logfile: string
        *_list.logファイルパス
        
    Return:
        list of float: エネルギーのリスト
    
    """
    logtext = _open_file(logfile)
    _, energy_idx = _read_mark_of_list(logtext)
    return _logtext2energies(logtext, energy_idx)

def read_connections(logfile):
    """logファイルからCONNECTIONSを取得する
    
    Parameters:
    
    logfile: string
        TS or PT_list.logファイルパス
        
    Return:
        list of tuple: CONNECTIONのリスト
    """
    logtext = _open_file(logfile)
    return _read_connections(logtext)

def _read_connections(logtext):
    connections = [
        [int(text.split()[2]) if text.split()[2].isdecimal() else text.split()[2],
         int(text.split()[4]) if text.split()[4].isdecimal() else text.split()[4]]
         for text in logtext
         if "CONNECTION" in text]
        
    return connections
    

def _read_mark_of_list(logtext):
    """logファイルから#とEnergy    =のある行数を返す"""
    hash_idx = [i for i,text in enumerate(logtext) if "#" in text]
    energy_idx = [i for i,text in enumerate(logtext) if "Energy    =" in text]
    if len(hash_idx)*2 == len(energy_idx):
        """ReEnergyの場合,`Energy    =`が2つ出るため"""
        energy_idx = [i for i in energy_idx if i%2 == 0] #偶数個目のEnergy    =の値だけ読み取る
    return hash_idx, energy_idx

def _logtext2energies(logtext, energy_idx):
    energies = np.array([float(logtext[i].split()[2]) for i in energy_idx])
    return energies

def _get_chemical_symbols(logtext, energy_idx):
    return [text.split()[0] for text in logtext[3:energy_idx[0]]] # 元素のリストを取得

def _get_natoms(logtext, energy_idx):
    return len(_get_chemical_symbols(logtext, energy_idx))

def _int2float(positions):
    return [[float(xyz) for xyz in atom.split()[1:]] for atom in positions]

def _open_file(logfile):
    with open(logfile,"r") as f:
        text = f.readlines()
    return text