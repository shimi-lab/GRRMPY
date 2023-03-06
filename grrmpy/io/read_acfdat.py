import numpy as np
import pandas as pd
from pathlib import Path
from ase import Atoms
from ase.units import Bohr

# User
from grrmpy.data import zval
"""
基本はase.io.dader.pyはBaderChargeの計算方法に誤りがある(Daderのバージョンの問題が原因??)ので,
修正したコードを書く.
"""

def read_acf(file="ACF.dat")->pd.DataFrame:
    """
    
    ACF.datファイルの中身をDataFrameにして返す.

    Parameters:
    
    fileobj: Union[str, Path]
        ACF.dataファイルパス

    Returns:
        np.DataFrame: ACF.datのデータフレーム
    """
    with open(file,"r") as f:
        firstline = f.readline().strip().split()
    if not firstline == ['#', 'X', 'Y', 'Z', 'CHARGE', 'MIN', 'DIST', 'ATOMIC', 'VOL']:
        raise IOError("ACF.datの形式が異なるため解析できませんでした")
    df = pd.read_table(file,skiprows=2,skipfooter=4, engine='python',delim_whitespace=True, header=None)
    df.columns = ["#","X","Y","Z","CHARGE","MIN DIST","ATOMIC VOL"]
    return df

def get_dader(atoms:Atoms,file="ACF.dat",displacement=1e-4)->np.ndarray:
    """
    
    | ACF.datファイルからBaderChargeを計算する.
    | BaderCharge = ZVAL - Charge(ACF.dat内に記載の)
    | で計算される.
    

    Parameters:
    
    atoms: Atoms
        Atomsオブジェクト
    file: Union[str, Path, DataFrame]
        ACF.dataファイルパス, DataFrame.
    displacement: float
        原子の順番がACF.datとAtomsとで同じか検証するため.

    Returns:
        np.ndarray: BaderChargeのリスト
    """
    if isinstance(file,pd.DataFrame):
        df = file
    else:
        df = read_acf(file)
    
    #  原子の順番がACF.datとAtomsとで同じか検証
    if displacement is not None:  # check if the atom positions match
        for atom,x,y,z in zip(atoms,df["X"].tolist(),df["Y"].tolist(),df["Z"].tolist()):
            xyz = np.array([x,y,z])
            # ACF.dat units could be Bohr or Angstrom
            norm1 = np.linalg.norm(atom.position - xyz)
            norm2 = np.linalg.norm(atom.position - xyz * Bohr)
            assert norm1 < displacement or norm2 < displacement
    
    return np.array([zval[atom.number]-charge for atom, charge in zip(atoms,df["CHARGE"])])
        
    


# def attach_charges(atoms, fileobj='ACF.dat', displacement=1e-4):
#     if isinstance(fileobj, str):
#         with open(fileobj) as fd:
#             lines = fd.readlines()
#     else:
#         lines = fileobj

#     sep = '---------------'
#     i = 0  # Counter for the lines
#     k = 0  # Counter of sep
#     assume6columns = False
#     for line in lines:
#         if line[0] == '\n':  # check if there is an empty line in the
#             i -= 1           # head of ACF.dat file
#         if i == 0:
#             headings = line
#             if 'BADER' in headings.split():
#                 j = headings.split().index('BADER')
#             elif 'CHARGE' in headings.split():
#                 j = headings.split().index('CHARGE')
#             else:
#                 print('Can\'t find keyword "BADER" or "CHARGE".'
#                       ' Assuming the ACF.dat file has 6 columns.')
#                 j = 4
#                 assume6columns = True
#         if sep in line:  # Stop at last separator line
#             if k == 1:
#                 break
#             k += 1
#         if not i > 1:
#             pass
#         else:
#             words = line.split()
#             if assume6columns is True:
#                 if len(words) != 6:
#                     raise IOError('Number of columns in ACF file incorrect!\n'
#                                   'Check that Bader program version >= 0.25')

#             atom = atoms[int(words[0]) - 1]
#             atom.charge = atomic_numbers[atom.symbol] - float(words[j])
#             if displacement is not None:  # check if the atom positions match
#                 xyz = np.array([float(w) for w in words[1:4]])
#                 # ACF.dat units could be Bohr or Angstrom
#                 norm1 = np.linalg.norm(atom.position - xyz)
#                 norm2 = np.linalg.norm(atom.position - xyz * Bohr)
#                 assert norm1 < displacement or norm2 < displacement
#         i += 1
