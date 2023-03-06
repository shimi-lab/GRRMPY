from pathlib import Path
import re
from natsort import natsorted
from ase.io import read,iread

def read_traj(folder,blank=True):
    """
    
    | フォルダ内のtrajファイルを検索し,Atomsのリストする.
    | trajファイル名には必ず数字が含まれている必要がある.
    | 例えば,0.traj, 3.traj, 2.trajファイルが存在する場合, 0.traj, 2.traj, 3.trajの順でファイルの読み込みを行ない,
    | blank=Trueの場合,[Atoms,None,Atoms,Atoms]のように数字の存在しないファイルはNoneの要素になる.
    | blank=Flaseの場合,None要素は作成されない.
    
    | ファイル名には数字を含む必要があるが,例えばsample_1_1.trajの場合,始めの数字である"1"のみに注目しソートされる.
    | 従って,blank=Trueの際にsample_1_1.traj, sample_1_2.trajの2つのファイルが存在するとき"1"のファイルが2つ存在すると判断されるためエラーとなる.
    | blank=Falseの場合このような問題は起きず,sample_1_1.traj, sample_1_2.trajの順で読み込まれる.

    Parameters:
    
    folder: str or Path
        フォルダ名

    Returns:
        tuple: Atomsのリスト,元々のファイルパス
    """
    return _read_traj_in_folder(folder,read,blank)

def iread_traj(folder,blank=True):
    """read_trajのimages版"""
    return _read_traj_in_folder(folder,lambda file: list(iread(file)),blank)

def _read_traj_in_folder(folder,read_func,blank=True):
    p = Path(folder)
    if not p.exists():
        raise IOError(f"{p.resolve()}フォルダが存在しません")

    traj_files = [p for p in p.glob("*traj")]
    if blank:
        file_num = [int(re.search(r'\d+', p.name).group()) for p in traj_files]
        max_i = max(file_num)
        file_dict = dict(zip(file_num,traj_files)) #{数字:ファイルパス}
        traj_files = [file_dict.get(i) for i in range(max_i+1)]
        atoms_list = [read_func(p) if not p is None else None for p in traj_files]
    else:
        traj_files = natsorted(traj_files,key=lambda p:p.name)
        atoms_list = [read_func(p) for p in traj_files]

    return atoms_list,traj_files