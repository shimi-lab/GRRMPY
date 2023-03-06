# coding: utf-8
from ase import Atoms
from natsort import natsorted
from pathlib import Path
import argparse

# user
import grrmpy.io.compressed_pickle as cpickle
from grrmpy.command.arg_formatter import CustomHelpFormatter

description = """
##########################################################
XXX_PTi.logのLUP構造を読みとり,pkl.bz2ファイルで出力する
##########################################################
    pkl.bz2ファイルはpickleのバイナリデータをbz2に圧縮したものである.
    """

parser = argparse.ArgumentParser(description=description, formatter_class=CustomHelpFormatter)

parser.add_argument('--index','-i', default=None, type=int,
                    help="指定したPT番号のLUPを作成する")
parser.add_argument('--outfile','-o', default=None, type=str,
                    help="指定した名前でファイルを作成する\n指定しない場合,LUP.pkl.bz2の名前で作成される")
parser.add_argument('--name', default=None, type=str,
                    help="ファイル名. XXX.comのXXXの部分.\n"+
                    "指定しない場合, フォルダにあるcomファイルから自動的に読みとる\n"+
                    "フォルダ内に複数のcomファイルが存在する時にはnameを設定する必要がある")
args = parser.parse_args()

def read_lup():
    """XXX_PTi.logのLUP構造を読みとり,pkl.bz2ファイルで出力する
    
    pkl.bz2ファイルはpickleのバイナリデータをbz2に圧縮したファイルである.

    Parameters:
    
    index(-i): int
        | 指定したPT番号のLUPを作成する
    outfile(-o):
        | 指定した名前でファイルを作成する
        | 指定しない場合, LUP.pkl.bzの名前で作成される
    name: str
        | ファイル名. XXX.comのXXXの部分.
        | 指定しない場合, フォルダにあるcomファイルから自動的に読みとる
        | フォルダ内に複数のcomファイルが存在する時にはnameを設定する必要がある
        
    Note:
        imagesとstepは,どちらかしか指定できない.
    """
    index = args.index
    outfile = args.outfile
    name = args.name
    return _read_lups(index,outfile,name)
    
def _read_lups(index=None,outfile=None,name=None):
    p = Path(".")
    if name is None:
        name = list(p.glob('*.com'))[0].stem
    if outfile is None:
        num = "" if index is None else index
        outfile = f"{name}_LUP{num}.pkl.bz2"
            
    logfile = natsorted(list(p.glob(f'{name}_PT*.log')))
    logfile = logfile[:-1] # 最後は*PT_list.logなので省く
    if type(index) == int:
        data = _read_lup(logfile[index])
    else:
        data = [_read_lup(log) for log in logfile]
    
    cpickle.dump(data, outfile)
        
def _read_lup(filename):
    """XXX_PTi.logファイルを読み込みAtoms情報のdictを要素とするリストを返す
    
    Parameters:
    
    filename: str
        XXX_PTi.logファイル
        
    Returns:
        list of dict: Atoms情報のdictを要素とするリスト
    """
    with open(filename) as f:
        l = f.readlines()
    node_idx = [i for i,text in enumerate(l) if "# NODE" in text]

    for i,text in enumerate(l[node_idx[0]+1:]):
        if len(text.split()) != 4:
            atoms_n = i # 原子の数
            break

    irc = [coordination2atoms_dict(l[i+1:i+atoms_n+1]) for i in node_idx]
    return irc
        
def coordination2atoms_dict(coordination:list):
    def int2float(text):
        text = text.split()
        return [float(text[1]),float(text[2]),float(text[3])]
    chemical_symbols = [text.split()[0] for text in coordination]
    cood = [int2float(text) for text in coordination]
    return Atoms(chemical_symbols,cood).todict()