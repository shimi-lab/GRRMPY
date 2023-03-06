# coding: utf-8
from ase import Atoms
from natsort import natsorted
from pathlib import Path
import argparse
import numpy as np

# user
import grrmpy.io.compressed_pickle as cpickle
from grrmpy.command.arg_formatter import CustomHelpFormatter

description = """
##########################################################
XXX_TSi.logのIRC構造を読みとり,pkl.bz2ファイルで出力する
##########################################################
    pkl.bz2ファイルはpickleのバイナリデータをbz2に圧縮したものである.
    
    imagesとstepは,どちらかしか指定できないので注意"""

parser = argparse.ArgumentParser(description=description, formatter_class=CustomHelpFormatter)
parser.add_argument('--images','-n', default=None, type=int,
                    help="始終構造を除くイメージの数.指定しない場合全ての構造を保存する")
parser.add_argument('--step','-s', default=None, type=int,
                    help="(おおよそ)step個飛ばしに構造を抽出する, 25辺りが最適?\n"+
                    "偶数個抽出されるように調整してある.\n"+
                    "min_n, max_nで最小,最大のイメージ数を決められる")
parser.add_argument('--index','-i', default=None, type=int,
                    help="指定したTS番号のIRCを作成する")
parser.add_argument('--outfile','-o', default=None, type=str,
                    help="指定した名前でファイルを作成する\n指定しない場合,IRC.pkl.bz2の名前で作成される")
parser.add_argument('--name', default=None, type=str,
                    help="ファイル名. XXX.comのXXXの部分.\n"+
                    "指定しない場合, フォルダにあるcomファイルから自動的に読みとる\n"+
                    "フォルダ内に複数のcomファイルが存在する時にはnameを設定する必要がある")
parser.add_argument('--min', default=8, type=int,
                    help="最小イメージ数")
parser.add_argument('--max', default=64, type=int,
                    help="最大イメージ数")
args = parser.parse_args()

def read_irc():
    """XXX_TSi.logのIRC構造を読みとり,pkl.bz2ファイルで出力する
    
    pkl.bz2ファイルはpickleのバイナリデータをbz2に圧縮したファイルである.

    Parameters:
    
    images(-n): int
        | 始終構造を除くイメージの数.指定しない場合全ての構造を保存する
    steps(-s): int
        | (おおよそ)step個飛ばしに構造を抽出する, 25辺りが最適?
        | 偶数個抽出されるように調整してある.
        | min_n, max_nで最小,最大のイメージ数を決められる
    index(-i): int
        | 指定したTS番号のIRCを作成する
    outfile(-o):
        | 指定した名前でファイルを作成する
        | 指定しない場合, IRC.pkl.bzの名前で作成される
    name: str
        | ファイル名. XXX.comのXXXの部分.
        | 指定しない場合, フォルダにあるcomファイルから自動的に読みとる
        | フォルダ内に複数のcomファイルが存在する時にはnameを設定する必要がある
    min_n: int
        最小イメージ数
    max_n: int
        最大イメージ数
        
    Note:
        imagesとstepは,どちらかしか指定できない.
    """
    images = args.images
    step = args.step
    index = args.index
    outfile = args.outfile
    name = args.name
    min_n = args.min
    max_n = args.max
    return _read_ircs(images,step,index,outfile,name,min_n,max_n)
    
def _read_ircs(images=None,step=None,index=None,outfile=None,name=None,min_n=8,max_n=64):
    p = Path(".")
    if name is None:
        name = list(p.glob('*.com'))[0].stem
    if outfile is None:
        num = "" if index is None else index
        outfile = f"{name}_IRC{num}.pkl.bz2"
            
    logfile = natsorted(list(p.glob(f'{name}_TS*.log')))
    logfile = logfile[:-1] # 最後は*TS_list.logなので省く
    if type(index) == int:
        data = _read_irc(logfile[index],images,step,min_n,max_n)
    else:
        data = [_read_irc(log,images,step,min_n,max_n) for log in logfile]
    
    cpickle.dump(data, outfile)
        
def _read_irc(filename,images=None,step=None,min_n=8,max_n=64):
    """XXX_TSi.logファイルを読み込みAtoms情報のdictを要素とするリストを返す
    
    Parameters:
    
    filename: str
        XXX_TSi.logファイル
        
    Returns:
        list of dict: Atoms情報のdictを要素とするリスト
    """
    with open(filename) as f:
        l = f.readlines()
    step_idx = [i for i,text in enumerate(l) if "# STEP" in text]
    separate_idx = [i for i,text in enumerate(l) if "Sphere optimization converged" in text][1]
    eq_idx = [i for i,text in enumerate(l) if "Optimized structure" in text]
    if len(eq_idx) != 2: # ini,finのどちらかが存在しない時
        return None
    ts_idx = 1
    for i,text in enumerate(l):
        if "ENERGY    =" in text:
            atoms_n = i-2 # 原子の数
            break
    for i in range(len(step_idx)):
        if step_idx[i] < separate_idx and separate_idx < step_idx[i+1]:
            forward_idx = step_idx[:i+1]
            reverse_idx = step_idx[i+1:]
            break
    forward_idx.reverse() #逆順にする
    idx = [eq_idx[0]] + forward_idx + [ts_idx] + reverse_idx + [eq_idx[1]]
    if images is not None:
        idx = n_sampling(idx,images-2)
    if step is not None:
        idx = step_sampling(idx,step,min_n,max_n)
    irc = [coordination2atoms_dict(l[i+1:i+atoms_n+1]) for i in idx]
    return irc
        
def coordination2atoms_dict(coordination:list):
    def int2float(text):
        text = text.split()
        return [float(text[1]),float(text[2]),float(text[3])]
    chemical_symbols = [text.split()[0] for text in coordination]
    cood = [int2float(text) for text in coordination]
    return Atoms(chemical_symbols,cood).todict()

def n_sampling(lst, n):
    """lst:listを与えたとき,いい感じに均等になるようにn個サンプリングする"""
    division = len(lst) / n
    return [lst[round(division * i):round(division * (i + 1))][0] for i in range(n)]

def step_sampling(lst,step,min_n,max_n):
    """n個飛ばしでサンプリングする"""
    n = int(len(lst)/step)
    n = n if n%2==0 else n+1
    if n < min_n+2:
        n = min_n+2
    elif n > max_n+2:
        n = max_n+2
    return n_sampling(lst, n)