import pickle
import bz2

def loads(comp):
    return pickle.loads(bz2.decompress(comp))

def dumps(obj):
    return bz2.decompress(pickle.dumps(obj))

def load(fname):
    """dumpで圧縮したbz2ファイルをpythonオブジェクトとして展開する
    
    Parameters
    
    fname:
        入力ファイル名(bz2ファイル)
    
    """
    fin = bz2.BZ2File(fname, 'rb')
    try:
        pkl = fin.read()
    finally:
        fin.close()
    return pickle.loads(pkl)

def dump(obj, fname, level=9):
    """pythonオブジェクトをbz2ファイルに圧縮しつつ,シリアライズ化する
    
    Parameter:
    
    obj: obj
        Pythonオブジェクト
    fname:
        出力ファイル名
    
    """
    pkl = pickle.dumps(obj)
    fout = bz2.BZ2File(fname, 'wb', compresslevel=level)
    try:
        fout.write(pkl)
    finally:
        fout.close()