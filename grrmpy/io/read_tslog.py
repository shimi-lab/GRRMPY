from ase import Atoms

def read_irc(filename,images=None,step=None,min_n=8,max_n=64):
    """XXX_TSi.logファイルを読み込みAtoms情報のdictを要素とするリストを返す
    
    Parameters:

    filename: str
        XXX_TSi.logファイルパス
    images(-n): int
        | 始終構造を除くイメージの数.指定しない場合全ての構造を保存する
    steps(-s): int
        | (おおよそ)step個飛ばしに構造を抽出する, 25辺りが最適?
        | 偶数個抽出されるように調整してある.
        | min_n, max_nで最小,最大のイメージ数を決められる
    min_n: int
        最小イメージ数
    max_n: int
        最大イメージ数

    Returns:
        list of dict: Atoms情報のdictを要素とするリスト
    """
    with open(filename) as f:
        l = f.readlines()
    step_idx = [i for i,text in enumerate(l) if "# STEP" in text]
    separate_idx = [i for i,text in enumerate(l) if "Sphere optimization converged" in text][1]
    eq_idx = [i for i,text in enumerate(l) if "Optimized structure" in text]
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
    irc = [coordination2atoms(l[i+1:i+atoms_n+1]) for i in idx]
    return irc

def coordination2atoms(coordination:list):
    """logファイルの座標テキストからAtomsオブジェクトに変換する"""
    def int2float(text):
        text = text.split()
        return [float(text[1]),float(text[2]),float(text[3])]
    chemical_symbols = [text.split()[0] for text in coordination]
    cood = [int2float(text) for text in coordination]
    return Atoms(chemical_symbols,cood)

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