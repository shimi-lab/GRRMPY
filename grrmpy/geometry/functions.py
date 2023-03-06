# USER
from grrmpy.conv.atoms2smiles import atomslist2smileses,atoms2smiles
from grrmpy.functions import get_bond


def grouping(atoms_list, target0=None, target1=[], target2=[],mult=1,rtype="2D", **kwargs):
    """同一構造をまとめる.
    
    | atoms_list中で同じ構造ごとにindex番号をまとめたリストを返す.
    | 例えば,戻り値が[[0,1,3],[2,5,6],[4]]の場合,0,1,3番目の構造が同一構造である.
    | 
    | 結合の有無で同一構造かどうか判定を行なう.厳密にはSMILESが同一かを判定している.
    | ただし,全て単結合としてSMILESを表現するため,2重結合と単結合の違い(角度の違い)は表現できない.
    |
    | ゼオライトなどを扱う場合,活性点以外の結合状態の違いはあまり意味をなさない.
    | 引数のtargetを指定することで注目したい結合部分のみで同一構造判定を行なうことができる.

    Parameters:
    
    atoms_list: list of Atoms
        Atomsのリスト
    target0: list of int
        | target0,target1とtarget2との結合状態を見る.
        | Noneの場合,全ての結合状態を見る.
    target1: list of int
        | target0,target1の結合状態を見る.
    target2: list of int
        | target0の原子との結合状態を見る.
    rtype: str
        | 戻り値のタイプ. "1D" または "2D"
        | "1D"の場合,[0,0,1,2,3,3,4...]のような1次元リストを返す.この場合0,1番目と4,5番目が同じ構造という意味.
        | "2D"の場合,[[0,1],[2],[3],[4,5],[6]...]のようなリストを返す.これも0,1番目と4,5番目が同じ構造という意味.
    mult: float
        結合の判定を行なうための係数.数字の大きい程,遠く離れていても結合していると判断される.
    kwargs:
        | 原子毎に結合半径を指定できる.
        | ex) Ag=1.4 のように指定する.
    
    Hint:
        | targetの指定の仕方として, 例えばゼオライトにAg4クラスターを担持した触媒の場合,次のように指定すると良い.
        
        .. csv-table::
            :header: 計算モード, target0, target1, target2

            厳密な構造判定, None, [], []
            かなりラフな構造判定, 反応物, [], Ag
            ラフな構造判定(Agのコンフォメーションも考慮), 反応物+Ag, [], []
            ゼオライトケージへの吸着を考慮, 反応物, Ag, その他の原子
            ゼオライトケージへの吸着も考慮+Agのコンフォメーションも考慮, 反応物+Ag, [], その他の原子
    """
    smiles_list = atomslist2smileses(atoms_list,target0,target1,target2,mult,**kwargs)
    if rtype == "2D":
        group = [[i for i, _x in enumerate(smiles_list) if _x == smiles]
                for smiles 
                in sorted(set(smiles_list), key=smiles_list.index)]
    elif rtype == "1D":
        key = sorted(set(smiles_list), key=smiles_list.index)
        val = [i for i in range(len(key))]
        unique_dict = dict(zip(key,val))
        group = [unique_dict[smiles] for smiles in smiles_list]
    else:
        raise Exception("rtypeは'1D'または'2D'です")
    return group

def issame(atoms1,atoms2,target0=None,target1=[],target2=[],mult=1,**kwargs):
    """2つの構造が同一構造かどうかを判定する
    
    Atomsの__eq__とは異なり,ジオメトリが同じかどうかを判定する.
    """
    smiles1 = atoms2smiles(atoms1,target0,target1,target2,mult,**kwargs)
    smiles2 = atoms2smiles(atoms2,target0,target1,target2,mult,**kwargs)
    if smiles1 == smiles2:
        return True
    return False

def get_bond_diff(atoms1,atoms2,mult=1,**kwargs):
    """atoms1とatoms2の結合の差分
    
    Parameters:
    
    atoms1: Atoms
        Atomsオブジェクト
    atoms2: Atoms:
        Atomsオブジェクト
    
    retruns: Tuple(List[List[int,int]],List[List[int,int]])
        | 2要素のタプルを返す.
        | 1要素目は生成した結合(atoms1-->atoms2)のインデックス番号のリスト
        | 2要素目は開裂した結合(atoms1-->atoms2)のインデックス番号のリスト
        
    Note:
        | atoms1とatoms2が同じ組成であるかの検証は行わないので注意.
        | 当たり前だか異なる組成の構造を比較しても意味がない.
        
    """
    bond1 = get_bond(atoms1,mult=mult,**kwargs)
    bond2 = get_bond(atoms2,mult=mult,**kwargs)
    cleavage = []
    formation = []
    for i, (a,b) in enumerate(zip(bond1,bond2)):
        if not (set_a:=set(a))==(set_b:=set(b)):
            if len((dif:=set_a-set_b))>0:
                cleavage += [[i,int(d)] for d in dif]
            elif len((dif:=set_b-set_a))>0:
                formation += [[i,int(d)] for d in dif]
    return formation, cleavage