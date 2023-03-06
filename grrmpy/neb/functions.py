from ase.geometry import find_mic
from ase.units import mol, kJ, Hartree
from scipy import signal
import numpy as np

# User
from grrmpy.calculator import pfp_calculator
from grrmpy.functions import to_html_energy_diagram, get_energy_list
from grrmpy.functions import set_calc2images

def to_html_nebgraph(neb_obj,calc_func=pfp_calculator,
                     full_html=False,
                     unit="kJ/mol",
                     title="NEB Energy Diagram",
                     xaxis_title="Reaction Coordinate",
                     yaxis_title=None,
                     highlight=None,
                     annotation=[],
                     **kwargs):
    """NEBのエネルギーダイアグラムのhtmlテキストを作成する
    
    Parameters:
    
    unit: string
        'eV', 'kJ/mol', 'Hartree', 'kcal/mol'のいずれか
    yaxis_title: string
        Noneの場合,Energy({unit})
    kwargs:
        plotpyのto_htmlの引数を参照
    """
    return to_html_energy_diagram(neb_obj.images,
                                  calc_func=calc_func,
                                  full_html=full_html,
                                  unit=unit,
                                  title=title,
                                  xaxis_title=xaxis_title,
                                  yaxis_title=yaxis_title,
                                  highlight=highlight,
                                  annotation=annotation,
                                  **kwargs)
    
def insert_image(neb, mic=False, i=None, a=0.01, clac_func=pfp_calculator):
    """TS周りに新たにimageを挿入する.
    
    TS構造の両側に新たなImageを挿入する.
    
    Parameters:
    
    neb: NEB object
        NEBオブジェクト
    mic: bool
        周期境界条件で最小画像規則を使用する場合は,True
    i: int or None
        | i番目のimageの両側に新たなイメージを作成する.
        | Noneの場合,imax(TS)に挿入する.
    a: float or list of float
        | 挿入する構造のTS構造との変位
        | i番目の構造からa[Å]離れた構造を挿入する
        | 2要素のリストで与えた場合,i-1番目に1番目要素,i+1番目に2番目の要素を適用する
    clac_func: fnction object
        claculatorを返す関数. デフォルトはpfpのcalculator
    """
    if i is None:
        i = neb.imax
    if type(a) == int or type(a) == float:
        a = [a,a]
    elif type(a) == list or type(a) == tuple:
        a = a
    else:
        raise TypeError("aはfloatまたは要素数2のリストです")
    images = neb.images
    ts = neb.images[i].copy()
    def insert(idx,ins_idx,a):
        """
        idx: tsの隣の構造のidx番号
        ins_idx: 挿入する位置
        """
        pos1 = images[idx].get_positions()
        pos2 = ts.get_positions()
        d = pos2 - pos1
        if mic:
            d = find_mic(d, ts.get_cell(), ts.pbc)[0]
        n = np.linalg.norm(d)/a
        d /= (n - 1.0)
        new_pos = pos1 + (n-2) * d # n-2(ts付近)
        unconstrained_image = ts.copy()
        unconstrained_image.set_positions(new_pos,apply_constraint=True)
        unconstrained_image.calc = clac_func()
        images.insert(ins_idx,unconstrained_image)
    if i != len(images)-1:
        insert(i+1,i+1,a[1])
    if i != 0:
        insert(i-1,i,a[0])


def is_barrier_less(neb_images,barrier_less=0,calc_func=pfp_calculator):
    """NEBイメージがバリアレスかどうか
    
    Parameters:
    
    neb_images: list of Atoms
        NEBイメージ
    barrier_less: float
        barrier_less(kJ/mol)以下のTSであればバリアレスと見做す    
    """
    try:
        energy_list = np.array([i.get_potential_energy()*mol/kJ for i in neb_images])
    except:
        for image in neb_images:
            image.calc = calc_func()
        energy_list = np.array([i.get_potential_energy()*mol/kJ for i in neb_images])      
    max_idxes = signal.argrelmax(np.array(energy_list))[0]
    if len(max_idxes) == 0:
        # 単純に極大値がない時(純粋なバリアレスの場合)
        return True
    if max(energy_list)-energy_list[-1] < barrier_less:
        # 上昇系のバリアレスの場合
        energy_list = list(reversed(energy_list)) # 反転し,下降系のエネルギーダイアグラムにする.
    elif max(energy_list)-energy_list[0] < barrier_less:
        # 下降系のバリアレスの場合
        pass
    else:
        # 完全にバリアレスでない場合
        return False
    # これより下はバリアレスか怪しいものを判定する
    # 確実にバリアレス,バリレスでないとっ判断されたものは既にTrue,Falseを返している
    max_idxes = signal.argrelmax(np.array(energy_list))[0]
    min_idxes = signal.argrelmin(np.array(energy_list))[0]
    if len(max_idxes) == 1 and len(min_idxes) == 0:
        # 極大値が1つのみの時
        return energy_list[max_idxes[0]] - energy_list[0] < barrier_less
    elif len(max_idxes) == 0 and len(min_idxes) == 1:
        # 極小値が1つのみの時
        return energy_list[-1] - energy_list[max_idxes[0]] < barrier_less
    # これより下は極大値,極小値が複数ある,より複雑なものを判定する.
    indexes = sorted([0] + list(max_idxes) + list(min_idxes) + [len(energy_list)-1]) #極値のindex
    if max_idxes[0] > min_idxes[0]:
        indexes = indexes[:-1] if len(indexes) %2 != 0 else indexes
    else:
        indexes = indexes[1:] if len(indexes) %2 != 0 else indexes
    for i in indexes[:-1]:
        if energy_list[i+1]-energy_list[i] > barrier_less:
            return False
    return True

def separate_images(neb_images,tolerance=0,threshold=0):
    """
    作成中
    
    | NEBイメージに複数の山があった際,一番高い山のみを分離抽出する.
    | イメージは書き変わる(破壊的メソッド)なので注意する.

    Parameters:
    
    neb_images: list of Atoms
        NEBイメージ,実行後には書き変わっている.
    tolerance: float
        tolerance(kJ/mol)以下の谷間は分離しない.
    threshold: float
        一番高い山(TSの山)のみを抽出した後に,両端の裾野の部分をさらにTSの高さに対しthreshold%削る.
    """
    _separate_images(neb_images,tolerance,threshold)


def _separate_images(neb_images,tolerance=0,threshold=0,not_considered_inifin=True,calc_func=pfp_calculator):
    """
    作成中
    
    | NEBイメージに複数の山があった際,一番高い山のみを分離抽出する.
    | イメージは書き変わる(破壊的メソッド)なので注意する.

    Parameters:
    
    neb_images: list of Atoms
        NEBイメージ,実行後には書き変わっている.
    tolerance: float
        tolerance(kJ/mol)以下の谷間は分離しない.
    threshold: float
        一番高い山(TSの山)のみを抽出した後に,両端の裾野の部分をさらにTSの高さに対しthreshold%削る.
    not_considered_inifin: bool
        NEB計算の際と同じようにini,finを除いて最もエネルギーの高い点をTSとする.
    calc_func: function object
        Calculatorを返す関数
        
    Return:
        ini_idx, fin_idx
    """
    # エネルギーの取得
    try:
        energy_list = np.array([i.get_potential_energy()*mol/kJ for i in neb_images])
    except:
        for image in neb_images:
            image.calc = calc_func()
        energy_list = np.array([i.get_potential_energy()*mol/kJ for i in neb_images])
    # TSのindex番号の取得  
    if not_considered_inifin:
        imax = numpy.argmax(energy_list[1:-1])+1
    else:
        imax = numpy.argmax(energy_list)
    # TSを境にエネルギーリストを2つに分離する
    reverse_energy = list(reversed(energy_list[:imax+1])) # 下降系のダイアグラムにする.
    reverse_atoms = list(reversed(neb_images[:imax+1]))
    forward_energy = energy_list[imax:]
    forward_atoms = neb_images[imax:]
    
    max_idxes = signal.argrelmax(np.array(reverse_energy))[0]
    min_idxes = signal.argrelmax(np.array(reverse_energy))[0]
    idxes = sorted(list(max_idxes) + list(min_idxes) + [len(reverse_energy)-1])
    
def get_ea(images,reverse=None,calc_func=pfp_calculator,unit="kJ/mol"):
    """活性化エネルギーを算出する
    
    Parameters:
    
    images: list of Atoms
        | Atomsのリスト(NEBイメージ)
    reverse: bool
        | Trueの場合,逆反応のEaを返す. Noneの場合,どちらか高い方を返す.
    calc_func: function object
        | calculatorを返す関数
    unit: str
        | 戻り値の単位(kJ/mol, eV, Hartreeのいずれか)
        
    Return:
        float: Ea値
    """
    energy_list = np.array(get_energy_list(images,calc_func))
    imax = energy_list.argmax()
    ts_e = energy_list.max()
    if reverse is None:
        ini_e = energy_list.min()
    elif reverse: # reverse is Ture
        ini_e = energy_list[imax:].min()
    else: # reverse is False
        ini_e = energy_list[:imax].min()

    if unit == "eV":
        conv = 1
    elif unit == "kJ/mol":
        conv = mol/kJ
    elif unit == "Hartree":
        conv = 1/Hartree
    else:
        raise Exception("unitはeV,kJ/mol,Hartreeのいずれかです")
    return (ts_e-ini_e)*conv

def find_local_minimum_points_and_ea(images,minenergy=5,calc_func=pfp_calculator):
    """
    NEBイメージの中から極小値を見つけそのindex番号を返す
    
    Parameters:
    
    images: list of Atoms
        Atomsのリスト(NEBイメージ)
    minenergy: int
        | 極小点を見つけたとしてもその両側の山との高さの差がどちらもminenergy以下であればそれは極小値とは見做さない.
    
    """
    energy_list = get_energy_list(images,calc_func)
    energy_list = np.array([energy*mol/kJ for energy in energy_list])
    min_idxes = signal.argrelmin(energy_list)[0].tolist()
    max_idxes = signal.argrelmax(energy_list)[0].tolist()
    marge_idxes = sorted(min_idxes+max_idxes)
    if len(marge_idxes) == 0:
        return [],[]
    # 極大値,最小値,極大値,最小値......極大値になるようにする
    if not marge_idxes[0] in max_idxes:
        marge_idxes.insert(0,0)
    if not marge_idxes[-1] in max_idxes:
        marge_idxes.append(len(energy_list)-1)
    new_min_index = []
    ea_list = []
    for i in range(1,len(marge_idxes)-1,2):
        ea1 = energy_list[marge_idxes[i-1]]-energy_list[marge_idxes[i]]
        ea2 = energy_list[marge_idxes[i+1]]-energy_list[marge_idxes[i]]
        if ea1>minenergy or ea2>minenergy:
            new_min_index.append(marge_idxes[i])
            ea_list.append(max([ea1,ea2]))
    return new_min_index,ea_list
    
def get_appropriate_nimages(ini_atoms,fin_atoms,dist=0.35,nmax=24,nmin=8,mic=None):
    """最適なNEBのイメージ数を算出する.
    
    最も移動距離の大きい原子が,d(Å)ずつ移動する時のイメージ数を返す.
    maxとminで最大,最小のイメージ数の制限を与えることができる.

    Parameters:
    
    ini_atoms: Atoms
        | Atoms
    fin_atoms: Atoms
        | Atoms
    dist: float
        | d(Å)毎にイメージ切る
    nmax: int
        | 最大のimage数
    nmin:
        | 最小のimage数.
    mic: bool
        | 最小イメージ規則を適用するか.Noneの場合与えたAtomsが周期境界条件である時,自動でTureにする.
    """
    pos1 = ini_atoms.get_positions()
    pos2 = fin_atoms.get_positions()
    d = pos2 - pos1
    if mic is None:
        mic = True if any(ini_atoms.get_pbc()) else False
    if mic:
        d = find_mic(d, ini_atoms.get_cell(), ini_atoms.pbc)[0]
    nimages = int(np.linalg.norm(d,axis=1).max()/dist)
    if nimages > nmax:
        nimages = nmax
    elif nimages < nmin:
        nimages = nmin
    return nimages

def get_ts_image(images,barrierless_is_none=True,copy=True,calc_func=pfp_calculator):
    """
    
    NEBイメージからTSのイメージを抽出する.
    最もエネルギーの高いimageをTSとする.
    
    Parameters:
    
    images: list of Atoms
        | Atomsのリスト
    barrierless_is_none:bool
        | Trueの場合,バリアレスの時Noneを返す
        | Falseの場合,最もエネルギーの高いimageを返す
    copy: bool
        | コピーを作成して返す場合True.
    calc_func: functions object
        | calculatorを返す関数
    """
    energy_list = np.array(get_energy_list(images, calc_func))
    imax = energy_list.argmax()
    max_idxes = signal.argrelmax(energy_list)[0]
    min_idxes = signal.argrelmin(energy_list)[0]
    if barrierless_is_none and len(max_idxes)==0 and len(min_idxes==0):
        return None
    if copy:
        return images[imax].copy()
    else:
        return images[imax]
    
def get_imax(images,threshold=5,calc_func=pfp_calculator)->int:
    """TSのindex番号を返す
    
    | バリアレス(極大値・極小値: 0個)の場合Noneを返す.
    | 山が1つ(極大値: 1個,極小値: 0個)かつ,活性化障壁がthreshold[kJ/mol]以下であればバリアレスと考えNoneを返す.
    | (正反応,逆反応のとちらかがthreshold[kJ/mol]以下の場合)
    
    Parameters:
    
    images: list of Atoms
        | NEBイメージ
    threshold: float
        | 活性化障壁がthreshold(kJ/mol)以下の場合はバリアレスと判定する
    calc_func: function object
        | imagesにcalcuatorが設定されていない場合に必要
        
    Returns:
        int or None: TSのindex番号,バリアレスの場合None.
    """
    energy_list = np.array(get_energy_list(images, calc_func))
    max_idxes = signal.argrelmax(energy_list)[0]
    min_idxes = signal.argrelmin(energy_list)[0]
    if len(max_idxes)==0 and len(min_idxes)==0:
        return None
    imax = energy_list.argmax()
    if len(max_idxes)==1 and len(min_idxes)==0:
        ea1 = (energy_list[imax]-energy_list[0])*mol/kJ # 正反応Ea
        ea2 = (energy_list[imax]-energy_list[-1])*mol/kJ # 逆反応Ea
        if ea1 < threshold or ea2 < threshold:
            return None
    return imax

