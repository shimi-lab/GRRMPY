import math
from math import sqrt
import numpy as np
import plotly.graph_objects as go
import plotly.io as pyi
import networkx as nx
from ase.units import kJ,Hartree,mol,kcal
from ase import Atoms
from ase.neb import NEB
from ase.geometry import find_mic
from ase.geometry.analysis import Analysis
from ase.build.rotate import rotation_matrix_from_points
from ase.neighborlist import build_neighbor_list,natural_cutoffs
# User
from grrmpy.calculator import pfp_calculator
# 本来デフォルトのレンダラーを使用したいが，なぜかうまく作動しないため追加
import plotly.io as pio
pio.renderers.default = "sphinx_gallery"

def get_fmax(atoms):
    """fmaxを返す
    
    Parameters:
    
    atoms: Atoms or NEB
    
    """
    return sqrt((atoms.get_forces()**2).sum(axis=1).max())

def copy_images(images):
    """imagesをコピーする"""
    return [image.copy() for image in images]

def set_calc2images(images,calc_func=pfp_calculator):
    """imagesにcalculatorを設定する."""
    for image in images:
        image.calc = calc_func()
        
def set_const2images(images,constraints):
    """imagesにconstraintsをかける."""
    for image in images:
        image.set_constraint(constraints)
        
def get_energy_list(images,calc_func):
    try:
        energy_list = [image.get_potential_energy() for image in images]
    except:
        set_calc2images(images,calc_func)
        energy_list = [image.get_potential_energy() for image in images]
    return energy_list

def draw_graph(images,
               highlight=None,
               annotation=[], #[((x,y),"text",{kwargs})]
               unit="kJ/mol",
               title="Energy Diagram",
               xaxis_title="",
               yaxis_title=None,
               calc_func=pfp_calculator,):
    """エネルギーダイアグラムのグラフを作成する
    
    Parameters:
    
    images: list of Atoms
        | Atomsオブジェクトのリスト
    highlight: list of int
        指定したindex番号のプロットのみ赤色になる
    annotation:list
        | アノテーションを入れたい場合に設定
        | [((x座標,y座標),"テキスト",{その他}),]のタプルのリストで与える
        | y座標もなくても良い,[x座標,テキスト]で与える
        | {その他}はなくても良い.(2要素で与える)
        | その他には例えば
        | "showarrow":True 矢印で示す
        | "arrowsize":2 矢印の大きさ
        | "arrowhead":3 矢印の種類
        | "ax:0" 矢印の向き
        | "ay":-50 矢印の向きなど
    unit: string
        | 'eV', 'kJ/mol', 'Hartree', 'kcal/mol'のいずれか
    xaxis_title: str
        | x軸のタイトル  
    yaxis_title: string
        | Noneの場合,Energy({unit})
    
    """   
    def find_y_value(x,y_list):
        if float(x).is_integer():
            return y_list[int(x)]
        else:
            x_ini,x_fin = int(x),int(x)+1
            y_ini,y_fin = y_list[x_ini],y_list[x_fin]
            y = (x-x_ini)*(y_fin-y_ini)/(x_fin-x_ini)+y_ini
            return y
      
    x = [i for i in range(len(images))]
    try:
        y = [atoms.get_potential_energy() for atoms in images]
    except:
        for image in images:
            image.calc = calc_func()
        y = [atoms.get_potential_energy() for atoms in images]
           
    y = [i-y[0] for i in y] # iniのエネルギーを0スタートで表記
    if unit == "kJ/mol":
        y = [i*mol/kJ for i in y]
    elif unit == "Hartree":
        y = [i/Hartree for i in y]
    elif unit == "kcal/mol":
        y = [i*mol/kcal for i in y]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=y,
                   line=dict(width=2,color="black"))
    )
    if highlight:
        fig.add_trace(
            go.Scatter(x=highlight, y=[y[i] for i in highlight],
                       mode="markers",
                       marker=dict(size=12,color="red"))
        )
        
    for position,*args in annotation:
        if type(position)==tuple or type(position)==tuple:
            x_pos,y_pos = position
        else:
            x_pos = position
            y_pos = find_y_value(x_pos,y)
        if len(args) == 1:
            text = args[0]
            kwargs = {}
        else:
            text,kwargs = args
        fig.add_annotation(x=x_pos, y=y_pos, text=text, **kwargs)
        
    if yaxis_title is None:
        yaxis_title = f"Energy({unit})"
    fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=unit,
        showlegend=False,
    )
    return fig
    
def to_html_energy_diagram(images,calc_func=pfp_calculator,
                           full_html=False,
                           unit="kJ/mol",
                           title="Energy Diagram",
                           xaxis_title="",
                           yaxis_title=None,
                           highlight=None,
                           annotation=[],
                           **kwargs):
    """エネルギーダイアグラムのhtmlテキストを作成する
    
    Parameters:
    
    images: list of Atoms
        Atomsオブジェクトのリスト
    highlight: list of int
        指定したindex番号のプロットのみ赤色になる
    annotation:list
        | アノテーションを入れたい場合に設定
        | [((x座標,y座標),"テキスト",{その他}),]のタプルのリストで与える
        | y座標もなくても良い,[x座標,テキスト]で与える
        | {その他}はなくても良い.(2要素で与える)
        | その他には例えば
        | "showarrow":True 矢印で示す
        | "arrowsize":2 矢印の大きさ
        | "arrowhead":3 矢印の種類
        | "ax:0" 矢印の向き
        | "ay":-50 矢印の向きなど
    full_html: bool
        | <html>タグたか始まる完全なhtmlテキストを出力する場合,True
        | Falseの場合<div>タグのみ
    unit: string
        'eV', 'kJ/mol', 'Hartree', 'kcal/mol'のいずれか
    xaxis_title: str
        x軸のタイトル  
    yaxis_title: string
        Noneの場合,Energy({unit})
    
    """
    fig = draw_graph(
        images,
        highlight = highlight,
        annotation = annotation,
        calc_func=calc_func,
        unit=unit,
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title)
    return pyi.to_html(fig,full_html=full_html,**kwargs)


def circle_fitting(x, y):
    """Circle Fitting with least squared
        input: point x-y positions  
        output  cxe x center position
                cye y center position
                re  radius of circle 
    """

    sumx = sum(x)
    sumy = sum(y)
    sumx2 = sum([ix ** 2 for ix in x])
    sumy2 = sum([iy ** 2 for iy in y])
    sumxy = sum([ix * iy for (ix, iy) in zip(x, y)])

    F = np.array([[sumx2, sumxy, sumx],
                  [sumxy, sumy2, sumy],
                  [sumx, sumy, len(x)]])

    G = np.array([[-sum([ix ** 3 + ix * iy ** 2 for (ix, iy) in zip(x, y)])],
                  [-sum([ix ** 2 * iy + iy ** 3 for (ix, iy) in zip(x, y)])],
                  [-sum([ix ** 2 + iy ** 2 for (ix, iy) in zip(x, y)])]])

    try:
        T = np.linalg.inv(F).dot(G)
    except np.linalg.LinAlgError:
        return 0, 0, float("inf")

    cxe = float(T[0] / -2)
    cye = float(T[1] / -2)

    try:
        re = math.sqrt(cxe ** 2 + cye ** 2 - T[2])
    except np.linalg.LinAlgError:
        return cxe, cye, float("inf")
    return cxe, cye, re


def calc_curvature_circle_fitting(x, y, npo=1):
    """各点での曲率を求める
    Calc curvature
    x,y: x-y position list
    npo: the number of points using Calculation curvature
    ex) npo=1: using 3 point
        npo=2: using 5 point
        npo=3: using 7 point
    """
    cv = []
    n_data = len(x)

    for i in range(n_data):
        lind = i - npo
        hind = i + npo + 1

        if lind < 0:
            lind = 0
        if hind >= n_data:
            hind = n_data

        xs = x[lind:hind]
        ys = y[lind:hind]
        (cxe, cye, re) = circle_fitting(xs, ys)

        if len(xs) >= 3:
            # sign evaluation
            c_index = int((len(xs) - 1) / 2.0)
            sign = (xs[0] - xs[c_index]) * (ys[-1] - ys[c_index]) - (
                    ys[0] - ys[c_index]) * (xs[-1] - xs[c_index])

            # check straight line
            a = np.array([xs[0] - xs[c_index], ys[0] - ys[c_index]])
            b = np.array([xs[-1] - xs[c_index], ys[-1] - ys[c_index]])
            theta = math.degrees(math.acos(
                np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))))

            if theta == 180.0:
                cv.append(0.0)  # straight line
            elif sign > 0:
                cv.append(1.0 / -re)
            else:
                cv.append(1.0 / re)
        else:
            cv.append(0.0)

    return cv

def n_sampling(lst, n, end=False):
    """listを与えたとき,いい感じに均等になるようにn個サンプリングする
    
    Parameters:
        lst: list
            リスト
        n: int
            抽出する数
        end: bool
            始めと終わりの構造を必ず含める場合True
    
    Return
        list: サンプリング後のリスト
    """
    division = len(lst) / n
    a = [lst[round(division * i):round(division * (i + 1))][0] for i in range(n)]
    if end:
        original_idx = [i for i in range(len(lst))]
        idx = [original_idx[round(division * i):round(division * (i + 1))][0] for i in range(n)]
        if original_idx[-1] != idx[-1]:
            a[-1] = lst[-1]
    return a

def get_diff(atoms1,atoms2,mic=False):
    """atoms1とatoms2の距離の差分を算出する

    Parameters:
    
    atoms1: Atoms
        Atomsオブジェクト
    atoms2: Atoms
        Atomsオブジェクト
    mic: bool
        周期境界条件で最小移動規則を適用する場合.True

    Returns:
        float: 距離
    """
    pos1 = atoms1.get_positions()
    pos2 = atoms2.get_positions()
    d = pos2 - pos1
    if mic:
        d = find_mic(d, atoms2.get_cell(), atoms2.pbc)[0]
    d = np.linalg.norm(d)
    return d


def minimize_rotation_and_translation_for_specified_indices_only(target, atoms, indices=None):
    """atomsの原子の位置をtargetの位置と近くなるように(最小二乗法)配置する.
    
    ase.build.rotate.minimize_rotation_and_translationでindexを指定できるようにした関数
    indices=Noneの時はASEのminimize_rotation_and_translationと同じで全ての原子を動かす.
    """
    if indices is None:
        indices = [i for i in range(len(target))]
    p = atoms.get_positions()[indices]
    p0 = target.get_positions()[indices]

    # centeroids to origin
    c = np.mean(p, axis=0)
    p -= c
    c0 = np.mean(p0, axis=0)
    p0 -= c0

    # Compute rotation matrix
    R = rotation_matrix_from_points(p.T, p0.T)

    atoms.positions[indices] = np.dot(p, R.T) + c0
    
def get_bond(atoms, unique=True, mult=1.0, **kwargs):
    """結合している原子のindex番号のリストを返す
    
    Parameters:
    
    atoms: Atoms
        対象のAtoms
    unique: bool
        Trueの場合,A-B,B-Aのいずれかのみを返す
    mult: float
        大きい程,離れていても結合していると判定される.
    kwargs:
        | 元素毎で解離の共有結合半径を指定できる
        | ex) H=0.5 で水素の共有結合半径を0.5Åに変更できる.
    """
    cutoff = natural_cutoffs(atoms, mult=mult,**kwargs)
    nl = build_neighbor_list(atoms,cutoff)
    ana = Analysis(atoms, nl=nl)
    if unique:
        bonds = ana.unique_bonds[0]
    else:
        bonds = ana.all_bonds[0]
    return bonds


def connected_components(atoms,indices=None,mult=1.0,**kwargs):
    """atoms中に存在する分子をindex番号毎にまとめたジェネレーターを返す
    
    Exampleを参照
    
    Parameters:
    
    atoms: Atoms
        対象のAtoms
    indices: list of int
        特定の原子のみcomponentを調べる時に指定する.
        Noneの場合,全ての原子を対象とする.
    mult: float
        大きい程,離れていても結合していると判定される.
    kwargs:
        | 元素毎で解離の共有結合半径を指定できる
        | ex) H=0.5 で水素の共有結合半径を0.5Åに変更できる.
        
        
    Examples:
    
        例えばatomsの[5,6,7,8,9]がメタン分子,[10,11]がNO分子だった場合,
        分子ごとに分離することができる
        
        >>> gen = connected_components(atoms, indices=[5,6,7,8,9,10,11])
        >>> for i in gen:
        >>>     print(i)
        >>> #--> {5,6,7,8,9}
        >>> #--> {10,11}
    """
    if indices is None:
        indices = [i for i in range(len(atoms))]
    G = nx.Graph()
    bonds = get_bond(atoms,mult=mult,**kwargs)
    for i in indices:
        for j in bonds[i]:
            if j in indices:
                G.add_edge(i, j)
    return nx.connected_components(G)

def get_unique_connections(connections):
    """ユニークなCONNECTIONSを返す
    
    同じEQを繋ぐconneciotns,同じ組み合わせのconnections,'??'などを含むconnectionsを削除し
    新しいconnectionsのリストを作成する.
    
    戻り値は2要素のタプルで1要素目がユニークなconnections, 2要素目がそれが元々どのindex番号のconnectionだったのかを返す.
    
    Parameters:
    
    connections: list of lists of int
        | [[0,1],[1,0],[2,5],["??",2],[3,4]]のようなリスト
        | grrmpy.io.read_connecitonsで取得できる

    Return:
        tuple : ユニークなconnections, 元々のindex番号
    """
    unique = []
    original_idx = []
    for i,(ini,fin) in enumerate(connections):
        if isinstance(ini,str) or isinstance(fin,str):
            continue
        if ini == fin:
            continue
        if [ini,fin] in unique or [fin,ini] in unique:
            continue
        unique.append([ini,fin])
        original_idx.append(i)
    return unique,original_idx

def partially_overlapping_groupings(iterable,n,i,strict=False):
    iterable = iter(iterable)
    l = [next(iterable) for _ in range(n)]
    yield tuple(l)
    while True:
        try:
            li = []
            for _ in range(i):
                li.append(next(iterable))
            l = l[i:] + li
            yield tuple(l)
        except StopIteration:
            if not strict and li!=[] :
                yield tuple(li)
            break