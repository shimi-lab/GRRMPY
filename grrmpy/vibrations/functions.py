import pandas as pd
import plotly.graph_objects as go
import plotly.io as pyi
import ase.units as units
from ase.io import write
from scipy.signal import argrelmax,argrelmin
import numpy as np
from ase.units import mol, kJ

# User Modules
from grrmpy.calculator import pfp_calculator

def get_vibdf(vib_obj):
    """vib.summary()の結果をDataFrameで取得する

    Parameters:
    
    vib_obj: Vibrarions object
        Vibrarionsオブジェクト()
    """
    vib_data = vib_obj.get_vibrations()
    table = vib_data.tabulate()
    data = [{"meV":row.split()[1],"cm^-1":row.split()[2]} for row in table.split("\n")[3:-3]]
    df = pd.DataFrame(data)
    return df

def find_ivib(vib_obj):
    """

    Parameters:
    
    vib_obj: Viblations object
        Viblationsオブジェクト
        
    Returns:
        DataFrame: 虚振動のみを抽出したDataFrame
    """
    df = get_vibdf(vib_obj)
    df = df[df['cm^-1'].str.contains('i')] # 虚振動が含まれるデータのみを抽出
    return df
    
def get_imode(vib_obj):
    """虚振動が1つだった場合,虚振動の振動モード番号を返す

    Parameters:
    
    vib_obj: Vibrarions object
        Vibrationsオブジェクト
        
    Returns:
        int: 振動モード番号
    """
    df = find_ivib(vib_obj) # 虚振動部分のDataFrameを抽出する
    n = [i for i in df.index.tolist() if i not in [1,2,3,4]]# 並行移動のモード1,2,3,4が虚振動の場合は排除
    if len(n) == 1: # 虚振動が1つだった場合
        return n[0]
    else:
        None
        
def get_vib_images(vib_obj,n, outfile=None, kT=units.kB * 300, nimages=30):
    """IRC計算に必要なini,fin構造を取得する
    
    Parameters:
    
    vib_obj: Viblations object
        Viblationsオブジェクト
    n: int
        振動モード番号
    outfile: str
        trajファイルを出力する場合にはファイル名を設定
    nimages: int
        Imageの数
    """
    n %= len(vib_obj.get_energies())
    vib_images = [image 
                  for image 
                  in vib_obj.get_vibrations().iter_animated_mode(n,temperature=kT, frames=nimages)]    
    if outfile:
        """trajファイルに保存する"""
        write(outfile,vib_images)
    return vib_images
    
def has_ivib(vib_obj):
    """たった1つの虚振動があるかを調べる

    Parameters:
    
    vib_obj: Viblations object
        Viblationsオブジェクト

    Returns:
        bool: 虚振動が1つの時,True
    """
    n = get_imode(vib_obj)
    if not n is None:
        return True
    else:
        return False
    
def _to_table_and_imode(vib_obj):
    vib_data = vib_obj.get_vibrations()
    table = vib_data.tabulate()
    data = [
        {"<b>#</b>":row.split()[0],
         "<b>meV</b>":row.split()[1],
         "<b>cm^-1</b>":row.split()[2]} 
        for row in table.split("\n")[3:-3]
        ]
    df = pd.DataFrame(data)
    idf = df[df['<b>cm^-1</b>'].str.contains('i')] # 虚振動が含まれるデータのみを抽出
    n = [i for i in idf.index.tolist() if i not in [1,2,3,4]]
    if len(n) == 1: #虚振動が1つだげ存在するとき
        fill_color = [["rgb(233,235,245)" if i!=n[0] else "rgb(255,155,155)" for i in range(len(df))]]*3
        n_mode = n[0]
    else:
        fill_color = [["rgb(233,235,245)" for _ in range(len(df))]]*3
        n_mode = None     
    table = go.Table(
            header=dict(values=list(df.columns),
                    fill = {"color":"rgb(68,114,196)"},
                    font = {"color":"rgb(255,255,255)","family" : "Arial"},
                    align='center'),
        cells=dict(values=[df["<b>#</b>"], df["<b>meV</b>"], df["<b>cm^-1</b>"]],
                    fill = {"color":"rgb(233,235,245)"},
                    fill_color = fill_color,
                    font = {"family":"Arial"},
                    align='center')
        ) 
    return table, n_mode

def _to_fig_table_and_imode(vib_obj):
    table,n_mode = _to_table_and_imode(vib_obj)
    fig = go.Figure()
    fig.add_trace(table)
    return fig, n_mode
    
def to_html_table_and_imode(vib_obj,full_html:bool=False,**kwargs):
    """振動数の表をhtmlのstrを出力する&虚振動の振動モード番号を出力する

    Parameters:
    
    vib_obj: Viblations object
        Viblationsオブジェクト
        
    Returns:
    
    html_string: str
        <div>タグから始まるhtml文字列
    n: int
        虚振動の振動モード番号, 虚振動が存在しない場合,または複数ある場合はNone
    """
    fig,n_mode = _to_fig_table_and_imode(vib_obj)
    return pyi.to_html(fig,full_html=full_html,**kwargs), n_mode

def to_fig_table(vib_obj):
    fig,_ = _to_fig_table_and_imode(vib_obj)
    return fig

def to_html_table(vib_obj,full_html=False):
    html_txt,_ = to_html_table_and_imode(vib_obj,full_html)
    return html_txt

def to_fig_graph(vib_obj, n:int, outfile=None, calc_func=pfp_calculator, kT=units.kB * 300, nimages=30):
    n %= len(vib_obj.get_energies())
    vib_images = [image 
                  for image 
                  in vib_obj.get_vibrations().iter_animated_mode(n,temperature=kT, frames=nimages)]
    for image in vib_images:
        image.calc = calc_func()
    
    if outfile:
        """trajファイルに保存する"""
        write(outfile,vib_images)
    x = [i for i in range(len(vib_images))]
    y = [atoms.get_potential_energy() for atoms in vib_images]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=y,
                   line=dict(width=2,color="black"))
    )
    return fig
    
def to_html_graph(vib_obj, n:int, outfile=None, calc_func=pfp_calculator, kT=units.kB * 300, nimages=30,full_html=False,**kwargs):
    """振動エネルギーダイアグラムのhtmlテキストを出力
    
    outfileを設定した場合にはtrajファイルを作成する
    
    Parameters:
    
    vib_obj: Viblations object
        Viblationsオブジェクト
    n: integer
        振動モード番号
    outfile: str
        | 振動の構造をtrajで出力する場合,trajファイル名.
        | Noneの場合出力しない
        | (write_mode()でtrajファイルを出力する)
    calc_func: function object
        Claculatorを返す関数
    nimages: integer
        | イメージの数, デフォルトは30
        | 偶数で設定する事を推奨する
    full_html:
        | <html>タグから始まる,完全なhtmlを出力する場合True
        | Falseの場合<div>タグから始まるテキストを出力
    """
    fig = to_fig_graph(vib_obj, n, outfile, calc_func, kT, nimages)
    return pyi.to_html(fig,full_html=full_html,**kwargs)


def find_ts_idx(vib_images,dif=2,calc_func=pfp_calculator,**kwargs):
    """TSのindex番号とNewton法などで計算すべきか調べる.
    
    | vib_imagesの中心辺り(vib_imagesを5等分したうち真ん中)の極大値をTSと判断する.
    | 極大値が複数ある場合はより中心に近い点をTSとする.
    | 極大値が存在しない場合は中心の構造をTSと判断する
    | 
    | 極大値と極小値間のエネルギー差が小さい場合,Newton法(FIRE)で計算すべきと判断しTrueを返す.
    | Falseの場合,容易に収束できる可能性が高いのでBFGSなどでIRC計算しても問題ない.

    Parameters:
    
    vib_images: list of Atoms
        振動のAtomsリスト.
    dif:float
        TSと極小値の差がdif(kJ/mol)以下の時,Newton法で計算すべきと判断する.デフォルトは2kJ/mol
    calc_func: function
        | clalulatorを返す関数,事前にvib_imagesにclaculatorを付けている場合は設定は不要.
        | デフォルトはPFPのCalculator

    Returns:
    
    2つの要素のタプルで返す
    
    ts_idx: integer or boolean
        | TSと判断したindex番号.
        | info(下記)の3番目の要素が1以外の場合は得られたts_idxはあまり信頼できないので注意.
    info: tuple of 3 elements
        - 1要素目:
            Reverse計算でNewton法を用いるべきと判断された場合True.
        - 2要素目:
            forward計算でNewton法を用いるべきと判断された場合True.
        - 3要素目:
            極大値の数.(1の場合うまくTSの場所を判断できたと考える.0の場合はTSはなかったと考える)
    """
    nimages = len(vib_images)
    for image in vib_images:
        if not image.get_calculator():
            image.calc = calc_func()
    
    energy = [i.get_potential_energy() for i in vib_images]
    middle_img = list(np.array_split(energy, 5))[2] # TS付近の構造(5等分した内お3番目)
    middle_ini_idx = energy.index(middle_img[0])
    middle_fin_idx = energy.index(middle_img[-1])
    peak_idx = argrelmax(middle_img)
    
    if len(peak_idx) == 0:
        """極大値が見つからなかった場合"""
        ts_idx = int(nimages/2)
    elif len(peak_idx) == 1:
        """極大値が1つ見つかった時(理想的)"""
        ts_idx = energy.index(middle_img[peak_idx])
    else:
        """極大値が複数見つかった時"""
        ts_idx = idx_of_the_nearest(peak_idx,int(len(middle_img)/2))
        
    ini = min(energy[middle_ini_idx:ts_idx+1])
    fin = min(energy[ts_idx:middle_fin_idx+1])
    ts = energy[ts_idx]
    
    if abs((ts-ini)*mol/kJ)<dif:
        reverse = True
    else:
        reverse = False
        
    if abs((ts-fin)*mol/kJ)<dif:
        forward = True
    else:
        forward = False
    return (ts_idx,(reverse,forward,len(peak_idx)))

def idx_of_the_nearest(data, value):
    """data(1次元リスト)から最もvalueに近いindex番号を返す"""
    idx = np.argmin(np.abs(np.array(data) - value))
    return idx