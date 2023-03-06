from ase.neb import NEB
from ase import Atoms
import ase.units as units
from ase.vibrations import Vibrations
from pathlib import Path

#User
from grrmpy.calculator import pfp_calculator
from grrmpy.vibrations.functions import to_html_table_and_imode,to_html_graph
from grrmpy.neb.functions import to_html_nebgraph
from grrmpy.functions import to_html_energy_diagram

def write_html(html,obj,append=False,*args,**kwargs):
    """htmlファイルとして保存する
    
    Parameters:
    
    html: string
        保存ファイル名
    obj: object
        | objの引数には次のオブジェクトを指定できる
        | - NEB
        |     NEBのエネルギーダイアグラムを作成
        | - Vibrations
        |     エネルギーダイアグラムと振動数表を作成
        | - trajファイル
        |     エネルギーダイアグラムを作成
    """
    if type(obj) == NEB:
        write_neb_graph(html,obj, calc_func=pfp_calculator,append=append,**kwargs)
    elif type(obj) == Vibrations:
        write_vibtb_and_vibgraph(html,obj, calc_func=pfp_calculator,append=append,**kwargs)
    elif type(obj) == Path or type(obj) == str:
        write_vib_graph(html,obj,calc_func=pfp_calculator,append=append,**kwargs)
    elif type(obj) == list:
        if type(obj[0]) == Atoms:
            write_energy_diagram(html,obj,calc_func=pfp_calculator,append=append,**kwargs)
    
def write_vib_table(
    html:str,
    vib_obj,
    full_html:bool=True,
    append=False,):
    """振動数の表をhtmlのstrを出力する&虚振動の振動モード番号を出力する

    Parameters:
    
    html: string or None
        | 出力するhtmlファイル名.
        | Noneの場合,htmlテキストを戻り値として出力する
    vib_obj: Viblations object
        Viblationsオブジェクト
    full_html: boolean
        | <html>タグから始まる,完全なhtmlを出力する場合True
        | Falseの場合<div>タグから始まるテキストを出力
    append:bool
        Trueの場合追記モード
    """
    full_html = False if append else full_html
    mode = "a" if append else "w"
    html_text,_ = to_html_table_and_imode(vib_obj,full_html=full_html,include_plotlyjs="cdn")
    if not html:
        return html_text        
    with open(html,mode) as f:
        f.write(html_text)
        
def write_vib_graph(
    html,vib_obj,
    n:int,
    outfile=None,
    calc_func=pfp_calculator, 
    kT=units.kB * 300,
    nimages=30,
    full_html=True,
    append=False,):
    """エネルギーダイアグラムを作成する

    Parameters:
    
    html: string or None
        | 出力するhtmlファイル名.
        | Noneの場合,htmlテキストを戻り値として出力する
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
    append: bool
        追記する場合True
    """
    full_html = False if append else full_html
    mode = "a" if append else "w"
    html_text = to_html_graph(
        vib_obj,
        n,
        outfile,
        calc_func,
        kT,
        nimages,
        full_html,
        include_plotlyjs="cdn")
    if not html:
        return html_text  
    with open(html,mode) as f:
        f.write(html_text)
        

def write_vibtb_and_vibgraph(
    html,
    vib_obj,
    outfile=None,
    calc_func=pfp_calculator,
    kT=units.kB * 300,
    nimages=30,
    append=False):
    """振動数の表とエネルギーダイアグラムのhtmlテキストを出力する.

    Parameters:
    
    html: string or None
        | 出力するhtmlファイル名.
        | Noneの場合,htmlテキストを戻り値として出力する
    vib_obj: Viblations object
        Viblationsオブジェクト
    outfile: str
        | 振動の構造をtrajで出力する場合,trajファイル名.
        | Noneの場合出力しない
        | (write_mode()でtrajファイルを出力する)
    calc_func: function object
        Claculatorを返す関数
    nimages: integer
        | イメージの数, デフォルトは30
        | 偶数で設定する事を推奨する
    append: bool
        追記の場合True.
    Note:
        | 虚振動がないor複数ある場合は,エネルギーダイアグラムは出力しない
    """
    mode = "a" if append else "w"
    tb_txt,imode = to_html_table_and_imode(vib_obj,full_html=False,include_plotlyjs="cdn")
    if imode:
        fig_txt = write_vib_graph(None,vib_obj,imode,outfile,calc_func,kT,nimages,full_html=False)
    else:
        fig_txt = ""
    if not html:
        return tb_txt+fig_txt
    with open(html,mode) as f:
        f.write(tb_txt+fig_txt)
        
def write_neb_graph(html,
                    neb_obj,
                    append=False,
                    title="NEB Energy Diagram",
                    calc_func=pfp_calculator,
                    full_html=True,
                    unit="kJ/mol",
                    highlight=None,
                    annotation=[],
                    **kwargs):
    full_html = False if append else full_html
    mode = "a" if append else "w"
    html_txt = to_html_nebgraph(
        neb_obj,
        calc_func=calc_func,
        full_html=full_html,
        unit=unit,
        include_plotlyjs="cdn",
        highlight=highlight,
        annotation=annotation,
        **kwargs)
    if not html:
        return html_txt
    with open(html,mode) as f:
        f.write(html_txt)
        
def write_energy_diagram(html,
                         images,
                         append=False,
                         title="Energy Diagram",
                         calc_func=pfp_calculator,
                         full_html=True,
                         unit="kJ/mol",
                         highlight=None,
                         annotation=[],
                         **kwargs):
    full_html = False if append else full_html
    mode = "a" if append else "w"
    html_txt = to_html_energy_diagram(
        images,
        calc_func=calc_func,
        full_html=full_html,
        unit=unit,
        title=title,
        include_plotlyjs="cdn",
        highlight=highlight,
        annotation=annotation,
        **kwargs)
    if not html:
        return html_txt
    with open(html,mode) as f:
        f.write(html_txt)
        
