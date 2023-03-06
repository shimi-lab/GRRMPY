import xlsxwriter
from xlsxwriter.utility import xl_cell_to_rowcol
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go


deffault_x_axis = {
    'visible':False,
    'name': "Reaction Coordinate", #ラベル名設定
    'num_font': {'name': 'Arial', 'size': 14}, #数値のフォント設定
    'name_font':{'name': 'Arial','size': 14,'bold': False,'italic': False}, #ラベルのフォント設定
    'min': 0,
    'crossing' : -10000, #軸の交点設定
    'major_tick_mark' : 'inside', #主目盛の向き設定
    'line': {'color': "black", 'width': 1.5} #軸の色と線の太さ設定
    }

deffault_y_axis = {
    'name' : "Energy (kJ/mol)", 
    'num_font': {'name': 'Arial', 'size': 14},
    'name_font': {'name': 'Arial','size': 14,'bold': False,'italic': False},
    'crossing' : -100,
    'major_gridlines': {'visible': False}, 
    'major_tick_mark' : 'inside', 
    'minor_tick_mark': 'inside', 
    'line': {'color': "black", 'width': 1.5}
    }

deffault_legend = {'none': True}
deffault_chartarea = {'border': {'none': True}}
deffault_size = {}
deffault_title = {}
deffault_plotarea = {}


def to_excell(outfile,
              table_data,
              solid_data,
              dot_data,
              pos="D2",
              sheet_name = None,
              eq_color = "black",
              ts_color = "red",
              pt_color = "green",
              dot_color = "black",
              dash_type = 'round_dot',
              solid_width = 2.5,
              dot_width = 2,
              graph_pos = "F3",
              x_offset = 25,
              y_offset = 10,
              x_axis = None,
              y_axis = None,
              title = None,
              legend = None,
              chartarea = None,
              plotarea = None,
              size = None,): 
    """エクセルにデータとグラフを挿入する

    Parameters:
    
    outfile (_type_): _description_
    table_data (_type_): _description_
    graph_data (_type_): _description_
    pos (str, optional): _description_. Defaults to "D2".
    """
    if x_axis is None:
        x_axis = deffault_x_axis
    else:
        x_axis = {**deffault_x_axis, **x_axis}

    if y_axis is None:
        y_axis = deffault_y_axis
    else:
        y_axis = {**deffault_y_axis, **y_axis}
        
    if legend is None:
        legend = deffault_legend
    else:
        legend = {**deffault_legend,**legend}
        
    if chartarea is None:
        chartarea = deffault_chartarea
    else:
        chartarea = {**deffault_chartarea,**chartarea}
    
    if plotarea is None:
        plotarea = deffault_plotarea
    else:
        plotarea = {**deffault_plotarea,**plotarea}
        
    if title is None:
        title = deffault_title
    else:
        title = {**deffault_title,**title}
        
    if size is None:
        size = deffault_size
    else:
        size = {**deffault_size,**size}
    
    book = xlsxwriter.Workbook(outfile)
    if sheet_name is None:
        sheet = book.add_worksheet()
    else:
        sheet = book.add_worksheet(sheet_name)
        
    ## データ書き込み
    sheet.write(0, 0, "name")
    for row_num, data in enumerate(table_data["name"],start=1):
        sheet.write(row_num, 0, data)
    sheet.write(0, 1, "energy(kJ/mol)")
    for row_num, data in enumerate(table_data["energy"],start=1):
        sheet.write(row_num, 1, data)
        
    ## グラフ作成のためのデータを書き込む
    pos = xl_cell_to_rowcol(pos)
    for i,(col_name,col) in enumerate(solid_data.items(),start=pos[1]):
        sheet.write(pos[0], i, col_name) # カラム名書き込み
        for j,data in enumerate(col,start=pos[0]+1):
            sheet.write(j, i, data)
    for i,(col_name,col) in enumerate(dot_data.items(),start=pos[1]):
        sheet.write(pos[0], i+len(dot_data.columns), col_name) # カラム名書き込み
        for j,data in enumerate(col,start=pos[0]+1):
            sheet.write(j, i+len(dot_data.columns), data)

    # グラフ作成
    n_data = len(dot_data)
    sheet_name = sheet.get_name()
    chart = book.add_chart({'type': 'scatter', 'subtype': 'straight'})
    # solid
    for i,name in enumerate(solid_data["name"],start=1):
        if "TS" in name:
            color = ts_color
        elif "PT" in name:
            color = pt_color
        else:
            color = eq_color
        chart.add_series({
            'categories': [sheet_name, pos[0]+i, pos[1], pos[0]+i, pos[1]+1],
            'name':       [sheet_name, pos[0]+i, pos[1]+2],
            'values':     [sheet_name, pos[0]+i, pos[1]+3, pos[0]+i, pos[1]+4],
            'line' : {'color': color, 'width': solid_width},
        })
    # dot
    for i in range(1,n_data+1):
        chart.add_series({
            'categories': [sheet_name, pos[0]+i, pos[1]+5, pos[0]+i, pos[1]+6],
            'name':       [sheet_name, pos[0]+i, pos[1]+7],
            'values':     [sheet_name, pos[0]+i, pos[1]+8, pos[0]+i, pos[1]+9],
            'line' : {'color': dot_color, 'width': dot_width, 'dash_type': dash_type},
        })
    sheet.insert_chart(graph_pos, chart, {'x_offset': x_offset, 'y_offset': y_offset})

    #グラフエリアのサイズ
    chart.set_size(size)
    # 凡例
    chart.set_legend(legend)
    # タイトル
    chart.set_title(title)
    #X軸の設定
    chart.set_x_axis(x_axis)
    
    #Y軸の設定
    chart.set_y_axis(y_axis)
    
    #グラフエリアの設定
    chart.set_chartarea(chartarea)
    # プロットエリアの設定
    chart.set_plotarea(plotarea)
    
    book.close()
    
def to_fig(solid_data,dot_data,xlabel="Energy",ylabel="Reaction Coordinate"):
    """matplotlibのfigに変換する"""
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    for _,solid_xini,solid_xfin,name,solid_yini,solid_yfin in solid_data.itertuples():
        x = [solid_xini,solid_xfin]
        y = [solid_yini,solid_yfin]
        if "TS" in name:
            color = "red"
        elif "PT" in name:
            color = "green"
        else:
            color = "black"
        ax.plot(x,y, color=color, linestyle='-',linewidth=2)
        
    for _,dot_xini,dot_xfin,dot_name,dot_yini,dot_yfin in dot_data.itertuples():
        x = [dot_xini,dot_xfin]
        y = [dot_yini,dot_yfin]
        ax.plot(x,y, color="black", linestyle='dotted',linewidth=1.5)
        
    ax.set_ylabel(xlabel)
    ax.set_xlabel(ylabel)
    return fig


def to_plotly(solid_data,dot_data,annotation=[],title_list=None):
    """plotlyのfigに変換する

        Property:
        
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
    """
    def find_y_value(x):
        for _,(dot_xini,dot_xfin,dot_name,dot_yini,dot_yfin) in dot_data.iterrows():
            if x > dot_xini and x < dot_xfin:
                y = (x-dot_xini)*(dot_yfin-dot_yini)/(dot_xfin-dot_xini)+dot_yini
                return y
        for i,(solid_xini,solid_xfin,name,solid_yini,solid_yfin) in solid_data.iterrows():
            if x > solid_xini and x < solid_xfin:
                return solid_yini
        
    dot_plot = [
        go.Scatter(
            x=[dot_xini,dot_xfin], y=[dot_yini, dot_yfin],
            name=dot_name,
            line_dash="dot",
            hoverinfo='none',
            mode = "lines",
            line=dict(color="black")
        ) for _,(dot_xini,dot_xfin,dot_name,dot_yini,dot_yfin) in dot_data.iterrows()
    ]

    color_list = ["red" if "TS" in text else "green" if "PT" in text else "black" for text in solid_data["name"]]
    
    if title_list is None:
        solid_plot = [
            go.Scatter(
                x=[solid_xini,solid_xfin], y=[solid_yini, solid_yfin],
                name=name,
                line_dash="solid",
                mode = "lines",
                hovertemplate=f'{name}: {round(solid_yini,2)}',
                line=dict(color=color_list[i],width=3.5),
            ) for i,(solid_xini,solid_xfin,name,solid_yini,solid_yfin) in solid_data.iterrows()
        ]
    else:
        solid_plot = [
            go.Scatter(
                x=[solid_xini,solid_xfin], y=[solid_yini, solid_yfin],
                name=name,
                line_dash="solid",
                mode = "lines",
                hovertemplate=f'{name}: {round(solid_yini,2)} ({title})',
                line=dict(color=color_list[i],width=3.5),
            ) for (i,(solid_xini,solid_xfin,name,solid_yini,solid_yfin)),title in zip(solid_data.iterrows(),title_list)
        ]

    fig = go.Figure(data=dot_plot+solid_plot)
    
    for position,*args in annotation:
        if type(position)==tuple or type(position)==tuple:
            x,y = position
        else:
            x = position
            y = find_y_value(x)
        if len(args) == 1:
            text = args[0]
            kwargs = {}
        else:
            text,kwargs = args
        fig.add_annotation(x=x, y=y, text=text, **kwargs)
    
    fig.update_layout(showlegend=False,
                      xaxis=dict(title='Reaction Coordinate',
                                 range=(min(solid_data["solid_xini"])-1,
                                        max(solid_data["solid_xfin"])+1)),
                      yaxis=dict(title='Energy (kJ/mol)',),
                     )
    return fig

def to_html(solid_data,dot_data,full_html=False,include_plotlyjs="cdn",annotation=[],title_list=None,**kwargs):
    """htmlテキストに変換する"""
    fig = to_plotly(solid_data,dot_data,annotation=annotation,title_list=title_list)
    return fig.to_html(full_html=full_html,include_plotlyjs=include_plotlyjs,**kwargs)