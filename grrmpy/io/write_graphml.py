import networkx as nx
import pandas as pd
from pandas import DataFrame

def write_graphml(savefile:str,eq_df:DataFrame,ts_df:DataFrame,pt_df:DataFrame=None):
    """Graphmlを作成する

    Parameters:

    savefile: str or path obj
        保存ファイル名(.graphml)
    eq_df : DataFrame
        | DataFrame. GrrmData.eq.summaryで取得可能.
        | 次のカラムを有している必要がある
        | "node","name","group","E/Hartree",E/kJmol-1",(必要に応じて"USE")
    ts_df : DataFrame
        | DataFrame. GrrmData.ts.summaryで取得可能.
        | 次のカラムを有している必要がある
        | "edge","name","source","target","E/Hartree",E/kJmol-1","forward/kJmol-1","reverse/kJmol-1",(必要に応じて"USE")
    pt_df : DataFrame
        | DataFrame. GrrmData.pt.summaryで取得可能.
        | ts_dfと同様のカラムを有している必要がある
    """
    G = nx.Graph()
    node_data = _make_node_data(eq_df)
    G.add_nodes_from(node_data)
    if pt_df is not None:
        pt_data = _make_edge_data(pt_df)
        G.add_edges_from(pt_data)
    ts_data = _make_edge_data(ts_df)
    G.add_edges_from(ts_data)
    nx.write_graphml_lxml(G,savefile)
    return G

def _make_node_data(eq_df):
    nodes = eq_df["node"]
    names = eq_df["name"]
    groups = eq_df["group"]
    energies = eq_df["E/Hartree"]
    jules_es = eq_df["E/kJmol-1"]
    iter_objs = zip(nodes,names,groups,energies,jules_es)
    if "USE" in eq_df.columns:
        uses = eq_df["USE"]
    else:
        uses = [True for _ in range(len(eq_df))]
    iter_objs = zip(nodes,names,groups,energies,jules_es,uses)
    node_data = [(node,
                    {"name":name,
                    "id":node,
                    "group":f"{group}",
                    "E/Hartree":energy,
                    "E/kJmol-1":jules_e,
                    "USE":use}
                    ) for node,name,group,energy,jules_e,use in iter_objs
                    ]
    return node_data

def _make_edge_data(df):
    df = df.sort_values('E/Hartree') #昇順でソート
    # 欠損地の削除(無意味なTS,??,DCなどを削除)
    df[['source', 'target']] = df[['source', 'target']].apply(pd.to_numeric, errors="coerce")#??やDCなどの文字列があればNaNに変更(この時intでなくfloatになる)
    df = df.dropna() #NaNを削除
    df = df.astype({'source': 'int', 'target': 'int'}) # floatになったのでintに戻す
    df = df[df['source']!=df['target']] #self-loop(グラフ理論)削除
    df = df.reset_index(drop=True) #index番号を振り直す
    s = df.apply( lambda r: set((r['source'], r['target'])), axis=1) #同じ反応経路(parallel edges)を抽出する
    df = df.loc[s.drop_duplicates().index] #同じ反応経路(parallel edges)の重複を無くす
    names = df["name"]
    sources = df["source"]
    targets = df["target"]
    energies = df["E/Hartree"]
    jules_es = df["E/kJmol-1"]
    forwards = df["forward/kJmol-1"]
    reverses = df["reverse/kJmol-1"]
    if "USE" in df.columns:
        uses = df["USE"]
    else:
        uses = [True for _ in range(len(df))]
    iter_objs = zip(names,sources,targets,energies,jules_es,forwards,reverses,uses)
    edge_data = [(source,target,
                    {'name':name,
                    "interaction":f"{source}-{target}",
                    "E/Hartree":energy,
                    "E/kJmol-1":jules_e,
                    "forward":forward,
                    "reverse":reverse,
                    "USE":use
                    }
                    ) for name,source,target,energy,jules_e,forward,reverse,use in iter_objs
                    ]
    return edge_data