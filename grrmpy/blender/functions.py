import pandas as pd

def rgba2hex(rbga):
    R,G,B,_ = [i for i in rbga] # 1で正規化されたrgba
    r = format(hex(int(255*R)).replace('0x', ''), '0>2')
    g = format(hex(int(255*G)).replace('0x', ''), '0>2')
    b = format(hex(int(255*B)).replace('0x', ''), '0>2') 
    return "#"+r+g+b

def hex2rgba(color_code,a=1.0):
    R = int(color_code[1:3], 16)
    G = int(color_code[3:5], 16)
    B = int(color_code[5:7], 16)
    return (R,G,B,a)

def read_elementsini(elements_ini_path):
    """構造をViewerで可視化する際の原子の色を規定する辞書を作成する.
    
    Parameters:
    
    elements_ini_path: Union[str, Path]
        | elements.ini(VESTAで使われるやつ)のパス
    Returns:
        dict: 元素名,RGBAの辞書.{'H':(0.2,0.3,0.4,1)}
    """     
    # iniファイルの読み込み   
    df = pd.read_table(elements_ini_path,header=None,engine='python',delim_whitespace=True)
    element = df[1].tolist()
    r = df[5].tolist()
    g = df[6].tolist()
    b = df[7].tolist()
    color_dict = {e:(R,G,B,1) for e,R,G,B in zip(element,r,g,b)}
    return color_dict

def read_csv(filename):
    """原子の色を規定したcsvファイルを読み込みcolor_dictを返す

    Parameters:
    
    filename: str 
        | csvファイル
        
    Returns:
        dict: 元素名,RGBAの辞書.{'H':(0.2,0.3,0.4,1)}
    """
    df = pd.read_csv(filename)
    element = df.iloc[:,0].tolist()
    r = df.iloc[:,1].tolist()
    g = df.iloc[:,2].tolist()
    b = df.iloc[:,3].tolist()
    color_dict = {e:(R/255,G/255,B/255,1) for e,R,G,B in zip(element,r,g,b)}
    return color_dict

def parsestr2list(num_str):
    """
    
    ハイフン，コンマで指定された文字列を数字のリストに変換する.
    
    Parameters:
    
    num_str: str
        '0-10,23,25-26'のように指定する.
        
    Return: list of int
        数字のリストとして返す
    """
    num = []
    for part in num_str.split(','):
        p1 = part.split('-')
        if len(p1) == 1:
            num.append(int(p1[0]))
        else:
            num.extend(list(range(int(p1[0]),int(p1[1])+1)))
    return num

def get_unique_items(l):
    l = list(reversed(l))
    return_list = []
    sum_list = set([])
    for i in l:
        s = set(i)
        dif = s-sum_list
        return_list.append(sorted(list(dif)))
        sum_list |= s
    return_list = list(reversed(return_list))
    return return_list