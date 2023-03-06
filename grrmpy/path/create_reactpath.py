from configparser import ConfigParser, ExtendedInterpolation, BasicInterpolation
import json
import re
from pathlib import Path
from ase.io import read
from grrmpy.path.reaction_path import ReactPath,ReactPaths

def create_reactpath(file):
    if not Path(file).exists():
        raise FileNotFoundError(f"{file}が存在しません")
    # parserの設定
    config  = ConfigParser(interpolation= ExtendedInterpolation())
    config.optionxform = str # キーが小文字になることを防ぐ
    config.read(file, encoding='utf-8')
    # セクション
    title_list = config.sections()
    if "Dir" in config.sections():
        title_list.remove("Dir") # ["ReactPath1","ReactPath2"]
        replace_dict = dict(config.items('Dir')) #{'EQ': './EQ_List/', 'TS': './TS_List/'}
    path_txt_list = [config[title]["path"] for title in title_list]
    solid_dot_list = [[i.strip() for i in re.split("(;)|(~)|(\|)", path_txt) if not i is None]
                  for path_txt in path_txt_list]
    solid_list = [[j for i,j in enumerate(solid_dot) if i%2==0] for solid_dot in solid_dot_list]
    dot_list = [[j for i,j in enumerate(solid_dot) if i%2==1 ] for solid_dot in solid_dot_list]

    name_list = []
    solid_path_list = []
    for title,solid in zip(title_list,solid_list):
        solid_path = []
        if "name" in config[title].keys():
            solid_path =[species+".traj" for species in solid]
            name_list.append(json.loads(config.get(title,"name")))
        elif "Dir" in config.sections():
            for species in solid:
                for p,path in replace_dict.items():
                    if p in species:
                        solid_path.append(species.replace(p,replace_dict[p])+".traj")
            name_list.append(solid)
        else:
            solid_path =[species+".traj" for species in solid]
            name_list.append(solid)
        solid_path_list.append(solid_path)         

    atoms_list = [[read(path) for path in solid_path] for solid_path in solid_path_list]

            
    barrier_idx_list = []
    mask_list = []
    for dot in dot_list:
        barrier_idx = []
        mask = []
        for i,sign in enumerate(dot):
            if "~" in sign:
                barrier_idx.append(i)
            if "|" in sign:
                mask.append(False)
            else:
                mask.append(True)
        barrier_idx_list.append(barrier_idx)
        mask_list.append(mask)
    
    data_list = [{"atoms":atoms,"name":name} for atoms,name in zip(atoms_list,name_list)]
    
    react_path_list = [ReactPath(data,mask=mask,barrier_less_index=barrier_idx,title=title)
                      for data,title,mask,barrier_idx 
                       in zip(data_list,title_list,mask_list,barrier_idx_list)]
    
    if len(react_path_list) == 1:
        return react_path_list[0]
    else:
        return ReactPaths(react_path_list)

