from configparser import ConfigParser, ExtendedInterpolation, BasicInterpolation
import json
import re
from pathlib import Path
from ase.io import read
from grrmpy.path.reaction_path import ReactPath,ReactPaths

def pathtxt2connection_list(file="Path.dat"):
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