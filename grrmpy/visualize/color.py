from pathlib import Path

# USER
from grrmpy.io.read_elementsini import read_elementsini,read_csv


default= read_elementsini(Path(__file__).joinpath('../../', 'default_color.ini').resolve())
vesta= read_elementsini(Path(__file__).joinpath('../../', 'vesta_color.ini').resolve())
jmol = read_csv(Path(__file__).joinpath('../../', 'jmol_color.csv').resolve())


# grrmpy.visualize.functions.generate_js_code を参照するとよい
# grrmpy.io.read_elementsini.read_elementsini を参照するとよい