from ase import Atoms
from ase.io import read

"""'...END'を伴うパラメータ"""
WITH_END = ["Add Interaction","Bond Condition","PriorityPath"]
SECTION = ["External Atoms","Frozen Atoms","Product","Reactant"]

def read_comfile(comfile):
    with open(comfile,"r") as f:
        comtext = f.readlines()
    section_idx = []
    for i,text in enumerate(comtext):
        options_list = []
        if "%" in text:
            options_list.append(text)
        elif "Options" in text:
            options_idx = i
            section_idx.append(i)
        elif "#" in text:
            job_type, method = _read_jobtype_method(text)
            charge,multiplicity = [int(j) for j in comtext[i+2].split()]
            section_idx.append(i+3)
        elif text.strip() in SECTION:
            section_idx.append(i)
    options_list += comtext[options_idx+1:]
    options_dict = _options_to_dict(options_list)
    atoms_dict = _coordinations_to_dict(comtext, section_idx)
    return job_type, method, charge, multiplicity, atoms_dict, options_dict
        
def _read_jobtype_method(text):
    text = text.strip()
    text = text[1:].split("/")
    job_type = text[0].strip()
    if len(text) >=2:
        method = "/".join(text[1:])
    else:
        method = None
    return job_type, method

def _options_to_dict(option_text:list):
    def has_with_end(text):
        for param in WITH_END:
            if param in text:
                return True
        return False
    def convert_option_value(value):
        if not value.isdecimal():
            try:
                float(value)
            except ValueError:
                return value
        else:
            return int(value)
    options = {}
    n = 0
    while n < (len(option_text)):
        if "=" in option_text[n]:
            text = option_text[n].split("=")
            options[text[0]] = convert_option_value(text[1].strip())
            n += 1
        elif has_with_end(option_text[n]):
            key = option_text[n].strip()
            options[key] = {}
            n += 1
            while True: 
                if "=" in option_text[n]:
                    text = option_text[n].split("=")
                    options[key][text[0]] =convert_option_value(text[1].strip())
                    n += 1
                elif "END" in option_text[n]:
                    n += 1
                    break
                else:
                   options[key][option_text[n].strip()] = ""
                   n += 1
        else:
            options[option_text[n].strip()] = ""
            n += 1        
    return options


def _coordinations_to_dict(comtext, section_idx):
    atoms_dict = dict(zip(SECTION,[Atoms() for _ in range(len(SECTION))]))
    key = "NORMAL"
    atoms_dict[key] = _position_to_atoms(comtext[section_idx[0]:section_idx[1]])
    for i,s in enumerate(section_idx[1:-1],start=2):
        key = comtext[s].strip()
        atoms_dict[key] = _position_to_atoms(comtext[s+1:section_idx[i]])
    return atoms_dict

def _position_to_atoms(positions_text):
    chemical_symbols = []
    positions = []
    for p in positions_text:
        pos = p.split()
        chemical_symbols.append(pos[0])
        positions.append([float(i) for i in pos[1:4]])
    return Atoms(chemical_symbols,positions)

def frozen2atoms(comfile,poscar=None):
    """COMファイルのFrozenAtomsをAtomsオブジェクトに変換する

    Parameters:
    
    comfile: str or Path
        | comファイルパス.
    poscar: str or Path
        | POSCARパスを設定した場合.セル情報を読み取りAtomsオブジェクトに適用する.
        
    Returns:
        Atoms: FrozenAtomsのAtomsオブジェクト
    """
    with open(comfile,"r") as f:
        comtext = f.readlines()
        
    for i,text in enumerate(comtext):
        if "Frozen Atoms" in text:      
            frozen_idx = i
        elif "Options" in text:
            option_idx = i
            break
    frozrn_atoms = _position_to_atoms(comtext[frozen_idx+1:option_idx])
    
    if poscar:
        pos_atoms = read(poscar,format="vasp")
        cell = pos_atoms.get_cell()
        frozrn_atoms.set_cell(cell)
        frozrn_atoms.set_pbc(True)    
    return frozrn_atoms