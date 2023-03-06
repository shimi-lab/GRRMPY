from grrmpy.io.read_com import read_comfile
from grrmpy.io.read_poscar import get_cell_and_pbc,get_cell,get_pbc
from copy import deepcopy
from pathlib import Path
from ase import Atoms
from ase.io import read
"""
足し算はできない
引き算をする事でOptionの違いを見ることができる
self.atoms_dict.keys() = "External Atoms","Frozen Atoms","Product","Reactant"
"""

class COM():
    """COMファイルの情報をまとめたクラス
    
    Parameters:

    comfile:str or Path obj
        XXX.comのファイルパス
    poscar: str or Path obj
        周期境界とCellを設定する場合,POSCARファイルパス
        
    Property:
    
    file: str
        comファイルパス
    options: dir
        POTIONSの情報をまとめた辞書
    job_type: str
        "SCAFIR"などジョブの種類
    method: str
        ab initio method
    charge: int
        電荷
    multiplicity: int
        多重度
    external_atoms: Atoms
        EXTERNAL ATOMS
    frozen_atoms: Atoms
        FROZEN ATOMS
    product: Atoms
        PRODUCT
    reactant: Atoms
        REACTANT
    """
    def __init__(self,comfile:str=None,poscar=None):
        self.__file = comfile
        if comfile:
            (self._job_type,
            self._method,
            self._charge,
            self._multiplicity,
            self._atoms_dict, 
            self._options) = read_comfile(comfile)
            if poscar:
                self._set_cell_and_pbc(poscar)
            
    @property
    def file(self):
        return self.__file
        
    @property
    def options(self):
        return self._options
    
    @property
    def job_type(self):
        return self._job_type
    
    @property
    def method(self):
        return self._method
    
    @property
    def charge(self):
        return self._charge
    
    @property
    def multiplicity(self):
        return self._multiplicity
    
    @property
    def external_atoms(self):
        if hasattr(self,"_atoms_dict"):
            atoms = self._atoms_dict["External Atoms"]
            if atoms:
                return atoms
        return None
    
    @property
    def frozen_atoms(self):
        if hasattr(self,"_atoms_dict"):
            atoms = self._atoms_dict["Frozen Atoms"]
            if atoms:
                return atoms
        return None
    
    @property
    def product(self):
        if hasattr(self,"_atoms_dict"):
            atoms = self._atoms_dict["Product"]
            if atoms:
                return atoms
        return None
    
    @property
    def reactant(self):
        if hasattr(self,"_atoms_dict"):
            atoms = self._atoms_dict["Reactant"]
            if atoms:
                return atoms
        return None
    
    def _set_cell_and_pbc(self,poscar):
        cell, pbc = get_cell_and_pbc(poscar)
        for atoms in self._atoms_dict.values():
            atoms.cell = cell
            atoms.pbc = pbc
             
    def set_cell(self,cell, scale_atoms=False, apply_constraint=True):
        if type(cell) == str or type(cell) == Path:
            """cell=POSCARファイルを指定した場合"""
            cell = get_cell(cell)
            for atoms in self._atoms_dict.items():# in self.atoms_listではダメ(Atoms()でなくNoneのため)
                atoms.cell = cell
        else:
            for atoms in self._atoms_dict.items(): 
                atoms.set_cell(cell,scale_atoms,apply_constraint)
                
    def set_pbc(self,pbc):
        for atoms in self._atoms_dict.items(): 
            atoms.pbc = pbc
        
    def __repr__(self) -> str:
        
        return f"{self.__class__.__name__}(com={self.file})"
    
    def __str__(self) -> str:
        return self.__repr__()
    
    def __copy__(self):
        if self.file is None:
            return __class__()
        else:
            new_obj = __class__()
            new_obj.__file = self.file
            new_obj._job_type = self.job_type
            new_obj._method = self.method
            new_obj._charge = self.charge
            new_obj._multiplicity = self.multiplicity
            new_obj._atoms_dict = deepcopy(self._atoms_dict)
            new_obj._options = deepcopy(self.options)
            return new_obj
        
    def copy(self):
        return self.__copy__()
    
    def __bool__(self):
        if self.file is None:
            return False
        return True
    
    def __eq__(self, other):
        if self.frozen_atoms == other.frozen_atoms:
            return True
        return False
    
    def __ne__(self, other: object) -> bool:
        return not self == other
    
    @classmethod
    def fromdict(cls,dct):
        new_obj = cls()
        if dct:
            job_type = dct["job_type"]
            method = dct["method"]
            charge = dct["charge"]
            multiplicity = dct["multiplicity"]
            options = dct["options"]
            atoms_dict = {key:Atoms.fromdict(val) for key,val in dct["atoms_dict"].items()}
            new_obj._job_type = job_type
            new_obj._method = method  
            new_obj._charge = charge
            new_obj._multiplicity = multiplicity
            new_obj._options = options
            new_obj._atoms_dict = atoms_dict
        return new_obj

    def todict(self):
        if self:
            job_type = self._job_type
            method = self._method
            charge = self._charge
            multiplicity = self._multiplicity
            options = self._options
            atoms_dict = {key:val.todict() for key,val in self._atoms_dict.items()}
            dct = {"job_type":job_type,
                "method":method,
                "charge":charge,
                "multiplicity":multiplicity,
                "options":options,
                "atoms_dict":atoms_dict}
            return dct
        return {}