from audioop import mul
from types import new_class
from ase import Atoms
from rdkit import Chem
import copy
from grrmpy.conv.atoms2mol import atoms2mol
import grrmpy.geometry as gg


class Geometry():
    def __init__(self,atoms:Atoms=None,target0=None,target1=[],target2=[],mult=1.0,**kwargs):
        self.__smiles = None
        self.__mol = None
        self.__target0 = target0
        self.__target1 = target1
        self.__target2 = target2
        self.__mult = mult
        self.atoms = atoms # __atomsにしないsetterで設定するため(smileとmolを同時生成している)!!
        
    @property
    def target0(self):
        return self.__target0
    def set_target0(self,target0:list):
        self.__target0 = target0
        self.build_mols_and_smiles()
    
    @property
    def target1(self):
        return self.__target1
    def set_target1(self,target1:list):
        self.__target1 = target1
        self.build_mols_and_smiles()
    
    @property
    def target2(self):
        return self.__target2
    def set_target2(self,target2:list):
        self.__target2 = target2
        self.build_mols_and_smiles()
        
    @property
    def mult(self):
        return self.__mult
    def set_mult(self,mult:float):
        self.__mult = mult
        self.build_mols_and_smiles()
        
    @property
    def atoms(self):
        return self.__atoms
    @atoms.setter
    def atoms(self,atoms:Atoms):
        if atoms is None:
             self.__atoms = atoms
        elif type(atoms) == Atoms:
            self.__atoms = atoms
            self.build_mols_and_smiles()
            
    @property
    def smiles(self):
        return self.__smiles
    
    @property
    def mol(self):
        return self.__mol
        
    def build_mols_and_smiles(self):
        if self.atoms is not None:
            self.__mol = atoms2mol(self.atoms,self.target0,self.target1,self.target2,self.mult)
            self.__smiles = Chem.MolToSmiles(self.mol)
        
    def _has_same_parameters_as(self,other):
        if type(other) != Geometry and type(other) != gg.Geometries:
            raise TypeError("GeometryまたはGeometriesでないと比較できません")
        else:
            t0 = self.target0 == other.target0
            t1 = self.target1 == other.target1
            t2 = self.target2 == other.target2
            mult = self.mult == other.mult
            return all(t0,t1,t2,mult)

    def __eq__(self, other: object) -> bool:
        if type(other) == Geometry:
            return self.smiles == other.smiles
        else:
            raise RuntimeError("Geometry同士でないと比較できません")
        
    def __ne__(self, other: object) -> bool:
        return not self == other
    
    def __str__(self):
        if self.atoms is not None:
            text = f"{__class__.__name__}(smiles={self.smiles},atoms={self.atoms.__repr__()})"
            return text
        else:
            return f"{__class__.__name__}()"
        
    def __repr__(self):
        name = self.__class__.__name__     
        atoms = self.atoms
        t0 = self.target0
        t1 = self.target1
        t2 = self.target2
        mult = self.mult
        return f"{name}(atoms={atoms},target0={t0},target1={t1},target2={t2},mult={mult})"
        
    def __copy__(self):
        atoms = self.atoms.copy() if self.atoms else None
        target0 = copy.copy(self.target0) if self.target0 else None
        target1 = copy.copy(self.target1) if self.target1 != [] else []
        target2 = copy.copy(self.target2) if self.target2 != [] else []
        mult = self.mult
        return self.__class__(atoms,target0,target1,target2,mult)
    
    def copy(self):
        return self.__copy__()
    
    @classmethod
    def fromdict(cls,dct):
        atoms = dct["atoms"]
        if atoms is not None:
            atoms = Atoms.fromdict(atoms)
        t0 = dct["target0"]
        t1 = dct["target1"]
        t2 = dct["target2"]
        mult = dct["mult"]
        return cls(atoms,t0,t1,t2,mult)
    
    def todict(self):
        atoms = self.atoms.todict() if self.atoms else None
        target0 = self.target0
        target1 = self.target1
        target2 = self.target2
        mult = self.mult
        dct = {"atoms":atoms,"target0":target0,"target1":target1,"target2":target2,"mult":mult}
        return dct