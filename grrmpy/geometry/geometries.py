import pandas as pd
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem
from rdkit.Avalon import pyAvalonTools
from copy import deepcopy
from grrmpy.geometry import Geometry
from ase import Atoms

"""
Geometriesは足し算できない
groupとclusterは呼び出しに時間をかけたくないのでself.__group,self.__clusterに格納する
"""

class Geometries():
    def __init__(self,atoms_list=None,target0=None,target1=[],target2=[],mult=1.0,**kwargs):  
        self._geometry_list = None # Geometryを格納する
        self.__atoms_list = atoms_list
        self.__target0 = target0
        self.__target1 = target1
        self.__target2 = target2
        self.__mult = mult
        if not atoms_list is None:
            self._geometry_list = self.build_geometry_list(self.atoms_list,self.target0,self.target1,self.target2,self.mult)
            self.__group,self.__cluster = self._assign_group_and_cluster(self.smileses)
        
    @property
    def atoms_list(self):
        return self.__atoms_list
    
    def set_atoms_list(self,atoms_list):
        if type(atoms_list) == list:
            self.__atoms_list = atoms_list
            self._geometry_list = self.build_geometry_list(self.atoms_list,self.target0,self.target1,self.target2,self.mult)
        else:
            raise TypeError("atoms_listにはAtomsのリストを指定してください")
        
    @property
    def target0(self):
        return self.__target0

    def set_target0(self, target0):
        if target0 is not None and type(target0)!=list:
            raise TypeError("target0はリストで指定して下さい.")
        else:
            self.__target0 = target0
            if target0 is None: 
                self.__target1 = []
                self.__target2 = []
            self._geometry_list = self.build_geometry_list(self.atoms_list,self.target0,self.target1,self.target2,self.mult)
    
    @property
    def target1(self):
        return self.__target1
    
    def set_target1(self, target1):
        if type(target1) == list:
            self.__target1 = target1
            self._geometry_list = self.build_geometry_list(self.atoms_list,self.target0,self.target1,self.target2,self.mult)
        else:
            raise TypeError("target1はリストで指定して下さい")
    
    @property
    def target2(self):
        return self.__target2
    
    def set_target2(self, target2):
        if type(target2) == list:
            self.__target1 = target2
            self._geometry_list = self.build_geometry_list(self.atoms_list,self.target0,self.target1,self.target2,self.mult)
        else:
            raise TypeError("target2はリストで指定して下さい")
    
    @property
    def mult(self):
        return self.__mult
    
    def set_mult(self,mult):
        if type(mult) == float:
            self.__mult = mult
            self._geometry_list = self.build_geometry_list(self.atoms_list,self.target0,self.target1,self.target2,self.mult)
        else:
            raise TypeError("multはfloatで指定して下さい")
        
    @property
    def mols(self):
        if  self._geometry_list is not None:
            return [geometry.mol for geometry in self._geometry_list]
        else:
            raise RuntimeError(f"{__class__.__name__}のatoms_listがNoneであるためmolsを呼び出せません")
        
    @property
    def smileses(self):
        if  self._geometry_list is not None:
            return [geometry.smiles for geometry in self._geometry_list]
        else:
            raise RuntimeError(f"{__class__.__name__}のatoms_listがNoneであるためsmilesesを呼び出せません")
        
    @property
    def group(self):
        if  self._geometry_list is not None:
            return self.__group
        else:
            raise RuntimeError(f"{__class__.__name__}のatoms_listがNoneであるためgroupを呼び出せません")
        
    def set_group(self,group):
        if  self._geometry_list is None:
            raise RuntimeError(f"{__class__.__name__}のatoms_listがNoneであるため実行できません")
        elif type(group) == list:
            if len(group) == len(self.atoms_list):
                pass
            else:
                raise Exception("要素数が異なります")
        else:
            raise TypeError("型が異なります: groupはlistです")
        
    @property
    def cluster(self):
        if  self._geometry_list is not None:
            return self.__cluster
        else:
            raise RuntimeError(f"{__class__.__name__}のatoms_listがNoneであるためclusterを呼び出せません")
    
    def build_geometry_list(self,atoms_list,target0,target1,target2,mult,**kwargs):
        return [Geometry(atoms,target0,target1,target2,mult,**kwargs) for atoms in atoms_list]
    
    def _assign_group_and_cluster(self,smileses):
        df = pd.DataFrame({'similes': smileses})
        unique = df['similes'].unique().tolist()
        group = [unique.index(smile) for smile in smileses]
        culster = self.group2cluster(group)
        return group, culster
    
    def group2cluster(self,group):
        geoup_n = max(group) + 1 # groupの数
        culster = [[] for _ in range(geoup_n)]
        for i,group_name in enumerate(group):
            culster[group_name].append(i)
        return culster
    
    def same(self,index,other=None):
        if other is None:
            target_name = self.group[index]
            return [i for i,name in enumerate(self.group) if target_name==name]
        else:
            if type(other) != Geometries:
                raise TypeError("Geometry同士でないと比較できません")
            smiles = self.smileses[index]
            other_smileses = other.smileses
            return [i for i,other_smiles in enumerate(other_smileses) if other_smiles==smiles]
        
    def similar(self,index,method="maccs"):
        mols = deepcopy(self.mols)
        valid_mol_dict = {i:mol for i,mol in enumerate(mols) if mol}
        valid_mol = valid_mol_dict.values()
        if not index in valid_mol_dict.keys():
            raise TypeError("指定されたindexの構造が存在しません")
        for mol in valid_mol: 
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS|\
                Chem.SanitizeFlags.SANITIZE_KEKULIZE|\
                Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|\
                Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|\
                Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|\
                Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                catchErrors=True) #長いので改行
        valid_index = valid_mol_dict.keys()
        real_index = dict(zip(valid_index,[i for i in range(len(valid_index))]))
        if method == "maccs":
            similar_list = self._macss(real_index[index],valid_mol)
        elif method == "tanimoto":
            similar_list = self._tanimoto(real_index[index],valid_mol)
        elif method == "avalon":
            similar_list = self._avalon(real_index[index],valid_mol)
        similar_dict = dict(zip(valid_mol_dict.keys(),similar_list))
        similar_dict = sorted(similar_dict.items(),key=lambda x: x[1], reverse=True)
        return similar_dict
    
    def _macss(self,index,mols):
        maccs_fps = [AllChem.GetMACCSKeysFingerprint(mol) for mol in mols]
        maccs = DataStructs.BulkTanimotoSimilarity(maccs_fps[index], maccs_fps)
        return maccs
    
    def _tanimoto(self,index,mols):
        morgan_fp = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
        tanimoto = DataStructs.BulkTanimotoSimilarity(morgan_fp[index], morgan_fp)
        return tanimoto
    
    def _avalon(self,index,mols):
        avalon_fps = [pyAvalonTools.GetAvalonFP(mol) for mol in mols]
        avalon = DataStructs.BulkTanimotoSimilarity(avalon_fps[index], avalon_fps)
        return avalon
    
    def __len__(self):
        return len(self.atoms_list)
    
    def __iter__(self):
        return iter(self._geometry_list)
    
    def __getitem__(self,index):
        if type(index) == slice:
            geometry_list = self._geometry_list[index]
            return self._build_from_geometry(geometry_list)
        index = np.array(index)
        if index.ndim ==0 and index.dtype == int:
            return self._geometry_list[index]
        elif index.ndim ==1 and index.dtype == int:
            geometry_list = [self._geometry_list[i] for i in index]
            return self._build_from_geometry(geometry_list)
        elif index.ndim ==1 and index.dtype == bool:
            if len(self) == len(index):
                geometry_list = [self._geometry_list[i] for i,b in enumerate(index) if b]
                return self._build_from_geometry(geometry_list)
            else:
                raise IndexError("ブーリアンインデックスの要素数が一致しません")
        else:
            raise IndexError("Indexに指定できるのはint,list,sliceのいずれかです.")
            
    def _build_from_geometry(self,geometry_list):
        """copyではないので注意"""
        ref = geometry_list[0]
        if all([ref._has_same_parameters_as(geo) for geo in geometry_list]):
            """全ての要素のパラメータが一致していたら"""
            atoms_list = [geo.atoms for geo in geometry_list]
            t0 = ref.target0
            t1 = ref.target1
            t2 = ref.target2
            mult = ref.mult
            return self.__class__(atoms_list,t0,t1,t2,mult)
            
    def _has_same_parameters_as(self,other):
        if type(other) != Geometry and type(other) != Geometries:
            raise TypeError("GeometryまたはGeometriesでないと比較できません")
        else:
            t0 = self.target0 == other.target0
            t1 = self.target1 == other.target1
            t2 = self.target2 == other.target2
            mult = self.mult == other.mult
            return all(t0,t1,t2,mult)
    
    def __eq__(self, other):
        smls1 = self.smileses
        smls2 = other.smileses
        if len(smls1) == len(smls2):
            return np.array([sml1 == sml2 for sml1,sml2 in zip(smls1,smls2)])
        else:
            raise Exception("要素数が一致しないため比較できません")
    
    def __ne__(self, other):
        result = self == other
        return np.array([not b for b in result])
    
    def __str__(self):
        name = self.__class__.__name__
        t0 = self.target0
        t1 = self.target1
        t2 = self.target2
        mult = self.mult
        return f"{name}(target0={t0},target1={t1},target2={t2},mult={mult})"
    
    def __repr__(self) -> str:
        name = self.__class__.__name__
        atoms_list = self.atoms_list
        t0 = self.target0
        t1 = self.target1
        t2 = self.target2
        mult = self.mult
        return f"{name}(atoms={atoms_list},target0={t0},target1={t1},target2={t2},mult={mult})"
    
    @classmethod
    def fromdit(cls,dct):
        atoms_list = [Atoms.fromdict(atoms) for atoms in dct["atoms_list"]] if dct["atoms_list"] else None        
        target0 = dct["target0"]
        target1 = dct["target1"]
        target2 = dct["target2"]
        mult = dct["mult"]
        return cls(atoms_list,target0,target1,target2,mult)
    
    def todict(self):
        geometry_list = [geometry.todict() for geometry in self._geometry_list] if self._geometry_list else None
        atoms_list = [atoms.todict() for atoms in self.atoms_list] if self.atoms_list else None
        target0 = self.target0
        target1 = self.target1
        target2 = self.target2
        mult = self.mult
        dct = {"atoms_list":atoms_list,
               "target0":target0,
               "target1":target1,
               "target2":target2,
               "mult":mult,
               "geometry_list":geometry_list
               }
        return dct