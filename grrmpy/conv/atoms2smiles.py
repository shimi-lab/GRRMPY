from rdkit import Chem
from grrmpy.conv.atoms2mol import atoms2mol,atomslist2mols

def atoms2smiles(atoms,target0:list=None,target1:list=[],target2:list=[],mult:float=1,**kwargs) -> str:
    mol = atoms2mol(atoms,target0,target1,target2,mult,**kwargs)
    return Chem.MolToSmiles(mol)

def atomslist2smileses(atoms_list,target0=None,target1=[],target2=[],mult=1,**kwargs):
    mols = atomslist2mols(atoms_list,target0,target1,target2,mult,**kwargs)
    return [Chem.MolToSmiles(mol) if mol else None for mol in mols]