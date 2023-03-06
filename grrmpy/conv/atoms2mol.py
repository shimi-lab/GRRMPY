from ase.geometry.analysis import Analysis
from ase.neighborlist import build_neighbor_list,natural_cutoffs
from rdkit import Chem

def atoms2mol(atoms,target0=None,target1=[],target2=[],mult=1.0,**kwargs):
    """_summary_

    Parameters:
    
    atoms: Atoms
        Atomsオブジェクト
    target0: list of integers
        | target0とtarget0+target1+target2に指定された原子間の結合を見る.
        | Noneの場合全ての原子間の結合を見る
    target1: list of integers
        target1とtarget0+target1に指定した原子間の結合を見る
    target2: list of integers
        target2とtarget0に指定した原子間の結合を見る
    mult: float
        数字が大きい程,遠く離れていても結合していると判断する.
    kwargs:
        | 各元素の共有結合半径を設定する
        | ex) Al=1.2

    Returns:
        Mol: RDkidのMolオブジェクト
    """
    target0,target1,_,target3 = _renumbering_target(len(atoms),target0,target1,target2)
    return _atoms2mol(atoms,target0,target1,target3,mult,**kwargs)

def atomslist2mols(atoms_list,target0=None,target1=[],target2=[],mult=1,**kwargs):
    for atoms in atoms_list:
        if atoms:
            n_atoms = len(atoms)
            break
    target0,target1,_,target3 = _renumbering_target(n_atoms,target0,target1,target2)
    mols = [_atoms2mol(atoms,target0,target1,target3,mult,**kwargs) if atoms else None for atoms in atoms_list]
    return mols

def _renumbering_target(n_atoms,target0,target1,target2):
    if target0 is None:
        target0 = [i for i in range(n_atoms)]
        target1,target2,target3 = [],[],[]
        return target0,target1,target2,target3
    else:
        target3 = [i for i in range(n_atoms) if not i in target0+target1+target2] # どれとも結合しない原子
        new_target0,new_target1,new_target2 = [],[],[]
        n = 0
        for i in range(n_atoms):
            if not i in target3:
                if i in target0:
                    new_target0.append(n)
                elif i in target1:
                    new_target1.append(n)
                elif i in target2:
                    new_target2.append(n)
                n += 1
        return new_target0,new_target1,new_target2,target3
    
def _atoms2mol(atoms,target0,target1,target3,mult=1,**kwargs):
    atoms = atoms.copy()
    del atoms[target3]
    m = Chem.MolFromSmiles('')
    mol = Chem.RWMol(m) #空のmolオブジェクトを作成
    for sybmols in atoms.get_chemical_symbols():
        """原子を追加するこの時点では結合はない"""
        mol.AddAtom(Chem.Atom(sybmols))
    cutoff = natural_cutoffs(atoms, mult=mult,**kwargs) # cutoff半径(結合の許容度)を設定，multは係数，maltが大きいと結合しているとみなされる
    nl = build_neighbor_list(atoms, cutoff)
    ana = Analysis(atoms, nl=nl)
    bonds = ana.unique_bonds
    for a_idx,bond in enumerate(bonds[0]):
        """結合を追加する"""
        for b_idx in bond:
            if a_idx in target0 or b_idx in target0:
                mol.AddBond(int(a_idx),int(b_idx),Chem.BondType.SINGLE)
            elif a_idx in target1 and b_idx in target1:
                mol.AddBond(int(a_idx),int(b_idx),Chem.BondType.SINGLE)
    for atom in mol.GetAtoms(): # 勝手に水素を表示しないようにする
        atom.SetProp("atomLabel", atom.GetSymbol())
    return mol