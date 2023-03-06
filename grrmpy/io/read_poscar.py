from ase.io import read

def get_cell(poscar_file):
    atoms = read(poscar_file,format="vasp")
    cell = atoms.get_cell()
    return cell

def get_pbc(poscar_file):
    atoms = read(poscar_file,format="vasp")
    pbc = atoms.get_pbc()
    return pbc
    
def get_cell_and_pbc(poscar_file):
    atoms = read(poscar_file,format="vasp")
    cell = atoms.get_cell()
    pbc = atoms.get_pbc()
    return cell, pbc