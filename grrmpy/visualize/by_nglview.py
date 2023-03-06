import numpy as np
import nglview as nv
from ase import Atoms

def view_with_index(atoms , label_color="black", label_scale=1.0):
    if type(atoms) != Atoms:
        raise TypeError("ase.Atomを引数にとって下さい")
    v = nv.show_ase(atoms, viewer="ngl")
    v.add_label(
        color=label_color, labelType="text",
        labelText=[atoms[i].symbol + str(i) for i in range(atoms.get_global_number_of_atoms())],
        zOffset=1.0, attachment='middle_center', radius=label_scale
    )
    return v

def _get_standard_pos(atoms: Atoms) -> np.ndarray:
    # NGLViewer shows atoms in standard positions.
    pos = atoms.positions
    if atoms.get_pbc().any():
        rcell, rot_t = atoms.cell.standard_form()
        standard_pos = pos.dot(rot_t.T)
    else:
        standard_pos = pos
    return standard_pos

def view_with_coordinate(atoms: Atoms, radius=0.5):
    if type(atoms) != Atoms:
        raise TypeError("ase.Atomを引数にとって下さい")
    darkgrey = [0.6, 0.6, 0.6]
    v = nv.show_ase(atoms, viewer="ngl")
    pos = atoms.positions
    standard_pos = _get_standard_pos(atoms)
    for i in range(len(atoms)):
        v.shape.add_sphere(standard_pos[i].tolist(), darkgrey, radius, f"x:{pos[i, 0]}, y:{pos[i, 1]}, z:{pos[i, 2]}")
    return v
